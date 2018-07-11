/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  Clifford Wolf <clifford@symbioticeda.com>
 *  Copyright (C) 2018  David Shah <david@symbioticeda.com>
 *
 *  Simulated annealing implementation based on arachne-pnr
 *  Copyright (C) 2015-2018 Cotton Seed
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include "place_vpr.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <queue>
#include <set>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "log.h"
#include "place_common.h"
#include "place_legaliser.h"
#include "timing.h"
#include "util.h"
NEXTPNR_NAMESPACE_BEGIN

static double get_std_dev(int n, double sum_x_squared, double av_x) {

	/* Returns the standard deviation of data set x.  There are n sample points, *
	 * sum_x_squared is the summation over n of x^2 and av_x is the average x.   *
	 * All operations are done in double precision, since round off error can be *
	 * a problem in the initial temp. std_dev calculation for big circuits.      */

	double std_dev;

	if (n <= 1)
		std_dev = 0.;
	else
		std_dev = (sum_x_squared - n * av_x * av_x) / (double) (n - 1);

	if (std_dev > 0.) /* Very small variances sometimes round negative */
		std_dev = sqrt(std_dev);
	else
		std_dev = 0.;

	return (std_dev);
}

class VPRPlacer
{
  private:
    std::vector<std::vector<BelId>> free_locations;
    std::vector<std::vector<BelId>::const_iterator> free_locations_back;

  public:
    VPRPlacer(Context *ctx) : ctx(ctx)
    {
        int num_bel_types = 0;
        for (auto bel : ctx->getBels()) {
            int x, y;
            bool gb;
            ctx->estimatePosition(bel, x, y, gb);
            BelType type = ctx->getBelType(bel);
            int type_idx;
            if (bel_types.find(type) == bel_types.end()) {
                type_idx = num_bel_types++;
                bel_types[type] = type_idx;
            } else {
                type_idx = bel_types.at(type);
            }
            if (int(fast_bels.size()) < type_idx + 1) {
                fast_bels.resize(type_idx + 1);
                free_locations.resize(type_idx + 1);
            }
            if (int(fast_bels.at(type_idx).size()) < (x + 1))
                fast_bels.at(type_idx).resize(x + 1);
            if (int(fast_bels.at(type_idx).at(x).size()) < (y + 1))
                fast_bels.at(type_idx).at(x).resize(y + 1);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
            fast_bels.at(type_idx).at(x).at(y).push_back(bel);
            free_locations[type_idx].push_back(bel);
        }
        diameter = std::max(max_x, max_y) + 1;
    }

    bool place()
    {
        log_break();

        size_t placed_cells = 0;
        // Initial constraints placer
        for (auto &cell_entry : ctx->cells) {
            CellInfo *cell = cell_entry.second.get();
            auto loc = cell->attrs.find(ctx->id("BEL"));
            if (loc != cell->attrs.end()) {
                std::string loc_name = loc->second;
                BelId bel = ctx->getBelByName(ctx->id(loc_name));
                if (bel == BelId()) {
                    log_error("No Bel named \'%s\' located for "
                              "this chip (processing BEL attribute on \'%s\')\n",
                              loc_name.c_str(), cell->name.c_str(ctx));
                }

                BelType bel_type = ctx->getBelType(bel);
                if (bel_type != ctx->belTypeFromId(cell->type)) {
                    log_error("Bel \'%s\' of type \'%s\' does not match cell "
                              "\'%s\' of type \'%s\'",
                              loc_name.c_str(), ctx->belTypeToId(bel_type).c_str(ctx), cell->name.c_str(ctx),
                              cell->type.c_str(ctx));
                }

                ctx->bindBel(bel, cell->name, STRENGTH_USER);
                locked_bels.insert(bel);
                placed_cells++;
            }
        }
        int constr_placed_cells = placed_cells;
        log_info("Placed %d cells based on constraints.\n", int(placed_cells));

        // Sort to-place cells for deterministic initial placement
        for (auto &cell : ctx->cells) {
            CellInfo *ci = cell.second.get();
            if (ci->bel == BelId()) {
                autoplaced.push_back(cell.second.get());
            }
        }
        std::sort(autoplaced.begin(), autoplaced.end(), [](CellInfo *a, CellInfo *b) { return a->name < b->name; });
        ctx->shuffle(autoplaced);

        // Remove locked_bels from free_locations
        // TODO Make this more efficient
        for (auto& i : free_locations)
            for (auto j = i.begin(); j != i.end(); ) {
                if (locked_bels.count(*j))
                    j = i.erase(j);
                else
                    ++j;
            }

        // Place cells randomly initially
        log_info("Creating initial placement for remaining %d cells.\n", int(autoplaced.size()));

        vpr_initial_placement(placed_cells, constr_placed_cells);

        log_info("Running simulated annealing placer.\n");

        // Calculate wirelength after initial placement
        curr_wirelength = 0;
        curr_tns = 0;
        for (auto &net : ctx->nets) {
            wirelen_t wl = get_net_wirelength(ctx, net.second.get(), curr_tns);
            wirelengths[net.first] = wl;
            curr_wirelength += wl;
        }

    	num_swap_rejected = 0;
    	num_swap_accepted = 0;
//    	num_swap_aborted = 0;
    
        move_lim = int(inner_num * pow(autoplaced.size(), 1.3333));
    
    	/* Sometimes I want to run the router with a random placement.  Avoid *
    	 * using 0 moves to stop division by 0 and 0 length vector problems,  *
    	 * by setting move_lim to 1 (which is still too small to do any       *
    	 * significant optimization).                                         */
    	if (move_lim <= 0)
    		move_lim = 1;
    
    	rlim = float(std::max(max_x, max_y));
    
    	//first_rlim = rlim; /*used in timing-driven placement for exponent computation */
    	//final_rlim = 1;
    	//inverse_delta_rlim = 1 / (first_rlim - final_rlim);
    
    	t = vpr_starting_t(move_lim /*, rlim,
    			placer_opts.place_algorithm, placer_opts.timing_tradeoff,
    			inverse_prev_bb_cost, inverse_prev_timing_cost, &delay_cost*/);

#if 1
        int n_no_progress = 0;
        double avg_wirelength = curr_wirelength;
//        t = 10000;

        // Main simulated annealing loop
        for (int iter = 1;; iter++) {
            n_move = n_accept = 0;
            improved = false;

            if (iter % 5 == 0 || iter == 1)
                log_info("  at iteration #%d: temp = %f, wire length = "
                         "%.0f, est tns = %.02fns\n",
                         iter, t, double(curr_wirelength), curr_tns);

            for (int m = 0; m < 15; ++m) {
                // Loop through all automatically placed cells
                for (auto cell : autoplaced) {
                    // Find another random Bel for this cell
                    BelId try_bel = random_bel_for_cell(cell);
                    // If valid, try and swap to a new position and see if
                    // the new position is valid/worthwhile
                    if (try_bel != BelId() && try_bel != cell->bel)
                        try_swap_position(cell, try_bel);
                }
            }
            // Heuristic to improve placement on the 8k
            if (improved)
                n_no_progress = 0;
            else
                n_no_progress++;

            if (t <= 1e-3 && n_no_progress >= 5) {
                if (iter % 5 != 0)
                    log_info("  at iteration #%d: temp = %f, wire length = %f\n", iter, t, double(curr_wirelength));
                break;
            }

            double Raccept = double(n_accept) / double(n_move);

            int M = std::max(max_x, max_y) + 1;

            double upper = 0.6, lower = 0.4;

            if (curr_wirelength < 0.95 * avg_wirelength) {
                avg_wirelength = 0.8 * avg_wirelength + 0.2 * curr_wirelength;
            } else {
                if (Raccept >= 0.8) {
                    t *= 0.7;
                } else if (Raccept > upper) {
                    if (diameter < M)
                        diameter++;
                    else
                        t *= 0.9;
                } else if (Raccept > lower) {
                    t *= 0.95;
                } else {
                    // Raccept < 0.3
                    if (diameter > 1)
                        diameter--;
                    else
                        t *= 0.8;
                }
            }
            // Once cooled below legalise threshold, run legalisation and start requiring
            // legal moves only
            if (t < legalise_temp && !require_legal) {
                legalise_design(ctx);
                require_legal = true;
                autoplaced.clear();
                for (auto cell : sorted(ctx->cells)) {
                    if (cell.second->belStrength < STRENGTH_STRONG)
                        autoplaced.push_back(cell.second);
                }
                t = post_legalise_temp;
                diameter *= post_legalise_dia_scale;
                ctx->shuffle(autoplaced);
                assign_budget(ctx);
            }

            // Recalculate total wirelength entirely to avoid rounding errors
            // accumulating over time
            curr_wirelength = 0;
            curr_tns = 0;
            for (auto &net : ctx->nets) {
                wirelen_t wl = get_net_wirelength(ctx, net.second.get(), curr_tns);
                wirelengths[net.first] = wl;
                curr_wirelength += wl;
            }
        }
#else
	tot_iter = 0;
	//moves_since_cost_recompute = 0;

	/* Outer loop of the simmulated annealing begins */
	while (!vpr_exit_crit(temp, wirelength)) {

//		if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//			cost = 1;
//		}

//		outer_loop_recompute_criticalities(placer_opts, num_connections,
//			crit_exponent, bb_cost, &place_delay_value, &timing_cost, &delay_cost,
//			&outer_crit_iter_count, &inverse_prev_timing_cost, &inverse_prev_bb_cost,
//            netlist_pin_lookup,
//#ifdef ENABLE_CLASSIC_VPR_STA
//            slacks,
//            timing_inf,
//#endif
//            *timing_info);

		placement_inner_loop(t, rlim, placer_opts, inverse_prev_bb_cost, inverse_prev_timing_cost,
			move_lim, crit_exponent, inner_recompute_limit, &stats,
			&cost, &bb_cost, &timing_cost, &delay_cost,
#ifdef ENABLE_CLASSIC_VPR_STA
            slacks,
            timing_inf,
#endif
            netlist_pin_lookup,
            *timing_info);

//		/* Lines below prevent too much round-off error from accumulating *
//		 * in the cost over many iterations.  This round-off can lead to  *
//		 * error checks failing because the cost is different from what   *
//		 * you get when you recompute from scratch.                       */
//
//		moves_since_cost_recompute += move_lim;
//		if (moves_since_cost_recompute > MAX_MOVES_BEFORE_RECOMPUTE) {
//			new_bb_cost = recompute_bb_cost();
//			if (fabs(new_bb_cost - bb_cost) > bb_cost * ERROR_TOL) {
//				vpr_throw(VPR_ERROR_PLACE, __FILE__, __LINE__,
//						"in try_place: new_bb_cost = %g, old bb_cost = %g\n",
//						new_bb_cost, bb_cost);
//			}
//			bb_cost = new_bb_cost;
//
//			if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//				comp_td_costs(&new_timing_cost, &new_delay_cost);
//				if (fabs(new_timing_cost - timing_cost) > timing_cost * ERROR_TOL) {
//					vpr_throw(VPR_ERROR_PLACE, __FILE__, __LINE__,
//							"in try_place: new_timing_cost = %g, old timing_cost = %g, ERROR_TOL = %g\n",
//							new_timing_cost, timing_cost, ERROR_TOL);
//				}
//				if (fabs(new_delay_cost - delay_cost) > delay_cost * ERROR_TOL) {
//					vpr_throw(VPR_ERROR_PLACE, __FILE__, __LINE__,
//							"in try_place: new_delay_cost = %g, old delay_cost = %g, ERROR_TOL = %g\n",
//							new_delay_cost, delay_cost, ERROR_TOL);
//				}
//				timing_cost = new_timing_cost;
//			}
//
//			if (placer_opts.place_algorithm == BOUNDING_BOX_PLACE) {
//				cost = new_bb_cost;
//			}
//			moves_since_cost_recompute = 0;
//		}
//
		tot_iter += move_lim;
		success_rat = ((float) stats.success_sum) / move_lim;
//		if (stats.success_sum == 0) {
//			stats.av_cost = cost;
//			stats.av_bb_cost = bb_cost;
//			stats.av_timing_cost = timing_cost;
//			stats.av_delay_cost = delay_cost;
//		} else {
//			stats.av_cost /= stats.success_sum;
//			stats.av_bb_cost /= stats.success_sum;
//			stats.av_timing_cost /= stats.success_sum;
//			stats.av_delay_cost /= stats.success_sum;
//		}
//		std_dev = get_std_dev(stats.success_sum, stats.sum_of_squares, stats.av_cost);

//		oldt = t; /* for finding and printing alpha. */
		vpr_update_t();

//        if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//            critical_path = timing_info->least_slack_critical_path();
//            sTNS = timing_info->setup_total_negative_slack();
//            sWNS = timing_info->setup_worst_negative_slack();
//        }
//
//        vtr::printf_info("%7.3f "
//                         "%7.4f %10.4f %-10.5g %-10.5g "
//                         "%-10.5g %7.3f % 10.3g % 8.3f "
//                         "%7.4f %7.4f %7.4f %6.3f"
//                         "%9d %6.3f\n",
//                         oldt,
//                         stats.av_cost, stats.av_bb_cost, stats.av_timing_cost, stats.av_delay_cost,
//                         place_delay_value, 1e9*critical_path.delay(), 1e9*sTNS, 1e9*sWNS,
//                         success_rat, std_dev, rlim, crit_exponent,
//                         tot_iter, t / oldt);
//
//#ifdef ENABLE_CLASSIC_VPR_STA
//        if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//            float cpd_diff_ns = std::abs(get_critical_path_delay() - 1e9*critical_path.delay());
//            if(cpd_diff_ns > ERROR_TOL) {
//                print_classic_cpds();
//                print_tatum_cpds(timing_info->critical_paths());
//
//                vpr_throw(VPR_ERROR_TIMING, __FILE__, __LINE__, "Classic VPR and Tatum critical paths do not match (%g and %g respectively)", get_critical_path_delay(), 1e9*critical_path.delay());
//            }
//        }
//#endif

		sprintf(msg, "Cost: %g  BB Cost %g  TD Cost %g  Temperature: %g",
				cost, bb_cost, timing_cost, t);
//		update_screen(ScreenUpdatePriority::MINOR, msg, PLACEMENT, timing_info);
		update_rlim(&rlim, success_rat, device_ctx.grid);

//		if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//			crit_exponent = (1 - (rlim - final_rlim) * inverse_delta_rlim)
//					* (placer_opts.td_place_exp_last - placer_opts.td_place_exp_first)
//					+ placer_opts.td_place_exp_first;
//		}
//#ifdef VERBOSE
//		if (getEchoEnabled()) {
//			print_clb_placement("first_iteration_clb_placement.echo");
//		}
//#endif
	}
	/* Outer loop of the simmulated annealing ends */


	outer_loop_recompute_criticalities(placer_opts, num_connections,
			crit_exponent, bb_cost, &place_delay_value, &timing_cost, &delay_cost,
			&outer_crit_iter_count, &inverse_prev_timing_cost, &inverse_prev_bb_cost,
            netlist_pin_lookup,
#ifdef ENABLE_CLASSIC_VPR_STA
            slacks,
            timing_inf,
#endif
            *timing_info);

	t = 0; /* freeze out */

	/* Run inner loop again with temperature = 0 so as to accept only swaps
	 * which reduce the cost of the placement */
	placement_inner_loop(t, rlim, placer_opts, inverse_prev_bb_cost, inverse_prev_timing_cost,
			move_lim, crit_exponent, inner_recompute_limit, &stats,
			&cost, &bb_cost, &timing_cost, &delay_cost,
#ifdef ENABLE_CLASSIC_VPR_STA
            slacks,
            timing_inf,
#endif
            netlist_pin_lookup,
            *timing_info);
#endif

        // Final post-pacement validitiy check
        for (auto bel : ctx->getBels()) {
            IdString cell = ctx->getBoundBelCell(bel);
            if (!ctx->isBelLocationValid(bel)) {
                std::string cell_text = "no cell";
                if (cell != IdString())
                    cell_text = std::string("cell '") + cell.str(ctx) + "'";
                if (ctx->force) {
                    log_warning("post-placement validity check failed for Bel '%s' "
                                "(%s)\n",
                                ctx->getBelName(bel).c_str(ctx), cell_text.c_str());
                } else {
                    log_error("post-placement validity check failed for Bel '%s' "
                              "(%s)\n",
                              ctx->getBelName(bel).c_str(ctx), cell_text.c_str());
                }
            }
        }
        return true;
    }

  private:
#if 0
    // Initial random placement
    void place_initial(CellInfo *cell)
    {
        bool all_placed = false;
        int iters = 25;
        while (!all_placed) {
            BelId best_bel = BelId();
            uint64_t best_score = std::numeric_limits<uint64_t>::max(),
                     best_ripup_score = std::numeric_limits<uint64_t>::max();
            CellInfo *ripup_target = nullptr;
            BelId ripup_bel = BelId();
            if (cell->bel != BelId()) {
                ctx->unbindBel(cell->bel);
            }
            BelType targetType = ctx->belTypeFromId(cell->type);
            for (auto bel : ctx->getBels()) {
                if (ctx->getBelType(bel) == targetType && (ctx->isValidBelForCell(cell, bel) || !require_legal)) {
                    if (ctx->checkBelAvail(bel)) {
                        uint64_t score = ctx->rng64();
                        if (score <= best_score) {
                            best_score = score;
                            best_bel = bel;
                        }
                    } else {
                        uint64_t score = ctx->rng64();
                        if (score <= best_ripup_score) {
                            best_ripup_score = score;
                            ripup_target = ctx->cells.at(ctx->getBoundBelCell(bel)).get();
                            ripup_bel = bel;
                        }
                    }
                }
            }
            if (best_bel == BelId()) {
                if (iters == 0 || ripup_bel == BelId())
                    log_error("failed to place cell '%s' of type '%s'\n", cell->name.c_str(ctx), cell->type.c_str(ctx));
                --iters;
                ctx->unbindBel(ripup_target->bel);
                best_bel = ripup_bel;
            } else {
                all_placed = true;
            }
            ctx->bindBel(best_bel, cell->name, STRENGTH_WEAK);

            // Back annotate location
            cell->attrs[ctx->id("BEL")] = ctx->getBelName(cell->bel).str(ctx);
            cell = ripup_target;
        }
    }
#endif

    void vpr_initial_placement(size_t& placed_cells, int constr_placed_cells) {
    
    	/* Randomly places the blocks to create an initial placement. We rely on
    	 * the legal_pos array already being loaded.  That legal_pos[itype] is an
    	 * array that gives every legal value of (x,y,z) that can accomodate a block.
    	 * The number of such locations is given by num_legal_pos[itype].
    	 */
//    	int itype, x, y, z, ipos;
//    	int *free_locations; /* [0..device_ctx.num_block_types-1].
//    						  * Stores how many locations there are for this type that *might* still be free.
//    						  * That is, this stores the number of entries in legal_pos[itype] that are worth considering
//    						  * as you look for a free location.
//    						  */
//        auto& device_ctx = g_vpr_ctx.device();
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//        auto& place_ctx = g_vpr_ctx.mutable_placement();

        // free_locations is populated in constructor

//    	initial_placement_pl_macros(MAX_NUM_TRIES_TO_PLACE_MACROS_RANDOMLY, free_locations);
//    
//    	// All the macros are placed, update the legal_pos[][] array
//    	for (itype = 0; itype < device_ctx.num_block_types; itype++) {
//    		VTR_ASSERT(free_locations[itype] >= 0);
//    		for (ipos = 0; ipos < free_locations[itype]; ipos++) {
//    			x = legal_pos[itype][ipos].x;
//    			y = legal_pos[itype][ipos].y;
//    			z = legal_pos[itype][ipos].z;
//    
//    			// Check if that location is occupied.  If it is, remove from legal_pos
//    			if (place_ctx.grid_blocks[x][y].blocks[z] != EMPTY_BLOCK_ID && place_ctx.grid_blocks[x][y].blocks[z] != INVALID_BLOCK_ID) {
//    				legal_pos[itype][ipos] = legal_pos[itype][free_locations[itype] - 1];
//    				free_locations[itype]--;
//    
//    				// After the move, I need to check this particular entry again
//    				ipos--;
//    				continue;
//    			}
//    		}
//    	} // Finish updating the legal_pos[][] and free_locations[] array
    
    	vpr_initial_placement_blocks(placed_cells, constr_placed_cells);
    
//    	if (pad_loc_type == USER) {
//    		read_user_pad_loc(pad_loc_file);
//    	}
//    
//    	/* Restore legal_pos */
//    	load_legal_placements();
//    
//    #ifdef VERBOSE
//    	vtr::printf_info("At end of initial_placement.\n");
//    	if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_INITIAL_CLB_PLACEMENT)) {
//    		print_clb_placement(getEchoFileName(E_ECHO_INITIAL_CLB_PLACEMENT));
//    	}
//    #endif
//    	free(free_locations);
    }

    /* Place blocks that are NOT a part of any macro.
    * We'll randomly place each block in the clustered netlist, one by one. */
    void vpr_initial_placement_blocks(size_t &placed_cells, int& constr_placed_cells) {
//        int itype, ipos, x, y, z;
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//        auto& place_ctx = g_vpr_ctx.mutable_placement();
//        auto& device_ctx = g_vpr_ctx.device();
        
        size_t itype;
        BelId bel;
    
        // Shuffle all free locations once here, rather than picking a block at random
        free_locations_back.reserve(free_locations.size());
        for (auto& free_locations_type : free_locations) {
            ctx->shuffle(free_locations_type);
            free_locations_back.push_back(free_locations_type.cbegin());
        }

        for (auto& cell : autoplaced) {
//    		if (place_ctx.block_locs[blk_id].x != -1) { // -1 is a sentinel for an empty block
//    			// block placed.
//    			continue;
//    		}
    
//    		/* Don't do IOs if the user specifies IOs; we'll read those locations later. */
//    		if (!(is_io_type(cluster_ctx.clb_nlist.block_type(blk_id)) && pad_loc_type == USER)) {
//    
//    		    /* Randomly select a free location of the appropriate type for blk_id.
//    			 * We have a linearized list of all the free locations that can
//    			 * accomodate a block of that type in free_locations[itype].
//    			 * Choose one randomly and put blk_id there. Then we don't want to pick
//    			 * that location again, so remove it from the free_locations array.
//    			 */
//    			itype = cluster_ctx.clb_nlist.block_type(blk_id)->index;
//    			if (free_locations[itype] <= 0) {
//    				vpr_throw(VPR_ERROR_PLACE, __FILE__, __LINE__,
//    						"Initial placement failed.\n"
//    						"Could not place block %s (#%zu); no free locations of type %s (#%d).\n",
//    						cluster_ctx.clb_nlist.block_name(blk_id).c_str(), size_t(blk_id), device_ctx.block_types[itype].name, itype);
//    		    }

    			vpr_initial_placement_location(cell, itype, bel);

                NPNR_ASSERT(ctx->checkBelAvail(bel));
                ctx->bindBel(bel, cell->name, STRENGTH_WEAK);

//                //Mark IOs as fixed if specifying a (fixed) random placement
//                if(is_io_type(cluster_ctx.clb_nlist.block_type(blk_idctx->belTypeFromId(cell->type))) && pad_loc_type == RANDOM) {
//                    place_ctx.block_locs[blk_id].is_fixed = true;
//     			}

                if (ctx->isIO(cell)) {
                    // TODO: Add method to change bind strength without unbinding and re-binding
                    ctx->unbindBel(bel);
                    ctx->bindBel(bel, cell->name, STRENGTH_LOCKED);
                    cell = nullptr;
                }

//    			/* Ensure randomizer doesn't pick this location again, since it's occupied. Could shift all the
//    			* legal positions in legal_pos to remove the entry (choice) we just used, but faster to
//    			* just move the last entry in legal_pos to the spot we just used and decrement the
//    			* count of free_locations. */
//    			legal_pos[itype][ipos] = legal_pos[itype][free_locations[itype] - 1]; /* overwrite used block position */
//    			free_locations[itype]--;

                ++free_locations_back[itype];

//    		}

            ++placed_cells;
            if ((placed_cells - constr_placed_cells) % 500 == 0)
                log_info("  initial placement placed %d/%d cells\n", int(placed_cells - constr_placed_cells),
                         int(autoplaced.size()));
    	}

        // Linear complexity removal of all IOs cells that were set to nullptr
        autoplaced.erase(std::remove(autoplaced.begin(), autoplaced.end(), nullptr), autoplaced.end());

        //if ((placed_cells - constr_placed_cells) % 500 != 0)
            log_info("  initial placement finished with %d unconstrained cells\n", int(autoplaced.size()));

    } 
    void vpr_initial_placement_location(CellInfo *cell, size_t& itype, BelId& bel) {
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//    
        auto type = ctx->belTypeFromId(cell->type);
    	itype = bel_types.at(type);

        bel = *free_locations_back[itype];
    }

    float vpr_starting_t(int max_moves) {
    
    	/* Finds the starting temperature (hot condition).              */
    
    	int i, num_accepted, move_lim;
    	double std_dev, av, sum_of_squares; /* Double important to avoid round off */
    
//    	if (annealing_sched.type == USER_SCHED)
//    		return (annealing_sched.init_t);
//    
//        auto& cluster_ctx = g_vpr_ctx.clustering();
    
    	move_lim = std::min(max_moves, (int) autoplaced.size());
    
    	num_accepted = 0;
    	av = 0.;
    	sum_of_squares = 0.;
    
    	/* Try one move per block.  Set t high so essentially all accepted. */
    
    	for (i = 0; i < move_lim; i++) {
    		auto swap_result = vpr_try_swap(std::numeric_limits<float>::max()/*, cost_ptr, bb_cost_ptr, timing_cost_ptr, rlim,
    				place_algorithm, timing_tradeoff,
    				inverse_prev_bb_cost, inverse_prev_timing_cost, delay_cost_ptr*/);
    
    		if (swap_result/* == ACCEPTED*/) {
    			num_accepted++;
    			av += cost;
    			sum_of_squares += cost * cost;
    			num_swap_accepted++;
    		//} else if (swap_result == ABORTED) {
    		//	num_swap_aborted++;
    		} else {
    			num_swap_rejected++;
    		}
    	}
    
    	if (num_accepted != 0)
    		av /= num_accepted;
    	else
    		av = 0.;
    
    	std_dev = get_std_dev(num_accepted, sum_of_squares, av);
    
    	if (num_accepted != move_lim) {
    		log_warning("Starting t: %d of %d configurations accepted.\n", num_accepted, move_lim);
    	}
    
//    #ifdef VERBOSE
    	log_info("std_dev: %g, average cost: %g, starting temp: %g\n", std_dev, av, 20. * std_dev);
//    #endif
    
    	/* Set the initial temperature to 20 times the standard of deviation */
    	/* so that the initial temperature adjusts according to the circuit */
    	return (20. * std_dev);
    }

    bool vpr_try_swap(float t /*, float *cost, float *bb_cost, float *timing_cost,
    		float rlim,
    		enum e_place_algorithm place_algorithm, float timing_tradeoff,
    		float inverse_prev_bb_cost, float inverse_prev_timing_cost,
    		float *delay_cost*/) {
    
//    	/* Picks some block and moves it to another spot.  If this spot is   *
//    	 * occupied, switch the blocks.  Assess the change in cost function. *
//    	 * rlim is the range limiter.                                        *
//         * Returns whether the swap is accepted, rejected or aborted.        *
//    	 * Passes back the new value of the cost functions.                  */
//    
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//        auto& place_ctx = g_vpr_ctx.mutable_placement();
//    
//    	num_ts_called ++;
    
    	/* I'm using negative values of temp_net_cost as a flag, so DO NOT   *
    	 * use cost functions that can go negative.                          */
    
    	/*float*/ delta_c = 0; /* Change in cost due to this swap. */
    	/*float*/ bb_delta_c = 0;
//    	float timing_delta_c = 0;
//    	float delay_delta_c = 0.0;
    
    	/* Pick a random block to be swapped with another random block.   */
        auto cell_from = vpr_pick_from_block();
        if (!cell_from) {
            return false/*ABORTED*/; //No movable block found
        }
    
//    	int x_from = place_ctx.block_locs[b_from].x;
//    	int y_from = place_ctx.block_locs[b_from].y;
//    	int z_from = place_ctx.block_locs[b_from].z;
//    
//        int x_to = OPEN;
//        int y_to = OPEN;
//        int z_to = OPEN;

        BelId bel_to;
    
//        auto cluster_from_type = cluster_ctx.clb_nlist.block_type(b_from);
//        auto grid_from_type = g_vpr_ctx.device().grid[x_from][y_from].type;
//        VTR_ASSERT(cluster_from_type == grid_from_type);
    
    	if (!vpr_find_to(cell_from, bel_to /*cluster_ctx.clb_nlist.block_type(b_from), rlim, x_from, y_from, &x_to, &y_to, &z_to*/))
    		return false/*REJECTED*/;
    
    #if 0
        auto& grid = g_vpr_ctx.device().grid;
    	int b_to = place_ctx.grid_blocks[x_to][y_to].blocks[z_to];
    	vtr::printf_info( "swap [%d][%d][%d] %s \"%s\" <=> [%d][%d][%d] %s \"%s\"\n",
    		x_from, y_from, z_from, grid[x_from][y_from].type->name, (b_from != -1 ? cluster_ctx.blocks[b_from].name : ""),
    		x_to, y_to, z_to, grid[x_to][y_to].type->name, (b_to != -1 ? cluster_ctx.blocks[b_to].name : ""));
    #endif
    
//    	/* Make the switch in order to make computing the new bounding *
//    	 * box simpler.  If the cost increase is too high, switch them *
//    	 * back.  (place_ctx.block_locs data structures switched, clbs not switched   *
//    	 * until success of move is determined.)                       *
//    	 * Also check that whether those are the only 2 blocks         *
//    	 * to be moved - check for carry chains and other placement    *
//    	 * macros.                                                     */
//    
//    	/* Check whether the from_block is part of a macro first.      *
//    	 * If it is, the whole macro has to be moved. Calculate the    *
//    	 * x, y, z offsets of the swap to maintain relative placements *
//    	 * of the blocks. Abort the swap if the to_block is part of a  *
//    	 * macro (not supported yet).                                  */
//    
//    	bool abort_swap = find_affected_blocks(b_from, x_to, y_to, z_to);
//    
//    	if (abort_swap == false) {

            BelId bel_from = cell_from->bel;
            IdString other = ctx->getBoundBelCell(bel_to);
            CellInfo *cell_to = nullptr;
            if (other != IdString()) {
                cell_to = ctx->cells[other].get();
                if (cell_to->belStrength > STRENGTH_WEAK)
                    return false;
            }
            ctx->unbindBel(bel_from);
            if (other != IdString()) {
                ctx->unbindBel(bel_to);
            }

            for (const auto &port : cell_from->ports)
                if (port.second.net != nullptr)
                    affected_nets.insert(port.second.net);

            if (other != IdString()) {
                for (const auto &port : cell_to->ports)
                    if (port.second.net != nullptr)
                        affected_nets.insert(port.second.net);
            }

    		// Find all the nets affected by this swap and update thier bounding box
    		/*int num_nets_affected =*/ vpr_find_affected_nets_and_update_costs(cell_from, cell_to, bel_from, bel_to, /*place_algorithm,*/ bb_delta_c /*, timing_delta_c, delay_delta_c*/);

    
//    		if (place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//    			/*in this case we redefine delta_c as a combination of timing and bb.  *
//    			 *additionally, we normalize all values, therefore delta_c is in       *
//    			 *relation to 1*/
//    
//    			delta_c = (1 - timing_tradeoff) * bb_delta_c * inverse_prev_bb_cost
//    					+ timing_tradeoff * timing_delta_c * inverse_prev_timing_cost;
//    		} else 
            {
    			delta_c = bb_delta_c;
    		}
    
    		/* 1 -> move accepted, 0 -> rejected. */
    		auto keep_switch = vpr_assess_swap(t);
    
    		if (keep_switch) {
    			cost += delta_c;
    			bb_cost += bb_delta_c;
    
//    			if (place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//    				/*update the point_to_point_timing_cost and point_to_point_delay_cost
//    				 * values from the temporary values */
//    				*timing_cost = *timing_cost + timing_delta_c;
//    				*delay_cost = *delay_cost + delay_delta_c;
//    
//    				update_td_cost();
//    			}
    
//    			/* update net cost functions and reset flags. */
//    			for (int inet_affected = 0; inet_affected < num_nets_affected; inet_affected++) {
//    				ClusterNetId net_id = ts_nets_to_update[inet_affected];
//    
//    				bb_coords[net_id] = ts_bb_coord_new[net_id];
//    				if (cluster_ctx.clb_nlist.net_sinks(net_id).size() >= SMALL_NET)
//    					bb_num_on_edges[net_id] = ts_bb_edge_new[net_id];
//    
//    				net_cost[net_id] = temp_net_cost[net_id];
//    
//    				/* negative temp_net_cost value is acting as a flag. */
//    				temp_net_cost[net_id] =	-1;
//    				bb_updated_before[net_id] = NOT_UPDATED_YET;
//    			}
//    
//    			/* Update clb data structures since we kept the move. */
//    			/* Swap physical location */
//    			for (int iblk = 0; iblk < blocks_affected.num_moved_blocks; iblk++) {
//    
//    				x_to = blocks_affected.moved_blocks[iblk].xnew;
//    				y_to = blocks_affected.moved_blocks[iblk].ynew;
//    				z_to = blocks_affected.moved_blocks[iblk].znew;
//    
//    				x_from = blocks_affected.moved_blocks[iblk].xold;
//    				y_from = blocks_affected.moved_blocks[iblk].yold;
//    				z_from = blocks_affected.moved_blocks[iblk].zold;
//    
//    				b_from = blocks_affected.moved_blocks[iblk].block_num;
//    
//    				place_ctx.grid_blocks[x_to][y_to].blocks[z_to] = b_from;
//    
//    				if (blocks_affected.moved_blocks[iblk].swapped_to_was_empty) {
//    					place_ctx.grid_blocks[x_to][y_to].usage++;
//    				}
//    				if (blocks_affected.moved_blocks[iblk].swapped_from_is_empty) {
//    					place_ctx.grid_blocks[x_from][y_from].usage--;
//    					place_ctx.grid_blocks[x_from][y_from].blocks[z_from] = EMPTY_BLOCK_ID;
//    				}
//
//    			} // Finish updating clb for all blocks

                    for (auto new_wl : new_lengths)
                        wirelengths.at(new_wl.first) = new_wl.second;
    
    
    		} else { /* Move was rejected.  */
    
//    			/* Reset the net cost function flags first. */
//    			for (int inet_affected = 0; inet_affected < num_nets_affected; inet_affected++) {
//    				ClusterNetId net_id = ts_nets_to_update[inet_affected];
//    				temp_net_cost[net_id] = -1;
//    				bb_updated_before[net_id] = NOT_UPDATED_YET;
//    			}
//    
//    			/* Restore the place_ctx.block_locs data structures to their state before the move. */
//    			for (int iblk = 0; iblk < blocks_affected.num_moved_blocks; iblk++) {
//    				b_from = blocks_affected.moved_blocks[iblk].block_num;
//    
//    				place_ctx.block_locs[b_from].x = blocks_affected.moved_blocks[iblk].xold;
//    				place_ctx.block_locs[b_from].y = blocks_affected.moved_blocks[iblk].yold;
//    				place_ctx.block_locs[b_from].z = blocks_affected.moved_blocks[iblk].zold;
//    			}

                if (other != IdString())
                    ctx->unbindBel(bel_from);
                ctx->unbindBel(bel_to);
                ctx->bindBel(bel_from, cell_from->name, STRENGTH_WEAK);
                if (other != IdString())
                    ctx->bindBel(bel_to, other, STRENGTH_WEAK);
                return false;
    		}
    
//    		/* Resets the num_moved_blocks, but do not free blocks_moved array. Defensive Coding */
//    		blocks_affected.num_moved_blocks = 0;
//    
//    #if 0
//            //Check that each accepted swap yields a valid placement
//            check_place(*bb_cost, *timing_cost, place_algorithm, *delay_cost);
//    #endif
    
    		return (keep_switch);
//    	} else {
//    
//    		/* Restore the place_ctx.block_locs data structures to their state before the move. */
//    		for (int iblk = 0; iblk < blocks_affected.num_moved_blocks; iblk++) {
//    			b_from = blocks_affected.moved_blocks[iblk].block_num;
//    
//    			place_ctx.block_locs[b_from].x = blocks_affected.moved_blocks[iblk].xold;
//    			place_ctx.block_locs[b_from].y = blocks_affected.moved_blocks[iblk].yold;
//    			place_ctx.block_locs[b_from].z = blocks_affected.moved_blocks[iblk].zold;
//    		}
//    
//    		/* Resets the num_moved_blocks, but do not free blocks_moved array. Defensive Coding */
//    		blocks_affected.num_moved_blocks = 0;
//    
//    		return ABORTED;
//    	}
    }

    //Pick a random block to be swapped with another random block.
    //If none is found return ClusterBlockId::INVALID()
    CellInfo* vpr_pick_from_block() {
//        /* Some blocks may be fixed, and should never be moved from their *
//         * initial positions. If we randomly selected such a block try    *
//         * another random block.                                          *
//         *                                                                *
//         * We need to track the blocks we have tried to avoid an infinite *
//         * loop if all blocks are fixed.                                  */
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//        auto& place_ctx = g_vpr_ctx.mutable_placement();
//    
//        std::unordered_set<ClusterBlockId> tried_from_blocks;
//    
//        //So long as untried blocks remain
//        while (tried_from_blocks.size() < cluster_ctx.clb_nlist.blocks().size()) {
//    
//            //Pick a block at random
//            ClusterBlockId b_from = ClusterBlockId(vtr::irand((int) cluster_ctx.clb_nlist.blocks().size() - 1));
//    
//            //Record it as tried
//            tried_from_blocks.insert(b_from);
//    
//            if (place_ctx.block_locs[b_from].is_fixed) {
//                continue; //Fixed location, try again
//            }
//    
//            //Found a movable block
//            return b_from;
//        }
//    
//        //No movable blocks found
//        return ClusterBlockId::INVALID();

        // Assume that autoplaced only contains movable blocks
        return autoplaced.at(ctx->rng(int(autoplaced.size())));
    }

    bool vpr_find_to(CellInfo* cell, BelId& bel) {
    
    	/* Returns the point to which I want to swap, properly range limited.
    	 * rlim must always be between 1 and device_ctx.grid.width() - 2 (inclusive) for this routine
    	 * to work. Note -2 for no perim channels
    	 */
    
    	int min_x, max_x, min_y, max_y;
    	int num_tries;
    	int active_area;
    	bool is_legal;
//    	int itype;
//    
//        auto& grid = g_vpr_ctx.device().grid;
//        auto& place_ctx = g_vpr_ctx.placement();
//    
//        auto grid_type = grid[x_from][y_from].type;
//    	VTR_ASSERT(type == grid_type);
    
    	int rlx = std::min<float>(this->max_x, rlim);
    	int rly = std::min<float>(this->max_y, rlim); /* Added rly for aspect_ratio != 1 case. */
    	active_area = 4 * rlx * rly;

        int x_from, y_from;
        bool gb;
        ctx->estimatePosition(cell->bel, x_from, y_from, gb);
    
    	min_x = std::max<float>(0, x_from - rlx);
    	max_x = std::min<float>(this->max_x, x_from + rlx);
    	min_y = std::max<float>(0, y_from - rly);
    	max_y = std::min<float>(this->max_y, y_from + rly);
    
    	if (rlx < 1 || rlx > int(this->max_x)) {
            log_error("in find_to: rlx = %d out of range\n", rlx);
    	}
    	if (rly < 1 || rly > int(this->max_y)) {
    		log_error("in find_to: rly = %d out of range\n", rly);
    	}
    
    	num_tries = 0;
//    	itype = type->index;
        auto type = ctx->belTypeFromId(cell->type);
    	auto itype = bel_types.at(type);
    
        int px_to, py_to;

    	do { /* Until legal */
    		is_legal = true;
    
    		/* Limit the number of tries when searching for an alternative position */
    		if(num_tries >= 2 * std::min<int>(active_area /*/ (type->width * type->height)*/, free_locations[itype].size()) + 10) {
    			/* Tried randomly searching for a suitable position */
    			return false;
    		} else {
    			num_tries++;
    		}
    
    		vpr_find_to_location(cell, bel);
            ctx->estimatePosition(bel, px_to, py_to, gb);

    		if((x_from == px_to) && (y_from == py_to)) {
    			is_legal = false;
    		} else if(px_to > max_x || px_to < min_x || py_to > max_y || py_to < min_y) {
    			is_legal = false;
    		} else if(type != ctx->getBelType(bel)) {
    			is_legal = false;
    		} else {
//    			/* Find z_to and test to validate that the "to" block is *not* fixed */
//    			*pz_to = 0;
//    			if (grid[*px_to][*py_to].type->capacity > 1) {
//    				*pz_to = vtr::irand(grid[*px_to][*py_to].type->capacity - 1);
//    			}
//    			ClusterBlockId b_to = place_ctx.grid_blocks[*px_to][*py_to].blocks[*pz_to];
//    			if ((b_to != EMPTY_BLOCK_ID) && (place_ctx.block_locs[b_to].is_fixed == true)) {
//    				is_legal = false;
//    			}
    		}
    
    		NPNR_ASSERT(px_to >= 0 && px_to <= int(this->max_x));
    		NPNR_ASSERT(py_to >= 0 && py_to <= int(this->max_y));
    	} while (is_legal == false);
    
    	if (px_to < 0 || px_to > int(this->max_x) || py_to < 0 || py_to > int(this->max_y)) {
    		log_error("in routine find_to: (x_to,y_to) = (%d,%d)\n", px_to, py_to);
    	}
    
    	NPNR_ASSERT(type == ctx->getBelType(bel));
    	return true;
    }

    void vpr_find_to_location(CellInfo* cell, BelId& bel/*t_type_ptr type, float rlim,
    		int x_from, int y_from,
    		int *px_to, int *py_to, int *pz_to*/) {
    
//        auto& device_ctx = g_vpr_ctx.device();
//        auto& grid = device_ctx.grid;
//    
//    	int itype = type->index;
        auto type = ctx->belTypeFromId(cell->type);
    	auto itype = bel_types.at(type);

        int x_from, y_from;
        bool gb;
        ctx->estimatePosition(cell->bel, x_from, y_from, gb);

    	int rlx = std::min<float>(this->max_x, rlim);
    	int rly = std::min<float>(this->max_y, rlim); /* Added rly for aspect_ratio != 1 case. */
    	unsigned active_area = 4 * rlx * rly;
    
    	int min_x = std::max<float>(0, x_from - rlx);
    	int max_x = std::min<float>(this->max_x, x_from + rlx);
    	int min_y = std::max<float>(0, y_from - rly);
    	int max_y = std::min<float>(this->max_y, y_from + rly);
    
    	//*pz_to = 0;
    	if (int(max_x / 4) < rlx || int(max_y / 4) < rly || free_locations[itype].size() < active_area) {
    		int ipos = ctx->rng(free_locations[itype].size());
            bel = free_locations[itype][ipos];
//    		*px_to = [itype][ipos].x;
//    		*py_to = [itype][ipos].y;
//    		*pz_to = [itype][ipos].z;
    	} else {
    		int x_rel = ctx->rng(std::max(0, max_x - min_x)+1);
    		int y_rel = ctx->rng(std::max(0, max_y - min_y)+1);
//    		*px_to = min_x + x_rel;
//    		*py_to = min_y + y_rel;
//    		*px_to = (*px_to) - grid[*px_to][*py_to].width_offset; /* align it */
//    		*py_to = (*py_to) - grid[*px_to][*py_to].height_offset; /* align it */
            int px_to = min_x + x_rel;
            int py_to = min_y + y_rel;
            if (px_to >= int(fast_bels.at(itype).size())) {
                bel = BelId();
                return;
            }
            if (py_to >= int(fast_bels.at(itype).at(px_to).size())) {
                bel = BelId();
                return;
            }
            const auto &fb = fast_bels.at(itype).at(px_to).at(py_to);
            if (fb.size() == 0) {
                bel = BelId();
                return;
            }
            bel = fb.at(ctx->rng(int(fb.size())));
            // TODO: Remove locked_bels from fb
            if (locked_bels.find(bel) != locked_bels.end()) {
                bel = BelId();
                return;
            }
    	}
    }

    bool vpr_assess_swap(/*float delta_c,*/ float t) {
    
//    	/* Returns: 1 -> move accepted, 0 -> rejected. */
//    
    	bool accept;
    	float prob_fac, fnum;
    
    	if (delta_c <= 0) {
    
            /* Reduce variation in final solution due to round off */
            fnum = ctx->rng() / float(0x3fffffff);
    
            accept = true;
            return (accept);
    	}
    
    	if (t == 0.)
    		return false;
    
        fnum = ctx->rng() / float(0x3fffffff);
    	prob_fac = std::exp(-delta_c / t);
    	if (prob_fac > fnum) {
    		accept = true;
    	}
    	else {
    		accept = false;
    	}
    	return (accept);
    }

    bool vpr_exit_crit(float t, float cost) {
//     	/* Return 1 when the exit criterion is met.                        */
//    
//    	if (annealing_sched.type == USER_SCHED) {
//    		if (t < annealing_sched.exit_t) {
//    			return (1);
//    		} else {
//    			return (0);
//    		}
//    	}
//    
//        auto& cluster_ctx = g_vpr_ctx.clustering();

    	/* Automatic annealing schedule */
        float t_exit = 0.005 * cost / ctx->nets.size();
    
        if (t < t_exit) {
    		return (1);
        } else if (std::isnan(t_exit)) {
            //May get nan if there are no nets
            return (1);
    	} else {
    		return (0);
    	}
    }

    /* Update the temperature according to the annealing schedule selected. */
    void vpr_update_t(float success_rat) {
    
    	/*  float fac; */
    
//    	if (annealing_sched.type == USER_SCHED) {
//    		*t = annealing_sched.alpha_t * (*t);
//    	} else 
        { /* AUTO_SCHED */
    		if (success_rat > 0.96) {
    			t *= 0.5;
    		} else if (success_rat > 0.8) {
    			t *= 0.9;
    		} else if (success_rat > 0.15 || rlim > 1.) {
    			t *= 0.95;
    		} else {
    			t *= 0.8;
    		}
    	}
    }

    void vpr_find_affected_nets_and_update_costs(CellInfo* cell_from, CellInfo* cell_to, BelId bel_from, BelId bel_to, /*e_place_algorithm place_algorithm,*/ float& bb_delta_c /*, float& timing_delta_c, float& delay_delta_c*/) {
//        VTR_ASSERT_SAFE(bb_delta_c == 0.);
//        VTR_ASSERT_SAFE(timing_delta_c == 0.);
//        VTR_ASSERT_SAFE(delay_delta_c == 0.);
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//    
//    	int num_affected_nets = 0;
        affected_nets.clear();
    
//    	//Go through all the blocks moved
//    	for (int iblk = 0; iblk < blocks_affected.num_moved_blocks; iblk++) {
//    		ClusterBlockId blk = blocks_affected.moved_blocks[iblk].block_num;
//    
//    		//Go through all the pins in the moved block
//            for (ClusterPinId blk_pin : cluster_ctx.clb_nlist.block_pins(blk)) {
//    			ClusterNetId net_id = cluster_ctx.clb_nlist.pin_net(blk_pin);
//                VTR_ASSERT_SAFE_MSG(net_id, "Only valid nets should be found in compressed netlist block pins");
//    
//    			if (cluster_ctx.clb_nlist.net_is_global(net_id))
//    				continue; //Global nets are assumed to span the whole chip, and do not effect costs
//    
//                //Record effected nets
//                record_affected_net(net_id, num_affected_nets);
//    
//                //Update the net bounding boxes
//                //
//                //Do not update the net cost here since it should only be updated
//                //once per net, not once per pin.
//                update_net_bb(net_id, iblk, blk, blk_pin);
//    
//                if (place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
//                    //Determine the change in timing costs if required
//                    update_td_delta_costs(net_id, blk_pin, timing_delta_c, delay_delta_c);
//                }
//    		}
//    	}

        for (const auto &port : cell_from->ports)
            if (port.second.net != nullptr)
                affected_nets.insert(port.second.net);

        if (cell_to) {
            for (const auto &port : cell_to->ports)
                if (port.second.net != nullptr)
                    affected_nets.insert(port.second.net);
        }

        ctx->bindBel(bel_to, cell_from->name, STRENGTH_WEAK);
        if (cell_to) {
            ctx->bindBel(bel_from, cell_to->name, STRENGTH_WEAK);
        }

    
//        /* Now update the bounding box costs (since the net bounding boxes are up-to-date).
//         * The cost is only updated once per net.
//         */
//        for (int inet_affected = 0; inet_affected < num_affected_nets; inet_affected++) {
//            ClusterNetId net_id = ts_nets_to_update[inet_affected];
//    
//            temp_net_cost[net_id] = get_net_cost(net_id, &ts_bb_coord_new[net_id]);
//            bb_delta_c += temp_net_cost[net_id] - net_cost[net_id];
//        }
    
        auto new_wirelength = curr_wirelength;

        // Recalculate wirelengths for all nets touched by the peturbation
        for (auto net : affected_nets) {
            new_wirelength -= wirelengths.at(net->name);
            float temp_tns = 0;
            wirelen_t net_new_wl = get_net_wirelength(ctx, net, temp_tns);
            new_wirelength += net_new_wl;
            new_lengths.push_back(std::make_pair(net->name, net_new_wl));
        }
        bb_delta_c = new_wirelength - curr_wirelength;

//    	return num_affected_nets;
    }


    // Attempt a SA position swap, return true on success or false on failure
    bool try_swap_position(CellInfo *cell, BelId newBel)
    {
        static std::unordered_set<NetInfo *> update;
        static std::vector<std::pair<IdString, wirelen_t>> new_lengths;
        new_lengths.clear();
        update.clear();
        BelId oldBel = cell->bel;
        IdString other = ctx->getBoundBelCell(newBel);
        CellInfo *other_cell = nullptr;
        if (other != IdString()) {
            other_cell = ctx->cells[other].get();
            if (other_cell->belStrength > STRENGTH_WEAK)
                return false;
        }
        wirelen_t new_wirelength = 0, delta;
        ctx->unbindBel(oldBel);
        if (other != IdString()) {
            ctx->unbindBel(newBel);
        }

        for (const auto &port : cell->ports)
            if (port.second.net != nullptr)
                update.insert(port.second.net);

        if (other != IdString()) {
            for (const auto &port : other_cell->ports)
                if (port.second.net != nullptr)
                    update.insert(port.second.net);
        }

        ctx->bindBel(newBel, cell->name, STRENGTH_WEAK);

        if (other != IdString()) {
            ctx->bindBel(oldBel, other_cell->name, STRENGTH_WEAK);
        }
        if (require_legal) {
            if (!ctx->isBelLocationValid(newBel) || ((other != IdString() && !ctx->isBelLocationValid(oldBel)))) {
                ctx->unbindBel(newBel);
                if (other != IdString())
                    ctx->unbindBel(oldBel);
                goto swap_fail;
            }
        }

        new_wirelength = curr_wirelength;

        // Recalculate wirelengths for all nets touched by the peturbation
        for (auto net : update) {
            new_wirelength -= wirelengths.at(net->name);
            float temp_tns = 0;
            wirelen_t net_new_wl = get_net_wirelength(ctx, net, temp_tns);
            new_wirelength += net_new_wl;
            new_lengths.push_back(std::make_pair(net->name, net_new_wl));
        }
        delta = new_wirelength - curr_wirelength;
        n_move++;
        // SA acceptance criterea
        if (delta < 0 || (t > 1e-6 && (ctx->rng() / float(0x3fffffff)) <= std::exp(-delta / t))) {
            n_accept++;
            if (delta < 2)
                improved = true;
        } else {
            if (other != IdString())
                ctx->unbindBel(oldBel);
            ctx->unbindBel(newBel);
            goto swap_fail;
        }
        curr_wirelength = new_wirelength;
        for (auto new_wl : new_lengths)
            wirelengths.at(new_wl.first) = new_wl.second;

        return true;
    swap_fail:
        ctx->bindBel(oldBel, cell->name, STRENGTH_WEAK);
        if (other != IdString()) {
            ctx->bindBel(newBel, other, STRENGTH_WEAK);
        }
        return false;
    }

    // Find a random Bel of the correct type for a cell, within the specified
    // diameter
    BelId random_bel_for_cell(CellInfo *cell)
    {
        BelType targetType = ctx->belTypeFromId(cell->type);
        int x, y;
        bool gb;
        ctx->estimatePosition(cell->bel, x, y, gb);
        while (true) {
            int nx = ctx->rng(2 * diameter + 1) + std::max(x - diameter, 0);
            int ny = ctx->rng(2 * diameter + 1) + std::max(y - diameter, 0);
            int beltype_idx = bel_types.at(targetType);
            if (nx >= int(fast_bels.at(beltype_idx).size()))
                continue;
            if (ny >= int(fast_bels.at(beltype_idx).at(nx).size()))
                continue;
            const auto &fb = fast_bels.at(beltype_idx).at(nx).at(ny);
            if (fb.size() == 0)
                continue;
            BelId bel = fb.at(ctx->rng(int(fb.size())));
            if (locked_bels.find(bel) != locked_bels.end())
                continue;
            return bel;
        }
    }

    Context *ctx;
    std::unordered_map<IdString, wirelen_t> wirelengths;
    wirelen_t curr_wirelength = std::numeric_limits<wirelen_t>::max();
    float curr_tns = 0;
    float t = 1000;
    bool improved = false;
    int n_move, n_accept;
    int diameter = 35, max_x = 1, max_y = 1;
    std::unordered_map<BelType, int> bel_types;
    std::vector<std::vector<std::vector<std::vector<BelId>>>> fast_bels;
    std::unordered_set<BelId> locked_bels;
    bool require_legal = false;
    const float legalise_temp = 1;
    const float post_legalise_temp = 20;
    const float post_legalise_dia_scale = 2;
    std::vector<CellInfo *> autoplaced;

    const float inner_num = 1.0;
    int move_lim, tot_iter;
    float rlim;
    float cost, bb_cost;
    float delta_c, bb_delta_c;
    float success_sum;
    float success_rat;
    std::unordered_set<NetInfo *> affected_nets;
    std::vector<std::pair<IdString, wirelen_t>> new_lengths;
    int num_swap_accepted, num_swap_rejected;
};

bool place_design_vpr(Context *ctx)
{
    try {
        VPRPlacer placer(ctx);
        placer.place();
        log_info("Checksum: 0x%08x\n", ctx->checksum());
#ifndef NDEBUG
        ctx->check();
#endif
        return true;
    } catch (log_execution_error_exception) {
#ifndef NDEBUG
        ctx->check();
#endif
        return false;
    }
}

NEXTPNR_NAMESPACE_END
