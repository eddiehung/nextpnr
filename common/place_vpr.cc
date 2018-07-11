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

class VPRPlacer
{
  private:
    std::vector<std::vector<BelId>> free_locations;

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
        std::vector<CellInfo *> autoplaced;
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

        vpr_initial_placement(autoplaced, placed_cells, constr_placed_cells);
#if 0
        for (auto cell : autoplaced) {
            place_initial(cell);
            placed_cells++;
            if ((placed_cells - constr_placed_cells) % 500 == 0)
                log_info("  initial placement placed %d/%d cells\n", int(placed_cells - constr_placed_cells),
                         int(autoplaced.size()));
        }
        if ((placed_cells - constr_placed_cells) % 500 != 0)
            log_info("  initial placement placed %d/%d cells\n", int(placed_cells - constr_placed_cells),
                     int(autoplaced.size()));
#endif

        log_info("Running simulated annealing placer.\n");

        // Calculate wirelength after initial placement
        curr_wirelength = 0;
        curr_tns = 0;
        for (auto &net : ctx->nets) {
            wirelen_t wl = get_net_wirelength(ctx, net.second.get(), curr_tns);
            wirelengths[net.first] = wl;
            curr_wirelength += wl;
        }

        int n_no_progress = 0;
        double avg_wirelength = curr_wirelength;
        temp = 10000;

        // Main simulated annealing loop
        for (int iter = 1;; iter++) {
            n_move = n_accept = 0;
            improved = false;

            if (iter % 5 == 0 || iter == 1)
                log_info("  at iteration #%d: temp = %f, wire length = "
                         "%.0f, est tns = %.02fns\n",
                         iter, temp, double(curr_wirelength), curr_tns);

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

            if (temp <= 1e-3 && n_no_progress >= 5) {
                if (iter % 5 != 0)
                    log_info("  at iteration #%d: temp = %f, wire length = %f\n", iter, temp, double(curr_wirelength));
                break;
            }

            double Raccept = double(n_accept) / double(n_move);

            int M = std::max(max_x, max_y) + 1;

            double upper = 0.6, lower = 0.4;

            if (curr_wirelength < 0.95 * avg_wirelength) {
                avg_wirelength = 0.8 * avg_wirelength + 0.2 * curr_wirelength;
            } else {
                if (Raccept >= 0.8) {
                    temp *= 0.7;
                } else if (Raccept > upper) {
                    if (diameter < M)
                        diameter++;
                    else
                        temp *= 0.9;
                } else if (Raccept > lower) {
                    temp *= 0.95;
                } else {
                    // Raccept < 0.3
                    if (diameter > 1)
                        diameter--;
                    else
                        temp *= 0.8;
                }
            }
            // Once cooled below legalise threshold, run legalisation and start requiring
            // legal moves only
            if (temp < legalise_temp && !require_legal) {
                legalise_design(ctx);
                require_legal = true;
                autoplaced.clear();
                for (auto cell : sorted(ctx->cells)) {
                    if (cell.second->belStrength < STRENGTH_STRONG)
                        autoplaced.push_back(cell.second);
                }
                temp = post_legalise_temp;
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

    void vpr_initial_placement(const std::vector<CellInfo *>& autoplaced, size_t& placed_cells, int constr_placed_cells) {
    
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
//    

        // free_locations is populated in constructor

//    	free_locations = (int *) vtr::malloc(device_ctx.num_block_types * sizeof(int));
//    	for (itype = 0; itype < device_ctx.num_block_types; itype++) {
//    		free_locations[itype] = num_legal_pos[itype];
//    	}
//    
//    	/* We'll use the grid to record where everything goes. Initialize to the grid has no
//    	 * blocks placed anywhere.
//    	 */
//    	for (size_t i = 0; i < device_ctx.grid.width(); i++) {
//    		for (size_t j = 0; j < device_ctx.grid.height(); j++) {
//    			place_ctx.grid_blocks[i][j].usage = 0;
//    			itype = device_ctx.grid[i][j].type->index;
//    			for (int k = 0; k < device_ctx.block_types[itype].capacity; k++) {
//    				if (place_ctx.grid_blocks[i][j].blocks[k] != INVALID_BLOCK_ID) {
//    					place_ctx.grid_blocks[i][j].blocks[k] = EMPTY_BLOCK_ID;
//    				}
//    			}
//    		}
//    	}
//    
//    	/* Similarly, mark all blocks as not being placed yet. */
//    	for (auto blk_id : cluster_ctx.clb_nlist.blocks()) {
//    		place_ctx.block_locs[blk_id].x = OPEN;
//    		place_ctx.block_locs[blk_id].y = OPEN;
//    		place_ctx.block_locs[blk_id].z = OPEN;
//    	}
//    
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
    
    	vpr_initial_placement_blocks(autoplaced, placed_cells, constr_placed_cells);
    
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
    void vpr_initial_placement_blocks(const std::vector<CellInfo *>& autoplaced, size_t &placed_cells, const int constr_placed_cells) {
//    	int itype, ipos, x, y, z;
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//        auto& place_ctx = g_vpr_ctx.mutable_placement();
//        auto& device_ctx = g_vpr_ctx.device();
    
    size_t itype;
    BelId bel;

    // Shuffle all free locations once here, rather than picking a block at random
    for (auto& free_locations_type : free_locations) {
        ctx->shuffle(free_locations_type);
    }

    for (auto cell : autoplaced) {
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
//    			}

    			vpr_initial_placement_location(cell, itype, bel);

                NPNR_ASSERT(ctx->checkBelAvail(bel));
                ctx->bindBel(bel, cell->name, STRENGTH_WEAK);

//                //Mark IOs as fixed if specifying a (fixed) random placement
//                if(is_io_type(cluster_ctx.clb_nlist.block_type(blk_idctx->belTypeFromId(cell->type))) && pad_loc_type == RANDOM) {
//                    place_ctx.block_locs[blk_id].is_fixed = true;
//     			}

//    			/* Ensure randomizer doesn't pick this location again, since it's occupied. Could shift all the
//    			* legal positions in legal_pos to remove the entry (choice) we just used, but faster to
//    			* just move the last entry in legal_pos to the spot we just used and decrement the
//    			* count of free_locations. */
//    			legal_pos[itype][ipos] = legal_pos[itype][free_locations[itype] - 1]; /* overwrite used block position */
//    			free_locations[itype]--;

                free_locations[itype].pop_back();

//    		}

            ++placed_cells;
            if ((placed_cells - constr_placed_cells) % 500 == 0)
                log_info("  initial placement placed %d/%d cells\n", int(placed_cells - constr_placed_cells),
                         int(autoplaced.size()));
    	}
        if ((placed_cells - constr_placed_cells) % 500 != 0)
            log_info("  initial placement placed %d/%d cells\n", int(placed_cells - constr_placed_cells),
                     int(autoplaced.size()));
    }

    void vpr_initial_placement_location(CellInfo *cell, size_t& itype, BelId& bel) {
//        auto& cluster_ctx = g_vpr_ctx.clustering();
//    
        auto type = ctx->belTypeFromId(cell->type);
    	itype = bel_types.at(type);

        bel = free_locations[itype].back();
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
        if (delta < 0 || (temp > 1e-6 && (ctx->rng() / float(0x3fffffff)) <= std::exp(-delta / temp))) {
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
    float temp = 1000;
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
