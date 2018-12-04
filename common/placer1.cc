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

#include "placer1.h"
#include <algorithm>
#include <boost/lexical_cast.hpp>
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
#include "timing.h"
#include "util.h"
#include <boost/timer/timer.hpp>
#include <fstream>

#include "z3++.h"
using namespace z3;

NEXTPNR_NAMESPACE_BEGIN

class SAPlacer
{
  public:
    SAPlacer(Context *ctx, Placer1Cfg cfg) : ctx(ctx), cfg(cfg)
    {
        int num_bel_types = 0;
        for (auto bel : ctx->getBels()) {
            IdString type = ctx->getBelType(bel);
            if (bel_types.find(type) == bel_types.end()) {
                bel_types[type] = std::tuple<int, int>(num_bel_types++, 1);
            } else {
                std::get<1>(bel_types.at(type))++;
            }
        }
        for (auto bel : ctx->getBels()) {
            Loc loc = ctx->getBelLocation(bel);
            IdString type = ctx->getBelType(bel);
            int type_idx = std::get<0>(bel_types.at(type));
            int type_cnt = std::get<1>(bel_types.at(type));
            if (type_cnt < cfg.minBelsForGridPick)
                loc.x = loc.y = 0;
            if (int(fast_bels.size()) < type_idx + 1)
                fast_bels.resize(type_idx + 1);
            if (int(fast_bels.at(type_idx).size()) < (loc.x + 1))
                fast_bels.at(type_idx).resize(loc.x + 1);
            if (int(fast_bels.at(type_idx).at(loc.x).size()) < (loc.y + 1))
                fast_bels.at(type_idx).at(loc.x).resize(loc.y + 1);
            max_x = std::max(max_x, loc.x);
            max_y = std::max(max_y, loc.y);
            fast_bels.at(type_idx).at(loc.x).at(loc.y).push_back(bel);
        }
        diameter = std::max(max_x, max_y) + 1;

        costs.resize(ctx->nets.size());
        old_udata.reserve(ctx->nets.size());
        decltype(NetInfo::udata) n = 0;
        for (auto &net : ctx->nets) {
            old_udata.emplace_back(net.second->udata);
            net.second->udata = n++;
        }
    }

    ~SAPlacer()
    {
        for (auto &net : ctx->nets)
            net.second->udata = old_udata[net.second->udata];
    }

    bool place()
    {
        log_break();
        ctx->lock();

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

                IdString bel_type = ctx->getBelType(bel);
                if (bel_type != cell->type) {
                    log_error("Bel \'%s\' of type \'%s\' does not match cell "
                              "\'%s\' of type \'%s\'\n",
                              loc_name.c_str(), bel_type.c_str(ctx), cell->name.c_str(ctx), cell->type.c_str(ctx));
                }
                if (!ctx->isValidBelForCell(cell, bel)) {
                    log_error("Bel \'%s\' of type \'%s\' is not valid for cell "
                              "\'%s\' of type \'%s\'\n",
                              loc_name.c_str(), bel_type.c_str(ctx), cell->name.c_str(ctx), cell->type.c_str(ctx));
                }

                auto bound_cell = ctx->getBoundBelCell(bel);
                if (bound_cell) {
                    log_error("Cell \'%s\' cannot be bound to bel \'%s\' since it is already bound to cell \'%s\'\n",
                              cell->name.c_str(ctx), loc_name.c_str(), bound_cell->name.c_str(ctx));
                }

                ctx->bindBel(bel, cell, STRENGTH_USER);
                locked_bels.insert(bel);
                placed_cells++;
            }
        }
        int constr_placed_cells = placed_cells;
        log_info("Placed %d cells based on constraints.\n", int(placed_cells));
        ctx->yield();

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

        // Place cells randomly initially
        log_info("Creating initial placement for remaining %d cells.\n", int(autoplaced.size()));

        // Assume there are 
        //   i0...iN (unconstrained) cells
        //   j0...jN (available) bels

        context z3;
        solver s(z3);
        std::unordered_map<BelId, std::pair<expr_vector,std::vector<CellInfo*>>> placement_by_bel;
        std::unordered_map<Loc, std::pair<expr,expr>> clk_by_tile;
        std::unordered_map<Loc, std::pair<expr,expr>> cen_by_tile;
        std::unordered_map<Loc, std::pair<expr,expr>> sr_by_tile;
        std::unordered_map<Loc, expr> inputs_by_tile;
        const std::string delim = ".=>.";
        std::stringstream ss;
        std::unordered_map<int, int> clk2int;
        std::unordered_map<unsigned, int> cen2int;
        std::unordered_map<unsigned, int> sr2int;
        for (auto cell : autoplaced) {
            if (!cell->lcInfo.dffEnable) continue;

            if (cell->lcInfo.negClk)
                clk2int.emplace(-cell->lcInfo.clk->name.index, clk2int.size());
            else
                clk2int.emplace(cell->lcInfo.clk->name.index, clk2int.size());

            if (cell->lcInfo.cen)
                cen2int.emplace(cell->lcInfo.cen->name.index, cen2int.size()+1);

            if (cell->lcInfo.sr)
                sr2int.emplace(cell->lcInfo.sr->name.index, sr2int.size()+1);
        }
        // I build a matrix of bools sized iN * jN,
        // each bool (p_{i,j}) representing the placement
        // of cell_i into bel_j
        auto allBels = ctx->getBels();
        auto numBels = allBels.e.cursor-allBels.b.cursor;
        std::vector<int> all_ones(numBels, 1);
        auto bv_one = z3.bv_val(1, 1);
        expr_vector placement(z3);
        for (auto cell : autoplaced) {
            expr cell_placement = z3.bv_const(cell->name.c_str(ctx), numBels);
            expr_vector one_bel_per_cell(z3);
            for (auto bel : allBels) {
                auto e = cell_placement.extract(bel.index, bel.index);
                // Eliminate any invalid cell-bel pairs
                if (ctx->getBelType(bel) != cell->type) {
                    s.add(e != bv_one);
                    continue;
                }
                if (!ctx->isValidBelForCell(cell, bel)) {
                    s.add(e != bv_one);
                    continue;
                }
                // Add it to a vector that will later enforce the
                // constraint that each cell must be placed onto
                // one and only one bel
                one_bel_per_cell.push_back(e == bv_one);
                auto it = placement_by_bel.emplace(bel, std::make_pair(expr_vector{z3}, std::vector<CellInfo*>{})).first;
                it->second.first.push_back(e == bv_one);
                it->second.second.push_back(cell);
                // Now encode the tile constraints,
                // which are only relevant if DFFs are used
                if (cell->lcInfo.dffEnable) {
                    auto loc = ctx->getBelLocation(bel);
                    loc.z = 0;
                    // Constraint that all DFF cells in the same tile must
                    // have the same clock net and polarity
                    if (clk2int.size() > 1) {
                        auto jt = clk_by_tile.find(loc);
                        if (jt == clk_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".clk";
                            jt = clk_by_tile.emplace(loc, std::make_pair(z3.bv_const(ss.str().c_str(), static_cast<int>(ceil(log2(clk2int.size())))), z3.bool_const((ss.str()+"_b").c_str()))).first;
                        }
                        assert(cell->lcInfo.clk);
                        // I use an implication here to say that the "clock" variable of
                        // each tile must be set to the clock net's unique identifier,
                        // thus placing a cell from a different clock will create a 
                        // conflicting implication
                        if (cell->lcInfo.negClk)
                            s.add(implies(e == bv_one, jt->second.first == clk2int.at(-cell->lcInfo.clk->name.index)));
                        else
                            s.add(implies(e == bv_one, jt->second.first == clk2int.at(cell->lcInfo.clk->name.index)));
                        // I also set an implication that the clock is being used
                        // so that we can count it as a tile input later
                        if (cell->lcInfo.clk->is_global)
                            s.add(implies(e == bv_one, jt->second.second));
                    }
                    // Constraint that all DFF cells in the same tile must
                    // have the same clock-enable net (or none at all)
                    if (!cen2int.empty())
                    {
                        auto jt = cen_by_tile.find(loc);
                        if (jt == cen_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".cen";
                            jt = cen_by_tile.emplace(loc, std::make_pair(z3.bv_const(ss.str().c_str(), static_cast<int>(ceil(log2(cen2int.size()+1)))), z3.bool_const((ss.str()+"_b").c_str()))).first;
                        }
                        if (cell->lcInfo.cen) {
                            s.add(implies(e == bv_one, jt->second.first == cen2int.at(cell->lcInfo.cen->name.index)));
                            if (cell->lcInfo.cen->is_global)
                                s.add(implies(e == bv_one, jt->second.second));
                        }
                        else
                            s.add(implies(e == bv_one, jt->second.first == 0));
                    }
                    // Constraint that all DFF cells in the same tile must
                    // have the same set-reset net (or none at all)
                    if (!sr2int.empty())
                    {
                        auto jt = sr_by_tile.find(loc);
                        if (jt == sr_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".sr";
                            jt = sr_by_tile.emplace(loc, std::make_pair(z3.bv_const(ss.str().c_str(), static_cast<int>(ceil(log2(sr2int.size())))), z3.bool_const((ss.str()+"_b").c_str()))).first;
                        }
                        if (cell->lcInfo.sr) {
                            s.add(implies(e == bv_one, jt->second.first == sr2int.at(cell->lcInfo.sr->name.index)));
                            if (cell->lcInfo.sr->is_global)
                                s.add(implies(e == bv_one, jt->second.second));
                        }
                        else
                            s.add(implies(e == bv_one, jt->second.first == 0));
                    }
                }
            }
            placement.push_back(cell_placement);
            // Constraint that each cell_i must have
            // exactly one p_{i,j} set for all j
            s.add(pbeq(one_bel_per_cell, all_ones.data(), 1));
        }
        // Constraint that each bel_j must have
        // at most one p_{i,j} set for all i
        s.add(distinct(placement));
        //try {
        //    s.add(mk_and(placement) == z3.bv_val(0, numBels));
        //}
        //catch(const z3::exception& e) {
        //    std::cout << e << std::endl;
        //    throw;
        //}
        // Now constrain the number of local inputs
        // of each tile to be less than 32
        for (auto i : clk_by_tile) {
            const auto &loc = i.first;
            expr_vector tile_placement(z3);
            std::vector<int> input_count;
            for (auto b : ctx->getBelsByTile(loc.x, loc.y)) {
                for (unsigned j = 0; j < placement_by_bel.at(b).first.size(); ++j) {
                    // For every p_{i,j} in the same tile
                    auto &p = placement_by_bel.at(b);
                    auto e = p.first[j];
                    auto cell = p.second[j];
                    // Add the bool to all possible placements in the local tile
                    tile_placement.push_back(e);
                    // Add the number of inputs that placement will incur
                    input_count.push_back(cell->lcInfo.inputCount);
                }
            }
            // Now account for clk, cen, sr inputs
            tile_placement.push_back(clk_by_tile.at(i.first).second);
            tile_placement.push_back(cen_by_tile.at(i.first).second);
            tile_placement.push_back(sr_by_tile.at(i.first).second);
            input_count.push_back(1);
            input_count.push_back(1);
            input_count.push_back(1);
            // Check that the total is less than 33
            //   where total is cell_i.inputs     * p_{i,0}   + cell_i.inputs     * p_{i,1) + ...
            //                  cell_{i+1}.inputs * p_{i+1,0} + cell_{i+1}.inputs * p_{i+1,1} + ...
            //                  any_clk_used * 1 + any_cen_used * 1 + any_sr_used * 1
            //                  < 33
            s.add(pble(tile_placement, input_count.data(), 32));
        }
        //std::cout << "|cells| * |bels| = " << placement.size() << std::endl;
        //std::ofstream fs("last.smt2");
        //fs << "(set-logic ALL)" << std::endl;
        //fs << s.to_smt2() << std::endl;
        //fs.close();

        set_param("verbose", 10);
        boost::timer::cpu_timer timer;
        std::cout << s.check() << "\n";
        std::cout << timer.format() << std::endl;
        std::cout << s.statistics() << "\n";

        model m = s.get_model();
        //std::cout << m << "\n";
        // traversing the model
        for (unsigned i = 0; i < m.size(); i++) {
            func_decl v = m[i];
            // this problem contains only constants
            assert(v.arity() == 0); 
            auto vname = v.name().str();
            auto e = m.get_const_interp(v);
            CellInfo* cell;
            BelId bel;
            if (e.is_bv()) {
                auto it = ctx->cells.find(ctx->id(vname));
                if (it == ctx->cells.end()) continue;
                cell = it->second.get();
                for (int i = 0; i < numBels; ++i) {
                    if ((e.extract(i,i) == bv_one).simplify().is_true()) {
                        bel.index = i;
                        break;
                    }
                }
            }
            else if (e.is_true()) {
                auto it = vname.find(delim);
                if (it == std::string::npos) continue;
                cell = ctx->cells.at(ctx->id(vname.substr(0, it))).get();
                bel = ctx->getBelByName(ctx->id(vname.substr(it+delim.size(), std::string::npos)));
            }
            assert(ctx->isValidBelForCell(cell, bel));
            ctx->bindBel(bel, cell, STRENGTH_WEAK);
        }

        for (auto bel : ctx->getBels()) {
            CellInfo *cell = ctx->getBoundBelCell(bel);
            if (!ctx->isBelLocationValid(bel)) {
                std::string cell_text = "no cell";
                if (cell != nullptr)
                    cell_text = std::string("cell '") + ctx->nameOf(cell) + "'";
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

        if (ctx->slack_redist_iter > 0)
            assign_budget(ctx);
        ctx->yield();

        log_info("Running simulated annealing placer.\n");

        // Calculate metric after initial placement
        curr_metric = 0;
        curr_tns = 0;
        for (auto &net : ctx->nets) {
            wirelen_t wl = get_net_metric(ctx, net.second.get(), MetricType::COST, curr_tns);
            costs[net.second->udata] = CostChange{wl, -1};
            curr_metric += wl;
        }

        int n_no_progress = 0;
        wirelen_t min_metric = curr_metric;
        double avg_metric = curr_metric;
        temp = 10000;

        // Main simulated annealing loop
        for (int iter = 1;; iter++) {
            n_move = n_accept = 0;
            improved = false;

            if (iter % 5 == 0 || iter == 1)
                log_info("  at iteration #%d: temp = %f, cost = "
                         "%.0f, est tns = %.02fns\n",
                         iter, temp, double(curr_metric), curr_tns);

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

            if (curr_metric < min_metric) {
                min_metric = curr_metric;
                improved = true;
            }

            // Heuristic to improve placement on the 8k
            if (improved)
                n_no_progress = 0;
            else
                n_no_progress++;

            if (temp <= 1e-3 && n_no_progress >= 5) {
                if (iter % 5 != 0)
                    log_info("  at iteration #%d: temp = %f, cost = %f\n", iter, temp, double(curr_metric));
                break;
            }

            double Raccept = double(n_accept) / double(n_move);

            int M = std::max(max_x, max_y) + 1;

            double upper = 0.6, lower = 0.4;

            if (curr_metric < 0.95 * avg_metric) {
                avg_metric = 0.8 * avg_metric + 0.2 * curr_metric;
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
            if (temp < legalise_temp && require_legal) {
                if (legalise_relative_constraints(ctx)) {
                    // Only increase temperature if something was moved
                    autoplaced.clear();
                    for (auto cell : sorted(ctx->cells)) {
                        if (cell.second->belStrength < STRENGTH_STRONG)
                            autoplaced.push_back(cell.second);
                    }
                    temp = post_legalise_temp;
                    diameter *= post_legalise_dia_scale;
                    ctx->shuffle(autoplaced);

                    // Legalisation is a big change so force a slack redistribution here
                    if (ctx->slack_redist_iter > 0)
                        assign_budget(ctx, true /* quiet */);
                }
                require_legal = false;
            } else if (ctx->slack_redist_iter > 0 && iter % ctx->slack_redist_iter == 0) {
                assign_budget(ctx, true /* quiet */);
            }

            // Recalculate total metric entirely to avoid rounding errors
            // accumulating over time
            curr_metric = 0;
            curr_tns = 0;
            for (auto &net : ctx->nets) {
                wirelen_t wl = get_net_metric(ctx, net.second.get(), MetricType::COST, curr_tns);
                costs[net.second->udata] = CostChange{wl, -1};
                curr_metric += wl;
            }

            // Let the UI show visualization updates.
            ctx->yield();
        }
        // Final post-pacement validitiy check
        ctx->yield();
        for (auto bel : ctx->getBels()) {
            CellInfo *cell = ctx->getBoundBelCell(bel);
            if (!ctx->isBelLocationValid(bel)) {
                std::string cell_text = "no cell";
                if (cell != nullptr)
                    cell_text = std::string("cell '") + ctx->nameOf(cell) + "'";
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
        for (auto cell : sorted(ctx->cells))
            if (get_constraints_distance(ctx, cell.second) != 0)
                log_error("constraint satisfaction check failed for cell '%s' at Bel '%s'\n", cell.first.c_str(ctx),
                          ctx->getBelName(cell.second->bel).c_str(ctx));
        timing_analysis(ctx);
        ctx->unlock();
        return true;
    }

  private:
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
            IdString targetType = cell->type;
            for (auto bel : ctx->getBels()) {
                if (ctx->getBelType(bel) == targetType && ctx->isValidBelForCell(cell, bel)) {
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
                            ripup_target = ctx->getBoundBelCell(bel);
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
            ctx->bindBel(best_bel, cell, STRENGTH_WEAK);

            // Back annotate location
            cell->attrs[ctx->id("BEL")] = ctx->getBelName(cell->bel).str(ctx);
            cell = ripup_target;
        }
    }

    // Attempt a SA position swap, return true on success or false on failure
    bool try_swap_position(CellInfo *cell, BelId newBel)
    {
        static std::vector<NetInfo *> updates;
        updates.clear();
        BelId oldBel = cell->bel;
        CellInfo *other_cell = ctx->getBoundBelCell(newBel);
        if (other_cell != nullptr && other_cell->belStrength > STRENGTH_WEAK) {
            return false;
        }
        int old_dist = get_constraints_distance(ctx, cell);
        int new_dist;
        if (other_cell != nullptr)
            old_dist += get_constraints_distance(ctx, other_cell);
        wirelen_t new_metric = 0, delta;
        ctx->unbindBel(oldBel);
        if (other_cell != nullptr) {
            ctx->unbindBel(newBel);
        }

        for (const auto &port : cell->ports) {
            if (port.second.net != nullptr) {
                auto &cost = costs[port.second.net->udata];
                if (cost.new_cost == 0)
                    continue;
                cost.new_cost = 0;
                updates.emplace_back(port.second.net);
            }
        }

        if (other_cell != nullptr) {
            for (const auto &port : other_cell->ports)
                if (port.second.net != nullptr) {
                    auto &cost = costs[port.second.net->udata];
                    if (cost.new_cost == 0)
                        continue;
                    cost.new_cost = 0;
                    updates.emplace_back(port.second.net);
                }
        }

        ctx->bindBel(newBel, cell, STRENGTH_WEAK);

        if (other_cell != nullptr) {
            ctx->bindBel(oldBel, other_cell, STRENGTH_WEAK);
        }
        if (!ctx->isBelLocationValid(newBel) || ((other_cell != nullptr && !ctx->isBelLocationValid(oldBel)))) {
            ctx->unbindBel(newBel);
            if (other_cell != nullptr)
                ctx->unbindBel(oldBel);
            goto swap_fail;
        }

        new_metric = curr_metric;

        // Recalculate metrics for all nets touched by the peturbation
        for (const auto &net : updates) {
            auto &c = costs[net->udata];
            new_metric -= c.curr_cost;
            float temp_tns = 0;
            wirelen_t net_new_wl = get_net_metric(ctx, net, MetricType::COST, temp_tns);
            new_metric += net_new_wl;
            c.new_cost = net_new_wl;
        }

        new_dist = get_constraints_distance(ctx, cell);
        if (other_cell != nullptr)
            new_dist += get_constraints_distance(ctx, other_cell);
        delta = new_metric - curr_metric;
        delta += (cfg.constraintWeight / temp) * (new_dist - old_dist);
        n_move++;
        // SA acceptance criterea
        if (delta < 0 || (temp > 1e-6 && (ctx->rng() / float(0x3fffffff)) <= std::exp(-delta / temp))) {
            n_accept++;
        } else {
            if (other_cell != nullptr)
                ctx->unbindBel(oldBel);
            ctx->unbindBel(newBel);
            goto swap_fail;
        }
        curr_metric = new_metric;
        for (const auto &net : updates) {
            auto &c = costs[net->udata];
            c = CostChange{c.new_cost, -1};
        }

        return true;
    swap_fail:
        ctx->bindBel(oldBel, cell, STRENGTH_WEAK);
        if (other_cell != nullptr) {
            ctx->bindBel(newBel, other_cell, STRENGTH_WEAK);
        }
        for (const auto &net : updates)
            costs[net->udata].new_cost = -1;
        return false;
    }

    // Find a random Bel of the correct type for a cell, within the specified
    // diameter
    BelId random_bel_for_cell(CellInfo *cell)
    {
        IdString targetType = cell->type;
        Loc curr_loc = ctx->getBelLocation(cell->bel);
        while (true) {
            int nx = ctx->rng(2 * diameter + 1) + std::max(curr_loc.x - diameter, 0);
            int ny = ctx->rng(2 * diameter + 1) + std::max(curr_loc.y - diameter, 0);
            int beltype_idx, beltype_cnt;
            std::tie(beltype_idx, beltype_cnt) = bel_types.at(targetType);
            if (beltype_cnt < cfg.minBelsForGridPick)
                nx = ny = 0;
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
    wirelen_t curr_metric = std::numeric_limits<wirelen_t>::max();
    float curr_tns = 0;
    float temp = 1000;
    bool improved = false;
    int n_move, n_accept;
    int diameter = 35, max_x = 1, max_y = 1;
    std::unordered_map<IdString, std::tuple<int, int>> bel_types;
    std::vector<std::vector<std::vector<std::vector<BelId>>>> fast_bels;
    std::unordered_set<BelId> locked_bels;
    bool require_legal = true;
    const float legalise_temp = 1;
    const float post_legalise_temp = 10;
    const float post_legalise_dia_scale = 1.5;
    Placer1Cfg cfg;

    struct CostChange
    {
        wirelen_t curr_cost;
        wirelen_t new_cost;
    };
    std::vector<CostChange> costs;
    std::vector<decltype(NetInfo::udata)> old_udata;
};

Placer1Cfg::Placer1Cfg(Context *ctx) : Settings(ctx)
{
    constraintWeight = get<float>("placer1/constraintWeight", 10);
    minBelsForGridPick = get<int>("placer1/minBelsForGridPick", 64);
}

bool placer1(Context *ctx, Placer1Cfg cfg)
{
    try {
        SAPlacer placer(ctx, cfg);
        placer.place();
        log_info("Checksum: 0x%08x\n", ctx->checksum());
#ifndef NDEBUG
        ctx->lock();
        ctx->check();
        ctx->unlock();
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
