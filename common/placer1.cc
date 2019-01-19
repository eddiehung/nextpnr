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
        //   c0...cM (unconstrained) cells
        //   s0...sN (available) sites/bels

        context z3;
        if (getenv("Z3_PARALLEL"))
            set_param("parallel.enable", true);
        if (getenv("Z3_RANDOM_SEED"))
            set_param("smt.random_seed", boost::lexical_cast<int>(getenv("Z3_RANDOM_SEED")));
        auto z3_logic = getenv("Z3_LOGIC");
        solver s = z3_logic ? solver{z3, z3_logic} : solver{z3};
        //optimize s(z3);
        std::unordered_map<CellInfo*, expr_vector> placement_by_cell;
        struct Loc_expr {
            expr x, y, z;
        };
        std::unordered_map<const CellInfo*, Loc_expr> cell_to_loc;
        std::unordered_map<Loc, expr> clk_by_tile, cen_by_tile, sr_by_tile;
        std::unordered_map<const NetInfo*, int> clk2index, nclk2index, cen2index, sr2index;

        const unsigned x_bits = ceil(log2(ctx->getGridDimX()));
        const unsigned y_bits = ceil(log2(ctx->getGridDimY()));
        const unsigned z_bits = ceil(log2(ctx->getTileBelDimZ(0,0)));

        std::stringstream ss;
        for (auto &cell : ctx->cells) {
            CellInfo *ci = cell.second.get();
            if (ci->bel != BelId()) {
                auto loc = ctx->getBelLocation(ci->bel);
                cell_to_loc.emplace(ci, Loc_expr{z3.bv_val(loc.x,  x_bits), z3.bv_val(loc.y, y_bits), z3.bv_val(loc.z, z_bits)});
            }
        }

        // Assign indices to all possible clocks+polarity/clock-enable/set-reset
        for (auto cell : autoplaced) {
            if (!cell->lcInfo.dffEnable) {
                cen2index.emplace(nullptr, 0);
                sr2index.emplace(nullptr, 0);
                continue;
            }

            assert(cell->lcInfo.clk);
            if (cell->lcInfo.negClk)
                nclk2index.emplace(cell->lcInfo.clk, clk2index.size() + nclk2index.size());
            else
                clk2index.emplace(cell->lcInfo.clk, clk2index.size() + nclk2index.size());

            if (cell->lcInfo.cen)
                cen2index.emplace(cell->lcInfo.cen, cen2index.size());

            if (cell->lcInfo.sr)
                sr2index.emplace(cell->lcInfo.sr, sr2index.size());
        }

        auto exactly1 = [&z3](const expr_vector &v) {
            /*static*/ std::function<std::pair<expr,expr>(const expr_vector&,unsigned,unsigned)> f = [&z3,&f](const expr_vector &v, unsigned b, unsigned e) {
                int dist = e - b;
                assert(dist > 0);
                if (dist == 1) return std::make_pair(v[b], z3.bv_val(0,1));
                if (dist == 2) return std::make_pair(v[b] ^ v[b+1], v[b] & v[b+1]);
                auto h = b + dist / 2;
                assert(h > b);
                auto l = f(v, b, h);
                auto r = f(v, h, e);
                return std::make_pair(l.first ^ r.first, l.second | r.second | (l.first & r.first));
            };
            auto r = f(v, 0, v.size());
            return r.second == z3.bv_val(0,1) && r.first == z3.bv_val(1,1);
        };

        auto exactly0 = [&z3](const expr_vector &v) {
            /*static*/ std::function<expr(const expr_vector&,unsigned,unsigned)> f = [&f](const expr_vector &v, unsigned b, unsigned e) {
                int dist = e - b;
                assert(dist > 0);
                if (dist == 1) return v[b];
                if (dist == 2) return v[b] | v[b+1];
                auto h = b + dist / 2;
                assert(h > b);
                return f(v, b, h) | f(v, h, e);
            };
            return f(v, 0, v.size()) == z3.bv_val(0,1);
        };

        const std::string delim{".->."};

        // I build a matrix of bools sized cM * sN,
        // each bool (x_{i,j}) representing the placement
        // of cell_i into site_j
        auto all_bels = ctx->getBels();
        std::unordered_map<BelId, expr_vector> placement_by_bel;
        for (auto cell : autoplaced) {
            ss.str("");
            ss << cell->name.str(ctx) << ".x";
            auto x = z3.bv_const(ss.str().c_str(), x_bits);
            ss.str("");
            ss << cell->name.str(ctx) << ".y";
            auto y = z3.bv_const(ss.str().c_str(), y_bits);
            ss.str("");
            ss << cell->name.str(ctx) << ".z";
            auto z = z3.bv_const(ss.str().c_str(), z_bits);
            cell_to_loc.emplace(cell, Loc_expr{x,y,z});

            //s.add(/*uge(x, 0) &&*/ ult(x, ctx->getGridDimX()));
            //s.add(/*uge(y, 0) &&*/ ult(y, ctx->getGridDimY()));
            //s.add(/*uge(z, 0) &&*/ ult(z, ctx->getTileBelDimZ(0,0)));

            for (auto bel : all_bels) {
                // Eliminate any invalid cell-bel pairs
                if (ctx->getBelType(bel) != cell->type)
                    continue;
                if (!ctx->isValidBelForCell(cell, bel))
                    continue;
                // Create an int variable named cell
                ss.str("");
                ss << cell->name.str(ctx);
                ss << delim;
                ss << ctx->getBelName(bel).str(ctx);
                auto p = z3.bv_const(ss.str().c_str(), 1);
                auto e = p == z3.bv_val(1,1);
                auto it = placement_by_cell.emplace(cell, expr_vector{z3}).first;
                it->second.push_back(p);
                auto jt = placement_by_bel.emplace(bel, expr_vector{z3}).first;
                jt->second.push_back(p);
                // Add it to a vector that will later enforce the
                // constraint that each cell must be placed onto
                // one and only one bel
                auto loc = ctx->getBelLocation(bel);
                s.add(implies(e, x == loc.x));
                s.add(implies(e, y == loc.y));
//                if (cell->type == id_ICESTORM_LC)
//                    s.add(implies(e, z == loc.z));
                // Now encode the tile constraints,
                // which are only relevant if DFFs are used
                if (cell->lcInfo.dffEnable) {
                    loc.z = 0;
                    // Constraint that all DFF cells in the same tile must
                    // have the same clock net and polarity
                    if ((clk2index.size() + nclk2index.size()) > 1) {
                        auto jt = clk_by_tile.find(loc);
                        if (jt == clk_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".clk";
                            auto c = z3.bv_const(ss.str().c_str(), std::max<int>(4,ceil(log2(clk2index.size()+nclk2index.size()))));
                            s.add(/*c >= 0 &&*/ ult(c,  clk2index.size()+nclk2index.size()));
                            jt = clk_by_tile.emplace(loc, c).first;
                        }
                        assert(cell->lcInfo.clk);
                        if (cell->lcInfo.negClk)
                            s.add(implies(e, jt->second == nclk2index.at(cell->lcInfo.clk)));
                        else
                            s.add(implies(e, jt->second == clk2index.at(cell->lcInfo.clk)));
                    }
                    // Constraint that all DFF cells in the same tile must
                    // have the same clock-enable net (or none at all)
                    if (cen2index.size() > 1) {
                        auto jt = cen_by_tile.find(loc);
                        if (jt == cen_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".cen";
                            auto c = z3.bv_const(ss.str().c_str(), std::max<int>(4,ceil(log2(cen2index.size()))));
                            s.add(/*c >= 0 &&*/ ult(c, cen2index.size()));
                            jt = cen_by_tile.emplace(loc, c).first;
                        }
                        s.add(implies(e, jt->second == cen2index.at(cell->lcInfo.cen)));
                    }
                    // Constraint that all DFF cells in the same tile must
                    // have the same set-reset net (or none at all)
                    if (sr2index.size() > 1) {
                        auto jt = sr_by_tile.find(loc);
                        if (jt == sr_by_tile.end()) {
                            ss.str("");
                            ss << "x" << loc.x << "y" << loc.y << ".sr";
                            auto c = z3.bv_const(ss.str().c_str(), std::max<int>(4, ceil(log2(sr2index.size()))));
                            //s.add(/*c >= 0 &&*/ ult(c, sr2index.size()));
                            s.add(/*c >= 0 &&*/ c < sr2index.size());
                            jt = sr_by_tile.emplace(loc, c).first;
                        }
                        s.add(implies(e, jt->second == sr2index.at(cell->lcInfo.sr)));
                    }
                }
            }
            s.add(exactly1(placement_by_cell.at(cell)));
        }

        for (auto i : placement_by_bel)
            s.add(exactly0(i.second) || exactly1(i.second));

        // TODO: Enforce constraint that each tile can have at most 32 inputs
        // non-global inputs -- this can only be an issue when non-global clk/cen/sr
        // nets are used
        
        if (ctx->timing_driven) {
            struct model_params_t
            {
                int neighbourhood;
            
                int model0_offset;
                int model0_norm1;
            
                int model1_offset;
                int model1_norm1;
                int model1_norm2;
                int model1_norm3;
            
                int model2_offset;
                int model2_linear;
                int model2_sqrt;
            
                int delta_local;
                int delta_lutffin;
                int delta_sp4;
                int delta_sp12;
            
                static const model_params_t &get(const ArchArgs &args)
                {
                    static const model_params_t model_hx8k = {588,    129253, 8658, 118333, 23915, -73105, 57696,
                                                              -86797, 89,     3706, -316,   -575,  -158,   -296};
            
                    static const model_params_t model_lp8k = {867,     206236, 11043, 191910, 31074, -95972, 75739,
                                                              -309793, 30,     11056, -474,   -856,  -363,   -536};
            
                    static const model_params_t model_up5k = {1761,    305798, 16705, 296830, 24430, -40369, 33038,
                                                              -162662, 94,     4705,  -1099,  -1761, -418,   -838};
            
                    if (args.type == ArchArgs::HX1K || args.type == ArchArgs::HX8K)
                        return model_hx8k;
            
                    if (args.type == ArchArgs::LP384 || args.type == ArchArgs::LP1K || args.type == ArchArgs::LP8K)
                        return model_lp8k;
            
                    if (args.type == ArchArgs::UP5K)
                        return model_up5k;
            
                    NPNR_ASSERT(0);
                }
            };
            const model_params_t &p = model_params_t::get(ctx->args);
            auto period = ctx->getDelayFromNS(1.0e9 / ctx->target_freq).maxDelay();
            auto min_slack = z3.bv_const("min_slack", 32);
            for (auto &net : ctx->nets) {
                auto driver = net.second->driver;
                CellInfo *driver_cell = driver.cell;
                if (!driver_cell)
                    continue;
                auto driver_loc = cell_to_loc.at(driver_cell);
                for (auto load : net.second->users) {
                    CellInfo *load_cell = load.cell;
                    if (load_cell == nullptr || load_cell == driver_cell)
                        continue;
                    if (load.budget < 0 || load.budget >= period)
                        continue;
                    if (load_cell->type != id_ICESTORM_LC)
                        continue;
                    /*if (ctx->timing_driven)*/ {
                        auto sink_loc = cell_to_loc.at(load_cell);

                        if (driver.port == id_COUT) {
                            continue; // FIXME
                            //auto delay = ite(dy == 0, z3.bv_val(0,32), z3.bv_val(190,32));
                            //s.add(delay >= 0 && delay <= load.budget);
                        }
                        else {
                            assert(load_cell->type == id_ICESTORM_LC);
                            auto adx = zext(ite(uge(sink_loc.x, driver_loc.x), sink_loc.x - driver_loc.x, driver_loc.x - sink_loc.x), 32-x_bits);
                            auto ady = zext(ite(uge(sink_loc.y, driver_loc.y), sink_loc.y - driver_loc.y, driver_loc.y - sink_loc.y), 32-y_bits);
                            auto delay = ite(ule(adx, 1) && ule(ady, 1), z3.bv_val(p.neighbourhood * 128, 32), p.model0_offset + p.model0_norm1 * (adx+ady));
                            auto slack = load.budget * 128 - delay;
                            s.add(min_slack <= slack);
                        }
                    }
                }
            }

#if 1
            s.add(min_slack >= 0);
#else
            //s.maximize(min_slack);
            params sp(z3);
            //sp.set("timeout", 120000u);
            s.set(sp);
#endif
        }

        for (auto cell : autoplaced) {
            auto parent = cell->constr_parent;
            if (!parent) continue;

            auto parent_loc = cell_to_loc.at(parent);
            auto this_loc = cell_to_loc.at(cell);
            s.add(this_loc.x == parent_loc.x + cell->constr_x);
            s.add(this_loc.y == parent_loc.y + cell->constr_y);
            if (cell->constr_abs_z)
                s.add(this_loc.z == cell->constr_z);
            else
                s.add(this_loc.z == parent_loc.z + cell->constr_z);
        }

        if (getenv("Z3_OUTPUT_SMT2")) {
            std::ofstream f(getenv("Z3_OUTPUT_SMT2"));
            f << s.to_smt2();
            f.close();
        }

        if (getenv("Z3_VERBOSITY"))
            set_param("verbose", boost::lexical_cast<int>(getenv("Z3_VERBOSITY")));
        boost::timer::cpu_timer timer;
        std::cout << s.check() << "\n";
        std::cout << timer.format() << std::endl;
        std::cout << s.statistics() << "\n";

        model m = s.get_model();
        //std::cout << m << "\n";
        // traversing the model
        for (const auto &i : placement_by_cell) {
            auto cell = i.first;
            const auto &bels = i.second;
            for (auto j = 0u; j < bels.size(); ++j) {
                if (m.eval(bels[j] == z3.bv_val(1,1)).is_false()) continue;
                auto n = bels[j].decl().name().str();
                auto it = n.find(delim);
                assert(it != std::string::npos);
                auto bel = ctx->getBelByName(ctx->id(n.substr(it+delim.size(), std::string::npos)));
                assert(bel != BelId());
                if (ctx->debug)
                    std::cout << bels[j].decl().name() << std::endl;
                assert(ctx->isValidBelForCell(cell, bel));
                ctx->bindBel(bel, cell, STRENGTH_WEAK);
                break;
            }
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

            break;

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
                        CellInfo *bound_cell = ctx->getBoundBelCell(bel);
                        if (score <= best_ripup_score && bound_cell->belStrength < STRENGTH_STRONG) {
                            best_ripup_score = score;
                            ripup_target = bound_cell;
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
