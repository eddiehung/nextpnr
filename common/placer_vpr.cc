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

#include "placer_vpr.h"
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
#include "cells.h"
#include "design_utils.h"

namespace NEXTPNR_NAMESPACE {
    struct CellChain
    {
        std::vector<CellInfo *> cells;
        float mid_x = 0, mid_y = 0;
    };
    std::vector<CellChain> all_chains;
}

namespace vpr {
    using namespace NEXTPNR_NAMESPACE;

    static Context* npnr_ctx = NULL;
    std::vector<CellInfo *> npnr_cells;

    struct DeviceGrid {
        inline int width() const { return _bels.size(); }
        inline int height() const { return _bels.front().size(); }
        inline const std::vector<std::vector<BelId>>& operator[](size_t x) { return _bels.at(x); }
        std::vector<std::vector<std::vector<BelId>>> _bels;
    };
    enum PinType
    {
        DRIVER = PortType::PORT_OUT,
        SINK = PortType::PORT_IN,
    };
    static struct {
        struct t_clustering {
            struct {
                std::unordered_map<IdString, std::unique_ptr<nextpnr_ice40::NetInfo>>& nets() const { return npnr_ctx->nets; }
                bool net_is_global(const NetInfo* net) const { return npnr_ctx->isGlobalNet(net); }
                std::vector<PortRef>& net_sinks(NetInfo* net) const { return net->users; }
                struct t_net_pins {
                    t_net_pins(const NetInfo* net) : _net(net) {}
                    size_t size() const { return _net->users.size() + 1; }
                    const NetInfo* _net;
                };
                t_net_pins net_pins(const NetInfo* net) const { return t_net_pins(net); }
                std::unordered_map<IdString, PortInfo>& block_pins(CellInfo* cell) { return cell->ports; }
                PinType pin_type(const PortInfo& port) { return static_cast<PinType>(port.type); }
                CellInfo* net_driver_block(const NetInfo* net) const { return net->driver.cell; }
                size_t pin_net_index(const PortInfo& port, const CellInfo* cell) {
                    auto net = port.net;
                    for (auto it = net->users.begin(); it != net->users.end(); ++it) {
                        if (it->port == port.name && it->cell == cell)
                            return it - net->users.begin() + 1;
                    }
                    throw;
                }
                std::unordered_map<IdString, std::unique_ptr<nextpnr_ice40::CellInfo>>& blocks() const { return npnr_ctx->cells; }
                IdString block_type(const CellInfo* cell) const { return cell->type; }
                std::string block_name(const CellInfo* cell) const { return cell->name.str(npnr_ctx); }
            } clb_nlist;
            t_clustering& operator()() { return *this; }
        } clustering;
        struct t_device {
            DeviceGrid grid;
            t_device& operator()() { return *this; }
        } device;
    } g_vpr_ctx;
    static struct t_annealing_sched {
        const float inner_num = 10;
    } annealing_sched;
    static struct t_placer_opts {
        bool enable_timing_computations;
        const int inner_loop_recompute_divider = 0;
        const int recompute_crit_iter = 1;
        const float td_place_exp_first = 1.0;
        const float td_place_exp_last = 8.0;
        const float timing_tradeoff = 0.5;
    } placer_opts;
    struct SetupTimingInfo {
        void update() { update_budget(npnr_ctx); }
    };
    static SetupTimingInfo timing_info;
    struct t_pl_macro_member {
        t_pl_macro_member(CellInfo* blk_index, int x_offset, int y_offset, int z_offset) : blk_index(blk_index), x_offset(x_offset), y_offset(y_offset), z_offset(z_offset) {}
        CellInfo* blk_index;
        int x_offset;
        int y_offset;
        int z_offset;
    };
    typedef std::vector<t_pl_macro_member> t_pl_macro;

    // timing_util.cpp
    float calculate_clb_net_pin_criticality(/*timing_info, pin_lookup,*/ const PortRef& load, const NetInfo* net)
    {
        NPNR_ASSERT(npnr_ctx->timing_driven);

#if 1
        int driver_x, driver_y;
        bool driver_gb;
        CellInfo *driver_cell = net->driver.cell;
        if (!driver_cell)
            return 0;
        if (driver_cell->bel == BelId())
            return 0;
        npnr_ctx->estimatePosition(driver_cell->bel, driver_x, driver_y, driver_gb);
        WireId drv_wire = npnr_ctx->getWireBelPin(driver_cell->bel, npnr_ctx->portPinFromId(net->driver.port));
        if (driver_gb)
            return 0;
        if (load.cell == nullptr)
            return 0;
        CellInfo *load_cell = load.cell;
        if (load_cell->bel == BelId())
            return 0;
        WireId user_wire = npnr_ctx->getWireBelPin(load_cell->bel, npnr_ctx->portPinFromId(load.port));
        delay_t raw_wl = npnr_ctx->estimateDelay(drv_wire, user_wire);
        float slack = npnr_ctx->getDelayNS(load.budget) - npnr_ctx->getDelayNS(raw_wl);
        if (slack <= 0)
            return 1 - slack;
        return 1/std::max<float>(1, slack);
#else
        return 1;
#endif
    }

    // place_macro.cpp
    int alloc_and_load_placement_macros(/*t_direct_inf* directs, int num_directs, t_pl_macro ** */ std::unordered_map<CellInfo*, t_pl_macro> &macros)
    {
        for (auto &chain : all_chains) {
            t_pl_macro entry;
            for (int z = 0; z < int(chain.cells.size()); z++)
                entry.emplace_back(chain.cells.at(z), 0, z / 8, z % 8);
            macros.emplace(chain.cells.front(), std::move(entry));
        }

        assign_budget(npnr_ctx);

        return macros.size();
    }
    
    #define VTR_ASSERT NPNR_ASSERT
    #define VTR_ASSERT_SAFE NPNR_ASSERT
    #define VTR_ASSERT_SAFE_MSG NPNR_ASSERT_MSG
    #define vpr_throw(__a, __b, __c, ...) log_error(__VA_ARGS__)
    namespace vtr {
        constexpr size_t bufsize = 32768;
        template <typename ...Args>
    	inline void printf_warning(const char*, unsigned, const char* fmt, Args... args) {
            log_warning(fmt, std::forward<Args>(args)...);
        }
        inline void printf_info(const char* fmt) {
            log_info(fmt);
        }
        template <typename ...Args>
        inline void printf_info(const char* fmt, Args... args) {
            log_info(fmt, std::forward<Args>(args)...);
        }
        template <typename ...Args>
        inline void printf_error(const char*, int, const char* fmt, Args... args) {
            log_error(fmt, std::forward<Args>(args)...);
        }
        inline void printf(const char* fmt) {
            log_info(fmt);
        }
        template <typename ...Args>
        inline void printf(const char* fmt, Args... args) {
            log_info(fmt, std::forward<Args>(args)...);
        }

        int irand(int imax) { return npnr_ctx->rng(imax+1); }
    }
    
    #define OPEN -1

    #include "vpr_types.h"
    #include "vpr_timing_place.cpp.inc"
    #include "vpr_place.cpp.inc"
}

NEXTPNR_NAMESPACE_BEGIN

class VPRPlacer
{
  public:
    VPRPlacer(Context *ctx) : ctx(ctx)
    {
        vpr::npnr_ctx = ctx;
        int max_y = 0;
        auto &grid = vpr::g_vpr_ctx.device.grid;
        for (auto bel : ctx->getBels()) {
            auto loc = ctx->getBelLocation(bel);
            if (loc.x >= int(grid._bels.size()))
                grid._bels.resize(loc.x+1);
            max_y = std::max(loc.y, max_y);
            if (max_y >= int(grid._bels[loc.x].size()))
                grid._bels[loc.x].resize(max_y+1);
            if (loc.z >= int(grid._bels[loc.x][loc.y].size()))
                grid._bels[loc.x][loc.y].resize(loc.z+1);
            grid._bels[loc.x][loc.y][loc.z] = bel;
        }
        for (auto& c : grid._bels)
            c.resize(max_y+1);

        find_carries();

        int32_t cell_idx = 0;
        for (auto &cell : ctx->cells) {
            CellInfo *ci = cell.second.get();
            if (ci->bel == BelId()) {
                vpr::npnr_cells.push_back(cell.second.get());
            }
            ci->udata = cell_idx++;

            auto loc = ci->attrs.find(ctx->id("BEL"));
            if (loc != ci->attrs.end()) {
                const std::string& loc_name = loc->second;
                auto bel = ctx->getBelByName(ctx->id(loc_name));
                if (bel == BelId()) {
                    log_error("No Bel named \'%s\' located for "
                            "this chip (processing BEL attribute on \'%s\')\n",
                            loc_name.c_str(), ci->name.c_str(ctx));
                }

                auto bel_type = ctx->getBelType(bel);
                if (bel_type != ctx->belTypeFromId(ci->type)) {
                    log_error("Bel \'%s\' of type \'%s\' does not match cell "
                            "\'%s\' of type \'%s\'",
                            loc_name.c_str(), ctx->belTypeToId(bel_type).c_str(ctx), ci->name.c_str(ctx),
                            ci->type.c_str(ctx));
                }
                ctx->bindBel(bel, ci->name, STRENGTH_USER);
            }
        }
        int32_t net_idx = 0;
        for (auto &net : ctx->nets) {
            NetInfo *ni = net.second.get();
            ni->udata = net_idx++;
        }

        vpr::placer_opts.enable_timing_computations = ctx->timing_driven;
    }

    bool place()
    {
        log_break();

        vpr::try_place(vpr::placer_opts, vpr::annealing_sched);

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

    // Generic chain finder
    template <typename F1, typename F2, typename F3>
    std::vector<CellChain> find_chains(const Context *ctx, F1 cell_type_predicate, F2 get_previous, F3 get_next,
                                       size_t min_length = 2)
    {
        std::set<IdString> chained;
        std::vector<CellChain> chains;
        for (auto cell : sorted(ctx->cells)) {
            if (chained.find(cell.first) != chained.end())
                continue;
            CellInfo *ci = cell.second;
            if (cell_type_predicate(ctx, ci)) {
                CellInfo *start = ci;
                CellInfo *prev_start = ci;
                while (prev_start != nullptr) {
                    start = prev_start;
                    prev_start = get_previous(ctx, start);
                }
                CellChain chain;
                CellInfo *end = start;
                while (end != nullptr) {
                    chain.cells.push_back(end);
                    end = get_next(ctx, end);
                }
                if (chain.cells.size() >= min_length) {
                    chains.push_back(chain);
                    for (auto c : chain.cells)
                        chained.insert(c->name);
                }
            }
        }
        return chains;
    }

    // Insert a logic cell to legalise a CIN->fabric connection
    CellInfo *make_carry_feed_in(CellInfo *cin_cell, PortInfo &cin_port)
    {
        NPNR_ASSERT(cin_port.net != nullptr);
        std::unique_ptr<CellInfo> lc = create_ice_cell(ctx, ctx->id("ICESTORM_LC"));
        lc->params[ctx->id("CARRY_ENABLE")] = "1";
        lc->params[ctx->id("CIN_CONST")] = "1";
        lc->params[ctx->id("CIN_SET")] = "1";
        lc->ports.at(ctx->id("I1")).net = cin_port.net;
        cin_port.net->users.erase(std::remove_if(cin_port.net->users.begin(), cin_port.net->users.end(),
                                                 [cin_cell, cin_port](const PortRef &usr) {
                                                     return usr.cell == cin_cell && usr.port == cin_port.name;
                                                 }));

        PortRef i1_ref;
        i1_ref.cell = lc.get();
        i1_ref.port = ctx->id("I1");
        lc->ports.at(ctx->id("I1")).net->users.push_back(i1_ref);

        std::unique_ptr<NetInfo> out_net(new NetInfo());
        out_net->name = ctx->id(lc->name.str(ctx) + "$O");

        PortRef drv_ref;
        drv_ref.port = ctx->id("COUT");
        drv_ref.cell = lc.get();
        out_net->driver = drv_ref;
        lc->ports.at(ctx->id("COUT")).net = out_net.get();

        PortRef usr_ref;
        usr_ref.port = cin_port.name;
        usr_ref.cell = cin_cell;
        out_net->users.push_back(usr_ref);
        cin_cell->ports.at(cin_port.name).net = out_net.get();

        IdString out_net_name = out_net->name;
        NPNR_ASSERT(ctx->nets.find(out_net_name) == ctx->nets.end());
        ctx->nets[out_net_name] = std::move(out_net);

        IdString name = lc->name;
        ctx->assignCellInfo(lc.get());
        ctx->cells[lc->name] = std::move(lc);
        //createdCells.insert(name);
        return ctx->cells[name].get();
    }

    // Insert a logic cell to legalise a COUT->fabric connection
    CellInfo *make_carry_pass_out(PortInfo &cout_port)
    {
        NPNR_ASSERT(cout_port.net != nullptr);
        std::unique_ptr<CellInfo> lc = create_ice_cell(ctx, ctx->id("ICESTORM_LC"));
        lc->params[ctx->id("LUT_INIT")] = "65280"; // 0xff00: O = I3
        lc->params[ctx->id("CARRY_ENABLE")] = "1";
        lc->ports.at(ctx->id("O")).net = cout_port.net;
        std::unique_ptr<NetInfo> co_i3_net(new NetInfo());
        co_i3_net->name = ctx->id(lc->name.str(ctx) + "$I3");
        co_i3_net->driver = cout_port.net->driver;
        PortRef i3_r;
        i3_r.port = ctx->id("I3");
        i3_r.cell = lc.get();
        co_i3_net->users.push_back(i3_r);
        PortRef o_r;
        o_r.port = ctx->id("O");
        o_r.cell = lc.get();
        cout_port.net->driver = o_r;
        lc->ports.at(ctx->id("I3")).net = co_i3_net.get();
        cout_port.net = co_i3_net.get();

        IdString co_i3_name = co_i3_net->name;
        NPNR_ASSERT(ctx->nets.find(co_i3_name) == ctx->nets.end());
        ctx->nets[co_i3_name] = std::move(co_i3_net);
        IdString name = lc->name;
        ctx->assignCellInfo(lc.get());
        ctx->cells[lc->name] = std::move(lc);
        //createdCells.insert(name);
        return ctx->cells[name].get();
    }

    // Split a carry chain into multiple legal chains
    std::vector<CellChain> split_carry_chain(CellChain &carryc)
    {
        bool start_of_chain = true;
        std::vector<CellChain> chains;
        std::vector<const CellInfo *> tile;
        const int max_length = (ctx->chip_info->height - 2) * 8 - 2;
        auto curr_cell = carryc.cells.begin();
        while (curr_cell != carryc.cells.end()) {
            CellInfo *cell = *curr_cell;
            if (tile.size() >= 8) {
                tile.clear();
            }
            if (start_of_chain) {
                tile.clear();
                chains.emplace_back();
                start_of_chain = false;
                if (cell->ports.at(ctx->id("CIN")).net) {
                    // CIN is not constant and not part of a chain. Must feed in from fabric
                    CellInfo *feedin = make_carry_feed_in(cell, cell->ports.at(ctx->id("CIN")));
                    chains.back().cells.push_back(feedin);
                    tile.push_back(feedin);
                }
            }
            tile.push_back(cell);
            chains.back().cells.push_back(cell);
            bool split_chain = (!ctx->logicCellsCompatible(tile)) || (int(chains.back().cells.size()) > max_length);
            if (split_chain) {
                CellInfo *passout = make_carry_pass_out(cell->ports.at(ctx->id("COUT")));
                tile.pop_back();
                chains.back().cells.back() = passout;
                start_of_chain = true;
            } else {
                NetInfo *carry_net = cell->ports.at(ctx->id("COUT")).net;
                bool at_end = (curr_cell == carryc.cells.end() - 1);
                if (carry_net != nullptr && (carry_net->users.size() > 1 || at_end)) {
                    if (carry_net->users.size() > 2 ||
                        (net_only_drives(ctx, carry_net, is_lc, ctx->id("I3"), false) !=
                         net_only_drives(ctx, carry_net, is_lc, ctx->id("CIN"), false)) ||
                        (at_end && !net_only_drives(ctx, carry_net, is_lc, ctx->id("I3"), true))) {
                        CellInfo *passout = make_carry_pass_out(cell->ports.at(ctx->id("COUT")));
                        chains.back().cells.push_back(passout);
                        tile.push_back(passout);
                        start_of_chain = true;
                    }
                }
                ++curr_cell;
            }
        }
        return chains;
    }

    bool find_carries()
    {
        std::vector<CellChain> carry_chains =
                find_chains(ctx, [](const Context *ctx, const CellInfo *cell) { return is_lc(ctx, cell); },
                            [](const Context *ctx, const

                               CellInfo *cell) {
                                CellInfo *carry_prev =
                                        net_driven_by(ctx, cell->ports.at(ctx->id("CIN")).net, is_lc, ctx->id("COUT"));
                                if (carry_prev != nullptr)
                                    return carry_prev;
                                /*CellInfo *i3_prev = net_driven_by(ctx, cell->ports.at(ctx->id("I3")).net, is_lc,
                                ctx->id("COUT")); if (i3_prev != nullptr) return i3_prev;*/
                                return (CellInfo *)nullptr;
                            },
                            [](const Context *ctx, const CellInfo *cell) {
                                CellInfo *carry_next = net_only_drives(ctx, cell->ports.at(ctx->id("COUT")).net, is_lc,
                                                                       ctx->id("CIN"), false);
                                if (carry_next != nullptr)
                                    return carry_next;
                                /*CellInfo *i3_next =
                                        net_only_drives(ctx, cell->ports.at(ctx->id("COUT")).net, is_lc, ctx->id("I3"),
                                false); if (i3_next != nullptr) return i3_next;*/
                                return (CellInfo *)nullptr;
                            });
        std::unordered_set<IdString> chained;
        for (auto &base_chain : carry_chains) {
            for (auto c : base_chain.cells)
                chained.insert(c->name);
        }
        // Any cells not in chains, but with carry enabled, must also be put in a single-carry chain
        // for correct processing
        for (auto cell : sorted(ctx->cells)) {
            CellInfo *ci = cell.second;
            if (chained.find(cell.first) == chained.end() && is_lc(ctx, ci) &&
                bool_or_default(ci->params, ctx->id("CARRY_ENABLE"))) {
                CellChain sChain;
                sChain.cells.push_back(ci);
                chained.insert(cell.first);
                carry_chains.push_back(sChain);
            }
        }

        // Find midpoints for all chains, before we start tearing them up
        for (auto &base_chain : carry_chains) {
            /*if (ctx->verbose)*/ {
                log_info("Found carry chain: \n");
                for (auto entry : base_chain.cells)
                    log_info("     %s\n", entry->name.c_str(ctx));
                log_info("\n");
            }
            std::vector<CellChain> split_chains = split_carry_chain(base_chain);
            for (auto &chain : split_chains) {
                //get_chain_midpoint(ctx, chain, chain.mid_x, chain.mid_y);
                all_chains.push_back(chain);
            }
        }
    }

  private:
    Context *ctx;
};

bool placer_vpr(Context *ctx)
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
