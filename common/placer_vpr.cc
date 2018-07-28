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
#include "place_legaliser.h"

namespace NEXTPNR_NAMESPACE {
    std::vector<std::vector<CellInfo*>> carries;
}

namespace vpr {
    using namespace NEXTPNR_NAMESPACE;

    static Context* npnr_ctx = NULL;
    std::vector<CellInfo *> npnr_cells;

    struct DeviceGrid {
        inline size_t width() const { return _bels.size(); }
        inline size_t height() const { return _bels.front().size(); }
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
        delay_t worst_slack;
        delay_t max_req;
        void update() { 
#if 1
            update_budget(npnr_ctx);

            worst_slack = std::numeric_limits<decltype(worst_slack)>::max();

            // Compute the delay for every pin on every net
            for (auto &n : npnr_ctx->nets) {
                auto net = n.second.get();

                bool driver_gb;
                CellInfo *driver_cell = net->driver.cell;
                if (!driver_cell)
                    continue;
                if (driver_cell->bel == BelId())
                    continue;
                driver_gb = npnr_ctx->getBelGlobalBuf(driver_cell->bel);
                WireId drv_wire = npnr_ctx->getBelPinWire(driver_cell->bel, npnr_ctx->portPinFromId(net->driver.port));
                if (driver_gb)
                    continue;
                for (auto& load : net->users) {
                    if (load.cell == nullptr)
                        continue;
                    CellInfo *load_cell = load.cell;
                    if (load_cell->bel == BelId())
                        continue;
                    WireId user_wire = npnr_ctx->getBelPinWire(load_cell->bel, npnr_ctx->portPinFromId(load.port));
                    delay_t raw_wl = npnr_ctx->estimateDelay(drv_wire, user_wire);
                    delay_t slack = load.budget - raw_wl;
                    worst_slack = std::min(worst_slack, slack);
                }
            }

            max_req = delay_t(1.0e12 / npnr_ctx->target_freq);
#endif
        }
    };
    static SetupTimingInfo timing_info;
    struct t_pl_macro_member {
        t_pl_macro_member(CellInfo* blk_index, int x_offset, int y_offset, int z_offset) : blk_index(blk_index), x_offset(x_offset), y_offset(y_offset), z_offset(z_offset) {}
        CellInfo* blk_index;
        int x_offset;
        int y_offset;
        int z_offset;
    };
    struct t_pl_macro {
        std::vector<t_pl_macro_member> members;
    };

    // timing_util.cpp
    float calculate_clb_net_pin_criticality(const SetupTimingInfo& timing_info, /*const ClusteredPinAtomPinsLookup& pin_lookup,*/ const PortRef& load, const NetInfo* net)
    {
        NPNR_ASSERT(npnr_ctx->timing_driven);

#if 1
        bool driver_gb;
        CellInfo *driver_cell = net->driver.cell;
        if (!driver_cell)
            return 0;
        if (driver_cell->bel == BelId())
            return 0;
        driver_gb = npnr_ctx->getBelGlobalBuf(driver_cell->bel);
        WireId drv_wire = npnr_ctx->getBelPinWire(driver_cell->bel, npnr_ctx->portPinFromId(net->driver.port));
        if (driver_gb)
            return 0;
        if (load.cell == nullptr)
            return 0;
        CellInfo *load_cell = load.cell;
        if (load_cell->bel == BelId())
            return 0;
        WireId user_wire = npnr_ctx->getBelPinWire(load_cell->bel, npnr_ctx->portPinFromId(load.port));
        delay_t raw_wl = npnr_ctx->estimateDelay(drv_wire, user_wire);
        delay_t slack = load.budget - raw_wl;
        delay_t shift = timing_info.worst_slack < 0 ? -timing_info.worst_slack : 0;
        float crit = 1 - (float(slack + shift) / (timing_info.max_req + shift));
        crit = std::max<float>(0., crit);
        crit = std::min<float>(1., crit);
        return crit;
#else
        return 1;
#endif
    }

    // place_macro.cpp
    int alloc_and_load_placement_macros(/*t_direct_inf* directs, int num_directs, t_pl_macro ** */ std::unordered_map<CellInfo*, t_pl_macro> &macros)
    {
        for (auto &chain : carries) {
            t_pl_macro entry;
            auto head = chain.front();
            for (int z = 0; z < int(chain.size()); z++) {
                auto cell = chain.at(z);
                entry.members.emplace_back(cell, 0, z / 8, z % 8);
                cell->attrs.emplace(npnr_ctx->id("carry_head"), head->name.str(npnr_ctx));
            }
            macros.emplace(head, std::move(entry));
        }
        return macros.size();
    }

    void get_imacro_from_iblk(/*int **/ CellInfo **imacro, /*ClusterBlockId*/ CellInfo *iblk /*, t_pl_macro *macros, int num_macros*/) {
        auto it = iblk->attrs.find(npnr_ctx->id("carry_head"));
        if (it != iblk->attrs.end())
            *imacro = npnr_ctx->cells.at(npnr_ctx->id(it->second)).get();
        else
            *imacro = nullptr;
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

        carries = prepare_and_find_carries(ctx);

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
        compute_fmax(ctx, true /* print_fmax */);
        return true;
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
