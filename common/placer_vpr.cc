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
    }

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
            int x, y;
            bool gb;
            ctx->estimatePosition(bel, x, y, gb);
            if (x >= int(grid._bels.size()))
                grid._bels.resize(x+1);
            max_y = std::max(y, max_y);
            if (max_y >= int(grid._bels[x].size()))
                grid._bels[x].resize(max_y+1);
            grid._bels[x][y].push_back(bel);
        }
        for (auto& c : grid._bels)
            c.resize(max_y+1);

        int32_t cell_idx = 0;
        for (auto &cell : ctx->cells) {
            CellInfo *ci = cell.second.get();
            if (ci->bel == BelId()) {
                vpr::npnr_cells.push_back(cell.second.get());
            }
            ci->udata = cell_idx++;
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
