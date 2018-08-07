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
#include "timing.h"
#include "util.h"

namespace vpr {
    using namespace NEXTPNR_NAMESPACE;

    static Context* npnr_ctx = NULL;
    std::vector<CellInfo *> npnr_cells;

    // base/device_grid.h
    struct DeviceGrid {
        inline size_t width() const { return npnr_ctx->chip_info->width; }
        inline size_t height() const { return npnr_ctx->chip_info->height; }
        inline const std::vector<std::vector<BelId>>& operator[](size_t x) { return _bels.at(x); }
        std::vector<std::vector<std::vector<BelId>>> _bels;
    };
    // base/netlist_fwd.h
    enum PinType
    {
        DRIVER = PortType::PORT_OUT,
        SINK = PortType::PORT_IN,
    };
    // base/clustered_netlist_fwd.h
    typedef size_t ClusterBlockId;
    // base/globals.h
    static struct VprContext {
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
    // libtatum
    namespace tatum {
        struct TimingPathInfo {
            TimingPathInfo(float delay=0) : _delay(delay) {}
            float delay() { return _delay; }
            float _delay;
        };
    };
    // timing/timing_info.h
    struct SetupTimingInfo {
        delay_t sWNS, sTNS;
        delay_t max_req;
        delay_t worst_path_slack;
        void update() { 
            if (npnr_ctx->slack_redist_iter > 0)
                worst_path_slack = assign_budget(npnr_ctx, true /* quiet */);
            else
                worst_path_slack = timing_analysis(npnr_ctx, false /* print_fmax */, false /* print_histogram */);

            sWNS = std::numeric_limits<decltype(sWNS)>::max();
            sTNS = 0;
            max_req = delay_t(1.0e12 / npnr_ctx->target_freq);

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
                if (driver_gb)
                    continue;
                for (auto& load : net->users) {
                    if (!load.cell)
                        continue;
                    CellInfo *load_cell = load.cell;
                    if (load_cell->bel == BelId())
                        continue;
                    auto net_delay = npnr_ctx->getNetinfoRouteDelay(net, load);
                    auto slack = load.budget - net_delay;
                    sWNS = std::min(sWNS, slack);
                    if (slack < 0)
                        sTNS += slack;
                }
            }
        }
        tatum::TimingPathInfo least_slack_critical_path()
        {
            return tatum::TimingPathInfo(npnr_ctx->getDelayNS(max_req - worst_path_slack));
        }
        float setup_total_negative_slack() { return npnr_ctx->getDelayNS(sTNS); }
        float setup_worst_negative_slack() { return npnr_ctx->getDelayNS(sWNS); }
    };
    std::unique_ptr<SetupTimingInfo> make_setup_timing_info(/*std::shared_ptr<DelayCalc> delay_calculator*/) {
        return std::unique_ptr<SetupTimingInfo>(new SetupTimingInfo);
    }
    // timing/timing_util.h
    float calculate_clb_net_pin_criticality(const SetupTimingInfo& timing_info, /*const ClusteredPinAtomPinsLookup& pin_lookup,*/ const PortRef& load, const NetInfo* net)
    {
        NPNR_ASSERT(npnr_ctx->timing_driven);

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
        delay_t shift = std::min(timing_info.sWNS, 0);
        float crit = 1 - (float(slack + shift) / (timing_info.max_req + shift));
        crit = std::max<float>(0., crit);
        crit = std::min<float>(1., crit);
        return crit;
    }
    // draw/draw.h
    void update_screen(/*ScreenUpdatePriority priority, const char *msg, enum pic_type pic_on_screen_val,
                    std::shared_ptr<SetupTimingInfo> timing_info*/)
    {
        npnr_ctx->yield();
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
        float frand() { return npnr_ctx->rng() / float(0x3fffffff); }
    }
    
    // libarchfpga
    #define OPEN -1

    #include "vpr/place/timing_place.cpp"
    #include "vpr/place/place.cpp"
    #include "vpr/place/place_macro.cpp"
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
                ctx->bindBel(bel, ci, STRENGTH_USER);
            }
        }
        int32_t net_idx = 0;
        for (auto &net : ctx->nets) {
            NetInfo *ni = net.second.get();
            ni->udata = net_idx++;
        }
    }

    bool place()
    {
        log_break();
        ctx->lock();

        vpr::t_placer_opts placer_opts;
        placer_opts.place_algorithm = vpr::PATH_TIMING_DRIVEN_PLACE;
        placer_opts.enable_timing_computations = ctx->timing_driven;
        placer_opts.inner_loop_recompute_divider = 0;
        placer_opts.recompute_crit_iter = 1;
        placer_opts.td_place_exp_first = 1.0;
        placer_opts.td_place_exp_last = 8.0;
        placer_opts.timing_tradeoff = 0.5;

        vpr::t_annealing_sched annealing_sched;
        annealing_sched.type = vpr::AUTO_SCHED;
        annealing_sched.inner_num = 10;

        vpr::try_place(placer_opts, annealing_sched);

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
