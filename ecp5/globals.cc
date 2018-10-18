/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  David Shah <david@symbioticeda.com>
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

#include "globals.h"
#include <algorithm>
#include <iomanip>
#include <queue>
#include "cells.h"
#include "log.h"
#include "nextpnr.h"
#include "place_common.h"
#define fmt_str(x) (static_cast<const std::ostringstream &>(std::ostringstream() << x).str())

NEXTPNR_NAMESPACE_BEGIN

static std::string get_quad_name(GlobalQuadrant quad)
{
    switch (quad) {
    case QUAD_UL:
        return "UL";
    case QUAD_UR:
        return "UR";
    case QUAD_LL:
        return "LL";
    case QUAD_LR:
        return "LR";
    }
    return "";
}

class Ecp5GlobalRouter
{
  public:
    Ecp5GlobalRouter(Context *ctx) : ctx(ctx){};

  private:
    bool is_clock_port(const PortRef &user)
    {
        if (user.cell->type == id_TRELLIS_SLICE && (user.port == id_CLK || user.port == id_WCK))
            return true;
        return false;
    }

    std::vector<NetInfo *> get_clocks()
    {
        std::unordered_map<IdString, int> clockCount;
        for (auto &net : ctx->nets) {
            NetInfo *ni = net.second.get();
            clockCount[ni->name] = 0;
            for (const auto &user : ni->users) {
                if (is_clock_port(user))
                    clockCount[ni->name]++;
            }
            // log_info("clkcount %s: %d\n", ni->name.c_str(ctx),clockCount[ni->name]);
        }
        std::vector<NetInfo *> clocks;
        while (clocks.size() < 16) {
            auto max = std::max_element(clockCount.begin(), clockCount.end(),
                                        [](const decltype(clockCount)::value_type &a,
                                           const decltype(clockCount)::value_type &b) { return a.second < b.second; });
            if (max == clockCount.end() || max->second < 5)
                break;
            clocks.push_back(ctx->nets.at(max->first).get());
            clockCount.erase(max->first);
        }
        return clocks;
    }

    PipId find_tap_pip(WireId tile_glb)
    {
        std::string wireName = ctx->getWireBasename(tile_glb).str(ctx);
        std::string glbName = wireName.substr(2);
        TapDirection td = ctx->globalInfoAtLoc(tile_glb.location).tap_dir;
        WireId tap_wire;
        Location tap_loc;
        tap_loc.x = ctx->globalInfoAtLoc(tile_glb.location).tap_col;
        tap_loc.y = tile_glb.location.y;
        if (td == TAP_DIR_LEFT) {
            tap_wire = ctx->getWireByLocAndBasename(tap_loc, "L_" + glbName);
        } else {
            tap_wire = ctx->getWireByLocAndBasename(tap_loc, "R_" + glbName);
        }
        NPNR_ASSERT(tap_wire != WireId());
        return *(ctx->getPipsUphill(tap_wire).begin());
    }

    PipId find_spine_pip(WireId tap_wire)
    {
        std::string wireName = ctx->getWireBasename(tap_wire).str(ctx);
        Location spine_loc;
        spine_loc.x = ctx->globalInfoAtLoc(tap_wire.location).spine_col;
        spine_loc.y = ctx->globalInfoAtLoc(tap_wire.location).spine_row;
        WireId spine_wire = ctx->getWireByLocAndBasename(spine_loc, wireName);
        return *(ctx->getPipsUphill(spine_wire).begin());
    }

    void route_logic_tile_global(NetInfo *net, int global_index, PortRef user)
    {
        WireId userWire = ctx->getBelPinWire(user.cell->bel, user.port);
        WireId globalWire;
        IdString global_name = ctx->id(fmt_str("G_HPBX" << std::setw(2) << std::setfill('0') << global_index << "00"));
        std::queue<WireId> upstream;
        std::unordered_map<WireId, PipId> backtrace;
        upstream.push(userWire);
        bool already_routed = false;
        WireId next;
        // Search back from the pin until we reach the global network
        while (true) {
            next = upstream.front();
            upstream.pop();

            if (ctx->getBoundWireNet(next) == net) {
                already_routed = true;
                globalWire = next;
                break;
            }

            if (ctx->getWireBasename(next) == global_name) {
                globalWire = next;
                break;
            }
            if (ctx->checkWireAvail(next)) {
                for (auto pip : ctx->getPipsUphill(next)) {
                    WireId src = ctx->getPipSrcWire(pip);
                    backtrace[src] = pip;
                    upstream.push(src);
                }
            }
            if (upstream.size() > 30000) {
                log_error("failed to route HPBX%02d00 to %s.%s\n", global_index,
                          ctx->getBelName(user.cell->bel).c_str(ctx), user.port.c_str(ctx));
            }
        }
        // Set all the pips we found along the way
        WireId cursor = next;
        while (true) {
            auto fnd = backtrace.find(cursor);
            if (fnd == backtrace.end())
                break;
            ctx->bindPip(fnd->second, net, STRENGTH_LOCKED);
            cursor = ctx->getPipDstWire(fnd->second);
        }
        // If the global network inside the tile isn't already set up,
        // we also need to bind the buffers along the way
        if (!already_routed) {
            ctx->bindWire(next, net, STRENGTH_LOCKED);
            PipId tap_pip = find_tap_pip(next);
            NetInfo *tap_net = ctx->getBoundPipNet(tap_pip);
            if (tap_net == nullptr) {
                ctx->bindPip(tap_pip, net, STRENGTH_LOCKED);
                PipId spine_pip = find_spine_pip(ctx->getPipSrcWire(tap_pip));
                NetInfo *spine_net = ctx->getBoundPipNet(spine_pip);
                if (spine_net == nullptr) {
                    ctx->bindPip(spine_pip, net, STRENGTH_LOCKED);
                } else {
                    NPNR_ASSERT(spine_net == net);
                }
            } else {
                NPNR_ASSERT(tap_net == net);
            }
        }
    }

    bool is_global_io(CellInfo *io, std::string &glb_name)
    {
        std::string func_name = ctx->getPioFunctionName(io->bel);
        if (func_name.substr(0, 5) == "PCLKT") {
            func_name.erase(func_name.find('_'), 1);
            glb_name = "G_" + func_name;
            return true;
        }
        return false;
    }

    WireId get_global_wire(GlobalQuadrant quad, int network)
    {
        return ctx->getWireByLocAndBasename(Location(0, 0),
                                            "G_" + get_quad_name(quad) + "PCLK" + std::to_string(network));
    }

    bool simple_router(NetInfo *net, WireId src, WireId dst, bool allow_fail = false)
    {
        std::queue<WireId> visit;
        std::unordered_map<WireId, PipId> backtrace;
        visit.push(src);
        WireId cursor;
        while (true) {

            if (visit.empty() || visit.size() > 50000) {
                if (allow_fail)
                    return false;
                log_error("cannot route global from %s to %s.\n", ctx->getWireName(src).c_str(ctx),
                          ctx->getWireName(dst).c_str(ctx));
            }
            cursor = visit.front();
            visit.pop();
            NetInfo *bound = ctx->getBoundWireNet(cursor);
            if (bound == net) {
            } else if (bound != nullptr) {
                continue;
            }
            if (cursor == dst)
                break;
            for (auto dh : ctx->getPipsDownhill(cursor)) {
                WireId pipDst = ctx->getPipDstWire(dh);
                if (backtrace.count(pipDst))
                    continue;
                backtrace[pipDst] = dh;
                visit.push(pipDst);
            }
        }
        while (true) {
            auto fnd = backtrace.find(cursor);
            if (fnd == backtrace.end())
                break;
            NetInfo *bound = ctx->getBoundWireNet(cursor);
            if (bound != nullptr) {
                NPNR_ASSERT(bound == net);
                break;
            }
            ctx->bindPip(fnd->second, net, STRENGTH_LOCKED);
            cursor = ctx->getPipSrcWire(fnd->second);
        }
        if (ctx->getBoundWireNet(src) == nullptr)
            ctx->bindWire(src, net, STRENGTH_LOCKED);
        return true;
    }

    bool route_onto_global(NetInfo *net, int network)
    {
        WireId glb_src;
        NPNR_ASSERT(net->driver.cell->type == id_DCCA);
        glb_src = ctx->getNetinfoSourceWire(net);
        for (int quad = QUAD_UL; quad < QUAD_LR + 1; quad++) {
            WireId glb_dst = get_global_wire(GlobalQuadrant(quad), network);
            NPNR_ASSERT(glb_dst != WireId());
            bool routed = simple_router(net, glb_src, glb_dst);
            if (!routed)
                return false;
        }
        return true;
    }

    // Attempt to place a DCC
    void place_dcc(CellInfo *dcc)
    {
        BelId best_bel;
        wirelen_t best_wirelen = 9999999;
        for (auto bel : ctx->getBels()) {
            if (ctx->getBelType(bel) == id_DCCA && ctx->checkBelAvail(bel)) {
                if (ctx->isValidBelForCell(dcc, bel)) {
                    ctx->bindBel(bel, dcc, STRENGTH_LOCKED);
                    float tns;
                    wirelen_t wirelen = get_net_metric(ctx, dcc->ports.at(id_CLKI).net, MetricType::WIRELENGTH, tns);
                    if (wirelen < best_wirelen) {
                        best_bel = bel;
                        best_wirelen = wirelen;
                    }
                    ctx->unbindBel(bel);
                }
            }
        }
        NPNR_ASSERT(best_bel != BelId());
        ctx->bindBel(best_bel, dcc, STRENGTH_LOCKED);
    }

    // Insert a DCC into a net to promote it to a global
    NetInfo *insert_dcc(NetInfo *net)
    {
        auto dcc = create_ecp5_cell(ctx, id_DCCA, "$gbuf$" + net->name.str(ctx));

        std::unique_ptr<NetInfo> glbnet = std::unique_ptr<NetInfo>(new NetInfo);
        glbnet->name = ctx->id("$glbnet$" + net->name.str(ctx));
        glbnet->driver.cell = dcc.get();
        glbnet->driver.port = id_CLKO;
        glbnet->is_global = true;
        dcc->ports[id_CLKO].net = glbnet.get();

        glbnet->users = net->users;
        for (auto user : net->users) {
            user.cell->ports.at(user.port).net = glbnet.get();
        }
        net->users.clear();

        dcc->ports[id_CLKI].net = net;
        PortRef clki_pr;
        clki_pr.port = id_CLKI;
        clki_pr.cell = dcc.get();
        net->users.push_back(clki_pr);

        place_dcc(dcc.get());

        ctx->cells[dcc->name] = std::move(dcc);
        NetInfo *glbptr = glbnet.get();
        ctx->nets[glbnet->name] = std::move(glbnet);
        return glbptr;
    }

    int global_route_priority(const PortRef &load)
    {
        if (load.port == id_WCK || load.port == id_WRE)
            return 90;
        return 99;
    }

    Context *ctx;

  public:
    void promote_and_route_globals()
    {
        log_info("Promoting and routing globals...\n");
        auto clocks = get_clocks();
        std::set<int> all_globals, fab_globals;
        for (int i = 0; i < 16; i++) {
            all_globals.insert(i);
            if (i < 8)
                fab_globals.insert(i);
        }
        for (auto clock : clocks) {
            bool drives_fabric = std::any_of(clock->users.begin(), clock->users.end(),
                                             [this](const PortRef &port) { return !is_clock_port(port); });
            int glbid;
            if (drives_fabric) {
                if (fab_globals.empty())
                    continue;
                glbid = *(fab_globals.begin());
            } else {
                glbid = *(all_globals.begin());
            }
            all_globals.erase(glbid);
            fab_globals.erase(glbid);
            log_info("    promoting clock net %s to global %d\n", clock->name.c_str(ctx), glbid);
            auto old_users = clock->users;
            NetInfo *global = insert_dcc(clock);
            bool routed = route_onto_global(global, glbid);
            NPNR_ASSERT(routed);

            // WCK must have routing priority
            auto sorted_users = global->users;
            std::sort(sorted_users.begin(), sorted_users.end(), [this](const PortRef &a, const PortRef &b) {
                return global_route_priority(a) < global_route_priority(b);
            });
            for (const auto &user : sorted_users) {
                route_logic_tile_global(global, glbid, user);
            }
        }
    }
};

void route_ecp5_globals(Context *ctx) { Ecp5GlobalRouter(ctx).promote_and_route_globals(); }

NEXTPNR_NAMESPACE_END
