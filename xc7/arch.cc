/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  Clifford Wolf <clifford@symbioticeda.com>
 *  Copyright (C) 2018  Serge Bazanski <q3k@symbioticeda.com>
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

#include <algorithm>
#include <cmath>
#include "cells.h"
#include "gfx.h"
#include "log.h"
#include "nextpnr.h"
#include "placer1.h"
#include "router1.h"
#include "util.h"

#include "torc/common/DirectoryTree.hpp"

NEXTPNR_NAMESPACE_BEGIN

std::unique_ptr<const TorcInfo> torc_info;
TorcInfo::TorcInfo(Arch *ctx, const std::string &inDeviceName, const std::string &inPackageName)
        : ddb(new DDB(inDeviceName, inPackageName)), sites(ddb->getSites()), tiles(ddb->getTiles()),
          segments(ddb->getSegments()), bel_to_site_index(construct_bel_to_site_index(ctx, sites)),
          num_bels(bel_to_site_index.size()),
          site_index_to_bel(construct_site_index_to_bel(ctx, sites, bel_to_site_index)),
          site_index_to_type(construct_site_index_to_type(ctx, sites)),
          bel_to_loc(construct_bel_to_loc(sites, tiles, num_bels, site_index_to_type)),
          wire_to_tilewire(construct_wire_to_tilewire(segments, tiles, segment_to_wire, trivial_to_wire)),
          num_wires(wire_to_tilewire.size()), wire_to_delay(construct_wire_to_delay(tiles, wire_to_tilewire, *ddb)),
          pip_to_arc(construct_pip_to_arc(wire_to_tilewire, *ddb, wire_to_pips_uphill, wire_to_pips_downhill)),
          num_pips(pip_to_arc.size())
{
    pip_to_dst_wire.reserve(num_pips);

    for (const auto &arc : pip_to_arc) {
        const auto &tw = arc.getSinkTilewire();
        pip_to_dst_wire.push_back(tilewire_to_wire(tw));
    }
}
std::vector<SiteIndex> TorcInfo::construct_bel_to_site_index(Arch *ctx, const Sites &sites)
{
    std::vector<SiteIndex> bel_to_site_index;
    bel_to_site_index.reserve(sites.getSiteCount());
    for (SiteIndex i(0); i < sites.getSiteCount(); ++i) {
        const auto &s = sites.getSite(i);
        const auto &pd = s.getPrimitiveDefPtr();
        const auto &type = pd->getName();
        if (type == "SLICEL" || type == "SLICEM") {
            bel_to_site_index.push_back(i);
            bel_to_site_index.push_back(i);
            bel_to_site_index.push_back(i);
            bel_to_site_index.push_back(i);
        } else
            bel_to_site_index.push_back(i);
    }
    return bel_to_site_index;
}
std::vector<BelId> TorcInfo::construct_site_index_to_bel(Arch *ctx, const Sites &sites, const std::vector<SiteIndex> &bel_to_site_index)
{
    std::vector<BelId> site_index_to_bel;
    site_index_to_bel.resize(sites.getSiteCount());
    BelId b;
    b.index = 0;
    for (auto i : bel_to_site_index) {
        site_index_to_bel[i] = b;
        ++b.index;
    }
    return site_index_to_bel;
}
std::vector<IdString> TorcInfo::construct_site_index_to_type(Arch *ctx, const Sites &sites)
{
    std::vector<IdString> site_index_to_type;
    site_index_to_type.resize(sites.getSiteCount());
    for (SiteIndex i(0); i < sites.getSiteCount(); ++i) {
        const auto &s = sites.getSite(i);
        const auto &pd = s.getPrimitiveDefPtr();
        const auto &type = pd->getName();
        if (type == "SLICEL" || type == "SLICEM")
            site_index_to_type[i] = id_SLICE_LUT6;
        else if (type == "IOB33S" || type == "IOB33M")
            site_index_to_type[i] = id_IOB33;
        else
            site_index_to_type[i] = ctx->id(type);
    }
    return site_index_to_type;
}
std::vector<Loc> TorcInfo::construct_bel_to_loc(const Sites &sites, const Tiles &tiles, const int num_bels,
                                                const std::vector<IdString> &site_index_to_type)
{
    std::vector<Loc> bel_to_loc;
    bel_to_loc.resize(num_bels);
    int32_t bel_index = 0;
    for (SiteIndex i(0); i < site_index_to_type.size(); ++i) {
        const auto &site = sites.getSite(i);
        const auto &tile_info = tiles.getTileInfo(site.getTileIndex());
        const auto x = (tile_info.getCol() + 1) / 2; // Divide by 2 because XDL coordinate space counts the INT tiles between CLBs
        const auto y = tile_info.getRow();

        if (site_index_to_type[i] == id_SLICE_LUT6) {
            const auto site_name = site.getName();
            const auto site_name_back = site_name.back();
            if (site_name_back == '0' || site_name_back == '2' || site_name_back == '4' || site_name_back == '6' ||
                site_name_back == '8') {
                bel_to_loc[bel_index++] = Loc(x, y, 0);
                bel_to_loc[bel_index++] = Loc(x, y, 1);
                bel_to_loc[bel_index++] = Loc(x, y, 2);
                bel_to_loc[bel_index++] = Loc(x, y, 3);
            } else {
                bel_to_loc[bel_index++] = Loc(x, y, 4);
                bel_to_loc[bel_index++] = Loc(x, y, 5);
                bel_to_loc[bel_index++] = Loc(x, y, 6);
                bel_to_loc[bel_index++] = Loc(x, y, 7);
            }
        } else
            bel_to_loc[bel_index++] = Loc(x, y, 0);
    }
    return bel_to_loc;
}
std::vector<Tilewire>
TorcInfo::construct_wire_to_tilewire(const Segments &segments, const Tiles &tiles,
                                     std::unordered_map<Segments::SegmentReference, int> &segment_to_wire,
                                     std::unordered_map<Tilewire, int> &trivial_to_wire)
{
    std::vector<Tilewire> wire_to_tilewire;

    Tilewire currentTilewire;
    for (TileIndex tileIndex(0); tileIndex < tiles.getTileCount(); tileIndex++) {
        // iterate over every wire in the tile
        const auto &tileInfo = tiles.getTileInfo(tileIndex);
        auto tileTypeIndex = tileInfo.getTypeIndex();
        auto wireCount = tiles.getWireCount(tileTypeIndex);
        currentTilewire.setTileIndex(tileIndex);
        for (WireIndex wireIndex(0); wireIndex < wireCount; wireIndex++) {
            currentTilewire.setWireIndex(wireIndex);
            const auto &currentSegment = segments.getTilewireSegment(currentTilewire);

            if (!currentSegment.isTrivial()) {
                if (currentSegment.getAnchorTileIndex() != tileIndex)
                    continue;
                segment_to_wire.emplace(currentSegment, wire_to_tilewire.size());
            } else
                trivial_to_wire.emplace(currentTilewire, wire_to_tilewire.size());

            wire_to_tilewire.push_back(currentTilewire);
        }
    }

    wire_to_tilewire.shrink_to_fit();
    return wire_to_tilewire;
}
std::vector<DelayInfo> TorcInfo::construct_wire_to_delay(const Tiles &tiles, const std::vector<Tilewire> &wire_to_tilewire, const DDB &ddb)
{
    std::vector<DelayInfo> wire_to_delay;
    wire_to_delay.reserve(wire_to_tilewire.size());

    const boost::regex re_124 = boost::regex("(.+_)?[NESW][NESWLR](\\d)((BEG(_[NS])?)|(END(_[NS])?)|[A-E])?\\d(_\\d)?");
    const boost::regex re_L = boost::regex("(.+_)?L(H|V|VB)(_L)?\\d+(_\\d)?");
    const boost::regex re_BYP = boost::regex("BYP(_ALT)?\\d");
    const boost::regex re_BYP_B = boost::regex("BYP_[BL]\\d");
    const boost::regex re_BOUNCE_NS = boost::regex("(BYP|FAN)_BOUNCE_[NS]3_\\d");
    const boost::regex re_FAN = boost::regex("FAN(_ALT)?\\d");
    const boost::regex re_CLB_I1_6 = boost::regex("CLBL[LM]_(L|LL|M)_[A-D]([1-6])");

    std::unordered_map</*TileTypeIndex*/unsigned, std::vector<delay_t>> delay_lookup;

    boost::cmatch what;
    for (const auto &tw : wire_to_tilewire) {
        const TileInfo& tileInfo = tiles.getTileInfo(tw.getTileIndex());
        auto tile_type_index = tileInfo.getTypeIndex();

        auto it = delay_lookup.find(tile_type_index);
        if (it == delay_lookup.end()) {
            auto wireCount = tiles.getWireCount(tile_type_index);
            std::vector<delay_t> tile_delays(wireCount);
            for (WireIndex wireIndex(0); wireIndex < wireCount; wireIndex++) {
                const WireInfo& wireInfo = tiles.getWireInfo(tile_type_index, wireIndex);
                auto wire_name = wireInfo.getName();
                if (boost::regex_match(wire_name, what, re_124)) {
                    switch (what.str(2)[0]) {
                    case '1': tile_delays[wireIndex] = 150; break;
                    case '2': tile_delays[wireIndex] = 170; break;
                    case '4': tile_delays[wireIndex] = 210; break;
                    case '6': tile_delays[wireIndex] = 210; break;
                    default: throw;
                    }
                } else if (boost::regex_match(wire_name, what, re_L)) {
                    std::string l(what[2]);
                    if (l == "H")
                        tile_delays[wireIndex] = 360;
                    else if (l == "VB")
                        tile_delays[wireIndex] = 300;
                    else if (l == "V")
                        tile_delays[wireIndex] = 350;
                    else
                        throw;
                } else if (boost::regex_match(wire_name, what, re_BYP)) {
                    tile_delays[wireIndex] = 190;
                } else if (boost::regex_match(wire_name, what, re_BYP_B)) {
                } else if (boost::regex_match(wire_name, what, re_FAN)) {
                    tile_delays[wireIndex] = 190;
                } else if (boost::regex_match(wire_name, what, re_CLB_I1_6)) {
                    switch (what.str(2)[0]) {
                        case '1': tile_delays[wireIndex] = 280; break;
                        case '2': tile_delays[wireIndex] = 280; break;
                        case '3': tile_delays[wireIndex] = 180; break;
                        case '4': tile_delays[wireIndex] = 180; break;
                        case '5': tile_delays[wireIndex] =  80; break;
                        case '6': tile_delays[wireIndex] =  40; break;
                        default: throw;
                    }
                }
            }
            it = delay_lookup.emplace(tile_type_index, std::move(tile_delays)).first;
        }
        assert(it != delay_lookup.end());
        DelayInfo d;
        d.delay = it->second[tw.getWireIndex()];
        wire_to_delay.emplace_back(std::move(d));
    }

    return wire_to_delay;
}
std::vector<Arc> TorcInfo::construct_pip_to_arc(const std::vector<Tilewire> &wire_to_tilewire, const DDB &ddb,
                                                std::vector<std::vector<int>> &wire_to_pips_uphill,
                                                std::vector<std::vector<int>> &wire_to_pips_downhill)
{
    const auto &tiles = ddb.getTiles();

    std::vector<Arc> pip_to_arc;
    wire_to_pips_downhill.resize(wire_to_tilewire.size());

    std::unordered_map<Arc, int> arc_to_pip;

    ArcVector arcs;
    ExtendedWireInfo ewi(ddb);
    for (auto i = 0u; i < wire_to_tilewire.size(); ++i) {
        const auto &tw = wire_to_tilewire[i];
        if (tw.isUndefined())
            continue;
        arcs.clear();

        const auto &tileInfo = tiles.getTileInfo(tw.getTileIndex());
        const auto tileTypeName = tiles.getTileTypeName(tileInfo.getTypeIndex());
        const bool clb = boost::starts_with(
                tileTypeName, "CLB"); // Disable all CLB route-throughs (i.e. LUT in->out, LUT A->AMUX, for now)

        const_cast<DDB &>(ddb).expandSegmentSinks(tw, arcs, DDB::eExpandDirectionNone, false /* inUseTied */,
                                                  true /*inUseRegular */, true /* inUseIrregular */,
                                                  !clb /* inUseRoutethrough */);

        auto index = pip_to_arc.size();
        pip_to_arc.insert(pip_to_arc.end(), arcs.begin(), arcs.end());

        const boost::regex bufg_i("(CMT|CLK)_BUFG_BUFGCTRL\\d+_I0");
        const boost::regex bufg_o("(CMT|CLK)_BUFG_BUFGCTRL\\d+_O");

        auto &pips = wire_to_pips_downhill[i];
        pips.reserve(arcs.size());
        const bool clk_tile = boost::starts_with(tileTypeName, "CMT") || boost::starts_with(tileTypeName, "CLK");
        for (const auto &a : arcs) {
            // Disable BUFG I0 -> O routethrough
            if (clk_tile) {
                ewi.set(a.getSourceTilewire());
                if (boost::regex_match(ewi.mWireName, bufg_i)) {
                    ewi.set(a.getSinkTilewire());
                    if (boost::regex_match(ewi.mWireName, bufg_o))
                        continue;
                }
            }
            pips.push_back(index);
            arc_to_pip.emplace(a, index);
            ++index;
        }
    }

    pip_to_arc.shrink_to_fit();

    wire_to_pips_uphill.resize(wire_to_tilewire.size());
    for (auto i = 0u; i < wire_to_tilewire.size(); ++i) {
        const auto &tw = wire_to_tilewire[i];
        if (tw.isUndefined())
            continue;
        arcs.clear();
        // TODO
        // const_cast<DDB&>(ddb).expandSegmentSources(tw, arcs, DDB::eExpandDirectionNone, false /* inUseTied */, true
        // /*inUseRegular */, true /* inUseIrregular */, false /* inUseRoutethrough */);

        auto &pips = wire_to_pips_uphill[i];
        pips.reserve(arcs.size());
        for (const auto &a : arcs)
            pips.push_back(arc_to_pip.at(a));
    }

    return pip_to_arc;
}

std::vector<int> construct_pip_to_dst_wire(const std::vector<Arc> &pip_to_arc)
{
    std::vector<int> pip_to_wire;
    return pip_to_wire;
}

// -----------------------------------------------------------------------

void IdString::initialize_arch(const BaseCtx *ctx)
{
#define X(t) initialize_add(ctx, #t, ID_##t);
#include "constids.inc"
#undef X
}

// -----------------------------------------------------------------------

Arch::Arch(ArchArgs args) : args(args)
{
    torc::common::DirectoryTree directoryTree("/opt/torc/src/torc");
    if (args.type == ArchArgs::Z020) {
        torc_info = std::unique_ptr<TorcInfo>(new TorcInfo(this, "xc7z020", args.package));
    } else {
        log_error("Unsupported XC7 chip type.\n");
    }

    //    package_info = nullptr;
    //    for (int i = 0; i < chip_info->num_packages; i++) {
    //        if (chip_info->packages_data[i].name.get() == args.package) {
    //            package_info = &(chip_info->packages_data[i]);
    //            break;
    //        }
    //    }
    //    if (package_info == nullptr)
    //        log_error("Unsupported package '%s'.\n", args.package.c_str());

    // bel_carry.resize(chip_info->num_bels);
    bel_to_cell.resize(torc_info->num_bels);
    wire_to_net.resize(torc_info->num_wires);
    pip_to_net.resize(torc_info->num_pips);
    // switches_locked.resize(chip_info->num_switches);
}

// -----------------------------------------------------------------------

std::string Arch::getChipName() const
{
    if (args.type == ArchArgs::Z020) {
        return "z020";
    } else {
        log_error("Unsupported XC7 chip type.\n");
    }
}

// -----------------------------------------------------------------------

IdString Arch::archArgsToId(ArchArgs args) const
{
    if (args.type == ArchArgs::Z020)
        return id("z020");
    return IdString();
}

// -----------------------------------------------------------------------

BelId Arch::getBelByName(IdString name) const
{
    auto it = torc_info->sites.findSiteIndex(name.str(this));
    if (it != SiteIndex(-1))
        return torc_info->site_index_to_bel.at(it);
    return BelId();
}

BelId Arch::getBelByLocation(Loc loc) const
{
    BelId bel;

    if (bel_by_loc.empty()) {
        for (int i = 0; i < torc_info->num_bels; i++) {
            BelId b;
            b.index = i;
            bel_by_loc[getBelLocation(b)] = b;
        }
    }

    auto it = bel_by_loc.find(loc);
    if (it != bel_by_loc.end())
        bel = it->second;

    return bel;
}

BelRange Arch::getBelsByTile(int x, int y) const
{
    BelRange br;

    br.b.cursor = Arch::getBelByLocation(Loc(x, y, 0)).index;
    br.e.cursor = br.b.cursor;

    if (br.e.cursor != -1) {
        while (br.e.cursor < chip_info->num_bels && chip_info->bel_data[br.e.cursor].x == x &&
               chip_info->bel_data[br.e.cursor].y == y)
            br.e.cursor++;
    }

    return br;
}

PortType Arch::getBelPinType(BelId bel, IdString pin) const
{
    NPNR_ASSERT(bel != BelId());

    int num_bel_wires = chip_info->bel_data[bel.index].num_bel_wires;
    const BelWirePOD *bel_wires = chip_info->bel_data[bel.index].bel_wires.get();

    if (num_bel_wires < 7) {
        for (int i = 0; i < num_bel_wires; i++) {
            if (bel_wires[i].port == pin.index)
                return PortType(bel_wires[i].type);
        }
    } else {
        int b = 0, e = num_bel_wires - 1;
        while (b <= e) {
            int i = (b + e) / 2;
            if (bel_wires[i].port == pin.index)
                return PortType(bel_wires[i].type);
            if (bel_wires[i].port > pin.index)
                e = i - 1;
            else
                b = i + 1;
        }
    }

    return PORT_INOUT;
}

WireId Arch::getBelPinWire(BelId bel, IdString pin) const
{
    WireId ret;

    auto pin_name = pin.str(this);
    auto bel_type = getBelType(bel);
    if (bel_type == id_SLICE_LUT6) {
        // For all LUT based inputs and outputs (I1-I6,O,OQ,OMUX) then change the I/O into the LUT
        if (pin_name[0] == 'I' || pin_name[0] == 'O') {
            switch (torc_info->bel_to_loc[bel.index].z) {
            case 0:
            case 4:
                pin_name[0] = 'A';
                break;
            case 1:
            case 5:
                pin_name[0] = 'B';
                break;
            case 2:
            case 6:
                pin_name[0] = 'C';
                break;
            case 3:
            case 7:
                pin_name[0] = 'D';
                break;
            default:
                throw;
            }
        }
    }
    else if (bel_type == id_PS7) {
        // e.g. Convert DDRARB[0] -> DDRARB0
        boost::erase_all(pin_name, "[");
        boost::erase_all(pin_name, "]");
    }

    auto site_index = torc_info->bel_to_site_index[bel.index];
    const auto &site = torc_info->sites.getSite(site_index);
    auto &tw = site.getPinTilewire(pin_name);

    if (tw.isUndefined())
        log_error("no wire found for site '%s' pin '%s' \n", torc_info->bel_to_name(bel.index).c_str(),
                  pin_name.c_str());

    ret.index = torc_info->tilewire_to_wire(tw);

    //    NPNR_ASSERT(bel != BelId());
    //
    //    int num_bel_wires = chip_info->bel_data[bel.index].num_bel_wires;
    //    const BelWirePOD *bel_wires = chip_info->bel_data[bel.index].bel_wires.get();
    //
    //    if (num_bel_wires < 7) {
    //        for (int i = 0; i < num_bel_wires; i++) {
    //            if (bel_wires[i].port == pin.index) {
    //                ret.index = bel_wires[i].wire_index;
    //                break;
    //            }
    //        }
    //    } else {
    //        int b = 0, e = num_bel_wires - 1;
    //        while (b <= e) {
    //            int i = (b + e) / 2;
    //            if (bel_wires[i].port == pin.index) {
    //                ret.index = bel_wires[i].wire_index;
    //                break;
    //            }
    //            if (bel_wires[i].port > pin.index)
    //                e = i - 1;
    //            else
    //                b = i + 1;
    //        }
    //    }

    return ret;
}

std::vector<IdString> Arch::getBelPins(BelId bel) const
{
    std::vector<IdString> ret;

    NPNR_ASSERT(bel != BelId());

    int num_bel_wires = chip_info->bel_data[bel.index].num_bel_wires;
    const BelWirePOD *bel_wires = chip_info->bel_data[bel.index].bel_wires.get();

    for (int i = 0; i < num_bel_wires; i++)
        ret.push_back(IdString(bel_wires[i].port));

    return ret;
}

// -----------------------------------------------------------------------

WireId Arch::getWireByName(IdString name) const
{
    WireId ret;

    if (wire_by_name.empty()) {
        for (int i = 0; i < chip_info->num_wires; i++)
            wire_by_name[id(chip_info->wire_data[i].name.get())] = i;
    }

    // auto it = wire_by_name.find(name);
    // if (it != wire_by_name.end())
    //    ret.index = it->second;

    return ret;
}

IdString Arch::getWireType(WireId wire) const
{
    NPNR_ASSERT(wire != WireId());
    //    switch (chip_info->wire_data[wire.index].type) {
    //    case WireInfoPOD::WIRE_TYPE_NONE:
    //        return IdString();
    //    case WireInfoPOD::WIRE_TYPE_GLB2LOCAL:
    //        return id("GLB2LOCAL");
    //    case WireInfoPOD::WIRE_TYPE_GLB_NETWK:
    //        return id("GLB_NETWK");
    //    case WireInfoPOD::WIRE_TYPE_LOCAL:
    //        return id("LOCAL");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_IN:
    //        return id("LUTFF_IN");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_IN_LUT:
    //        return id("LUTFF_IN_LUT");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_LOUT:
    //        return id("LUTFF_LOUT");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_OUT:
    //        return id("LUTFF_OUT");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_COUT:
    //        return id("LUTFF_COUT");
    //    case WireInfoPOD::WIRE_TYPE_LUTFF_GLOBAL:
    //        return id("LUTFF_GLOBAL");
    //    case WireInfoPOD::WIRE_TYPE_CARRY_IN_MUX:
    //        return id("CARRY_IN_MUX");
    //    case WireInfoPOD::WIRE_TYPE_SP4_V:
    //        return id("SP4_V");
    //    case WireInfoPOD::WIRE_TYPE_SP4_H:
    //        return id("SP4_H");
    //    case WireInfoPOD::WIRE_TYPE_SP12_V:
    //        return id("SP12_V");
    //    case WireInfoPOD::WIRE_TYPE_SP12_H:
    //        return id("SP12_H");
    //    }
    return IdString();
}

// -----------------------------------------------------------------------

PipId Arch::getPipByName(IdString name) const
{
    PipId ret;

    if (pip_by_name.empty()) {
        for (int i = 0; i < chip_info->num_pips; i++) {
            PipId pip;
            pip.index = i;
            pip_by_name[getPipName(pip)] = i;
        }
    }

    auto it = pip_by_name.find(name);
    if (it != pip_by_name.end())
        ret.index = it->second;

    return ret;
}

IdString Arch::getPipName(PipId pip) const
{
    NPNR_ASSERT(pip != PipId());

    ExtendedWireInfo ewi_src(*torc_info->ddb, torc_info->pip_to_arc[pip.index].getSourceTilewire());
    ExtendedWireInfo ewi_dst(*torc_info->ddb, torc_info->pip_to_arc[pip.index].getSinkTilewire());
    std::stringstream pip_name;
    pip_name << ewi_src.mTileName << "." << ewi_src.mWireName << ".->." << ewi_dst.mWireName;
    return id(pip_name.str());

//#if 1
//    int x = chip_info->pip_data[pip.index].x;
//    int y = chip_info->pip_data[pip.index].y;
//
//    std::string src_name = chip_info->wire_data[chip_info->pip_data[pip.index].src].name.get();
//    std::replace(src_name.begin(), src_name.end(), '/', '.');
//
//    std::string dst_name = chip_info->wire_data[chip_info->pip_data[pip.index].dst].name.get();
//    std::replace(dst_name.begin(), dst_name.end(), '/', '.');
//
//    return id("X" + std::to_string(x) + "/Y" + std::to_string(y) + "/" + src_name + ".->." + dst_name);
//#else
//    return id(chip_info->pip_data[pip.index].name.get());
//#endif
}

// -----------------------------------------------------------------------

BelId Arch::getPackagePinBel(const std::string &pin) const
{
    return getBelByName(id(pin));
}

std::string Arch::getBelPackagePin(BelId bel) const
{
    //    for (int i = 0; i < package_info->num_pins; i++) {
    //        if (package_info->pins[i].bel_index == bel.index) {
    //            return std::string(package_info->pins[i].name.get());
    //        }
    //    }
    return "";
}

// -----------------------------------------------------------------------

GroupId Arch::getGroupByName(IdString name) const
{
    for (auto g : getGroups())
        if (getGroupName(g) == name)
            return g;
    return GroupId();
}

IdString Arch::getGroupName(GroupId group) const
{
    std::string suffix;

    switch (group.type) {
    case GroupId::TYPE_FRAME:
        suffix = "tile";
        break;
    case GroupId::TYPE_MAIN_SW:
        suffix = "main_sw";
        break;
    case GroupId::TYPE_LOCAL_SW:
        suffix = "local_sw";
        break;
    case GroupId::TYPE_LC0_SW:
        suffix = "lc0_sw";
        break;
    case GroupId::TYPE_LC1_SW:
        suffix = "lc1_sw";
        break;
    case GroupId::TYPE_LC2_SW:
        suffix = "lc2_sw";
        break;
    case GroupId::TYPE_LC3_SW:
        suffix = "lc3_sw";
        break;
    case GroupId::TYPE_LC4_SW:
        suffix = "lc4_sw";
        break;
    case GroupId::TYPE_LC5_SW:
        suffix = "lc5_sw";
        break;
    case GroupId::TYPE_LC6_SW:
        suffix = "lc6_sw";
        break;
    case GroupId::TYPE_LC7_SW:
        suffix = "lc7_sw";
        break;
    default:
        return IdString();
    }

    return id("X" + std::to_string(group.x) + "/Y" + std::to_string(group.y) + "/" + suffix);
}

std::vector<GroupId> Arch::getGroups() const
{
    std::vector<GroupId> ret;

    for (int y = 0; y < chip_info->height; y++) {
        for (int x = 0; x < chip_info->width; x++) {
            TileType type = chip_info->tile_grid[y * chip_info->width + x];
            if (type == TILE_NONE)
                continue;

            GroupId group;
            group.type = GroupId::TYPE_FRAME;
            group.x = x;
            group.y = y;
            // ret.push_back(group);

            group.type = GroupId::TYPE_MAIN_SW;
            ret.push_back(group);

            group.type = GroupId::TYPE_LOCAL_SW;
            ret.push_back(group);

            if (type == TILE_LOGIC) {
                group.type = GroupId::TYPE_LC0_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC1_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC2_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC3_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC4_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC5_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC6_SW;
                ret.push_back(group);

                group.type = GroupId::TYPE_LC7_SW;
                ret.push_back(group);
            }
        }
    }
    return ret;
}

std::vector<BelId> Arch::getGroupBels(GroupId group) const
{
    std::vector<BelId> ret;
    return ret;
}

std::vector<WireId> Arch::getGroupWires(GroupId group) const
{
    std::vector<WireId> ret;
    return ret;
}

std::vector<PipId> Arch::getGroupPips(GroupId group) const
{
    std::vector<PipId> ret;
    return ret;
}

std::vector<GroupId> Arch::getGroupGroups(GroupId group) const
{
    std::vector<GroupId> ret;
    return ret;
}

// -----------------------------------------------------------------------

bool Arch::getBudgetOverride(const NetInfo *net_info, const PortRef &sink, delay_t &budget) const { return false; }

// -----------------------------------------------------------------------

bool Arch::place() { return placer1(getCtx(), Placer1Cfg(getCtx())); }

bool Arch::route()
{
    return router1(getCtx(), Router1Cfg(getCtx()));
}

// -----------------------------------------------------------------------

DecalXY Arch::getBelDecal(BelId bel) const
{
    DecalXY decalxy;
    decalxy.decal.type = DecalId::TYPE_BEL;
    decalxy.decal.index = bel.index;
    decalxy.decal.active = bel_to_cell.at(bel.index) != nullptr;
    return decalxy;
}

DecalXY Arch::getWireDecal(WireId wire) const
{
    DecalXY decalxy;
    decalxy.decal.type = DecalId::TYPE_WIRE;
    decalxy.decal.index = wire.index;
    decalxy.decal.active = wire_to_net.at(wire.index) != nullptr;
    return decalxy;
}

DecalXY Arch::getPipDecal(PipId pip) const
{
    DecalXY decalxy;
    decalxy.decal.type = DecalId::TYPE_PIP;
    decalxy.decal.index = pip.index;
    decalxy.decal.active = pip_to_net.at(pip.index) != nullptr;
    return decalxy;
};

DecalXY Arch::getGroupDecal(GroupId group) const
{
    DecalXY decalxy;
    decalxy.decal.type = DecalId::TYPE_GROUP;
    decalxy.decal.index = (group.type << 16) | (group.x << 8) | (group.y);
    decalxy.decal.active = true;
    return decalxy;
};

std::vector<GraphicElement> Arch::getDecalGraphics(DecalId decal) const
{
    std::vector<GraphicElement> ret;

    if (decal.type == DecalId::TYPE_GROUP) {
        int type = (decal.index >> 16) & 255;
        int x = (decal.index >> 8) & 255;
        int y = decal.index & 255;

        if (type == GroupId::TYPE_FRAME) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_LINE;
            el.style = GraphicElement::STYLE_FRAME;

            el.x1 = x + 0.01, el.x2 = x + 0.02, el.y1 = y + 0.01, el.y2 = y + 0.01;
            ret.push_back(el);
            el.x1 = x + 0.01, el.x2 = x + 0.01, el.y1 = y + 0.01, el.y2 = y + 0.02;
            ret.push_back(el);

            el.x1 = x + 0.99, el.x2 = x + 0.98, el.y1 = y + 0.01, el.y2 = y + 0.01;
            ret.push_back(el);
            el.x1 = x + 0.99, el.x2 = x + 0.99, el.y1 = y + 0.01, el.y2 = y + 0.02;
            ret.push_back(el);

            el.x1 = x + 0.99, el.x2 = x + 0.98, el.y1 = y + 0.99, el.y2 = y + 0.99;
            ret.push_back(el);
            el.x1 = x + 0.99, el.x2 = x + 0.99, el.y1 = y + 0.99, el.y2 = y + 0.98;
            ret.push_back(el);

            el.x1 = x + 0.01, el.x2 = x + 0.02, el.y1 = y + 0.99, el.y2 = y + 0.99;
            ret.push_back(el);
            el.x1 = x + 0.01, el.x2 = x + 0.01, el.y1 = y + 0.99, el.y2 = y + 0.98;
            ret.push_back(el);
        }

        if (type == GroupId::TYPE_MAIN_SW) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_BOX;
            el.style = GraphicElement::STYLE_FRAME;

            el.x1 = x + main_swbox_x1;
            el.x2 = x + main_swbox_x2;
            el.y1 = y + main_swbox_y1;
            el.y2 = y + main_swbox_y2;
            ret.push_back(el);
        }

        if (type == GroupId::TYPE_LOCAL_SW) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_BOX;
            el.style = GraphicElement::STYLE_FRAME;

            el.x1 = x + local_swbox_x1;
            el.x2 = x + local_swbox_x2;
            el.y1 = y + local_swbox_y1;
            el.y2 = y + local_swbox_y2;
            ret.push_back(el);
        }

        if (GroupId::TYPE_LC0_SW <= type && type <= GroupId::TYPE_LC7_SW) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_BOX;
            el.style = GraphicElement::STYLE_FRAME;

            el.x1 = x + lut_swbox_x1;
            el.x2 = x + lut_swbox_x2;
            el.y1 = y + logic_cell_y1 + logic_cell_pitch * (type - GroupId::TYPE_LC0_SW);
            el.y2 = y + logic_cell_y2 + logic_cell_pitch * (type - GroupId::TYPE_LC0_SW);
            ret.push_back(el);
        }
    }

    if (decal.type == DecalId::TYPE_WIRE) {
        int n = chip_info->wire_data[decal.index].num_segments;
        const WireSegmentPOD *p = chip_info->wire_data[decal.index].segments.get();

        GraphicElement::style_t style = decal.active ? GraphicElement::STYLE_ACTIVE : GraphicElement::STYLE_INACTIVE;

        for (int i = 0; i < n; i++)
            gfxTileWire(ret, p[i].x, p[i].y, GfxTileWireId(p[i].index), style);
    }

    if (decal.type == DecalId::TYPE_PIP) {
        const PipInfoPOD &p = chip_info->pip_data[decal.index];
        GraphicElement::style_t style = decal.active ? GraphicElement::STYLE_ACTIVE : GraphicElement::STYLE_HIDDEN;
        gfxTilePip(ret, p.x, p.y, GfxTileWireId(p.src_seg), GfxTileWireId(p.dst_seg), style);
    }

    if (decal.type == DecalId::TYPE_BEL) {
        BelId bel;
        bel.index = SiteIndex(decal.index);

        auto bel_type = getBelType(bel);

        if (bel_type == id_ICESTORM_LC) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_BOX;
            el.style = decal.active ? GraphicElement::STYLE_ACTIVE : GraphicElement::STYLE_INACTIVE;
            el.x1 = chip_info->bel_data[bel.index].x + logic_cell_x1;
            el.x2 = chip_info->bel_data[bel.index].x + logic_cell_x2;
            el.y1 = chip_info->bel_data[bel.index].y + logic_cell_y1 +
                    (chip_info->bel_data[bel.index].z) * logic_cell_pitch;
            el.y2 = chip_info->bel_data[bel.index].y + logic_cell_y2 +
                    (chip_info->bel_data[bel.index].z) * logic_cell_pitch;
            ret.push_back(el);
        }

        if (bel_type == id_SB_IO) {
            GraphicElement el;
            el.type = GraphicElement::TYPE_BOX;
            el.style = decal.active ? GraphicElement::STYLE_ACTIVE : GraphicElement::STYLE_INACTIVE;
            el.x1 = chip_info->bel_data[bel.index].x + logic_cell_x1;
            el.x2 = chip_info->bel_data[bel.index].x + logic_cell_x2;
            el.y1 = chip_info->bel_data[bel.index].y + logic_cell_y1 +
                    (4 * chip_info->bel_data[bel.index].z) * logic_cell_pitch;
            el.y2 = chip_info->bel_data[bel.index].y + logic_cell_y2 +
                    (4 * chip_info->bel_data[bel.index].z + 3) * logic_cell_pitch;
            ret.push_back(el);
        }

        if (bel_type == id_ICESTORM_RAM) {
            for (int i = 0; i < 2; i++) {
                GraphicElement el;
                el.type = GraphicElement::TYPE_BOX;
                el.style = decal.active ? GraphicElement::STYLE_ACTIVE : GraphicElement::STYLE_INACTIVE;
                el.x1 = chip_info->bel_data[bel.index].x + logic_cell_x1;
                el.x2 = chip_info->bel_data[bel.index].x + logic_cell_x2;
                el.y1 = chip_info->bel_data[bel.index].y + logic_cell_y1 + i;
                el.y2 = chip_info->bel_data[bel.index].y + logic_cell_y2 + i + 7 * logic_cell_pitch;
                ret.push_back(el);
            }
        }
    }

    return ret;
}

// -----------------------------------------------------------------------

bool Arch::getCellDelay(const CellInfo *cell, IdString fromPort, IdString toPort, DelayInfo &delay) const
{
    if (cell->type == id_SLICE_LUT6) {
        if (fromPort.index >= id_I1.index && fromPort.index <= id_I6.index) {
            if (toPort == id_O) {
                delay.delay = 124; // Tilo
                return true;
            }
            if (toPort == id_OQ) {
                delay.delay = 95; // Tas
                return true;
            }
        }
        if (fromPort == id_CLK) {
            if (toPort == id_OQ) {
                delay.delay = 456; // Tcko
                return true;
            }
        }
    } else if (cell->type == id_BUFGCTRL) {
        return true;
    }
    return false;
}

// Get the port class, also setting clockPort to associated clock if applicable
TimingPortClass Arch::getPortTimingClass(const CellInfo *cell, IdString port, int &clockInfoCount) const
{
    if (cell->type == id_SLICE_LUT6) {
        if (port == id_CLK)
            return TMG_CLOCK_INPUT;
        if (port == id_CIN)
            return TMG_COMB_INPUT;
        if (port == id_COUT)
            return TMG_COMB_OUTPUT;
        if (port == id_O) {
            // LCs with no inputs are constant drivers
            if (cell->lcInfo.inputCount == 0)
                return TMG_IGNORE;
            return TMG_COMB_OUTPUT;
        }
        if (cell->lcInfo.dffEnable) {
            clockInfoCount = 1;
            if (port == id_OQ)
                return TMG_REGISTER_OUTPUT;
            return TMG_REGISTER_INPUT;
        } else {
            return TMG_COMB_INPUT;
        }
        // TODO
        // if (port == id_OMUX)
    } else if (cell->type == id_IOB33) {
        if (port == id_I)
            return TMG_STARTPOINT;
        else if (port == id_O)
            return TMG_ENDPOINT;
    } else if (cell->type == id_BUFGCTRL) {
        if (port == id_O)
            return TMG_COMB_OUTPUT;
        return TMG_COMB_INPUT;
    }
    else if (cell->type == id_PS7) {
        // TODO
        return TMG_IGNORE;
    }
    log_error("no timing info for port '%s' of cell type '%s'\n", port.c_str(this), cell->type.c_str(this));
}

TimingClockingInfo Arch::getPortClockingInfo(const CellInfo *cell, IdString port, int index) const
{
    TimingClockingInfo info;
    if (cell->type == id_SLICE_LUT6) {
        info.clock_port = id_CLK;
        info.edge = cell->lcInfo.negClk ? FALLING_EDGE : RISING_EDGE;
        if (port == id_OQ) {
            bool has_clktoq = getCellDelay(cell, id_CLK, id_OQ, info.clockToQ);
            NPNR_ASSERT(has_clktoq);
        } else {
            info.setup.delay = 124; // Tilo
            info.hold.delay = 0;
        }
    } else {
        NPNR_ASSERT_FALSE("unhandled cell type in getPortClockingInfo");
    }
    return info;
}

bool Arch::isGlobalNet(const NetInfo *net) const
{
    if (net == nullptr)
        return false;
    return net->driver.cell != nullptr && net->driver.cell->type == id_BUFGCTRL && net->driver.port == id_O;
}

// Assign arch arg info
void Arch::assignArchInfo()
{
    for (auto &net : getCtx()->nets) {
        NetInfo *ni = net.second.get();
        if (isGlobalNet(ni))
            ni->is_global = true;
        ni->is_enable = false;
        ni->is_reset = false;
        for (auto usr : ni->users) {
            if (is_enable_port(this, usr))
                ni->is_enable = true;
            if (is_reset_port(this, usr))
                ni->is_reset = true;
        }
    }
    for (auto &cell : getCtx()->cells) {
        CellInfo *ci = cell.second.get();
        assignCellInfo(ci);
    }
}

void Arch::assignCellInfo(CellInfo *cell)
{
    cell->belType = cell->type;
    if (cell->type == id_SLICE_LUT6) {
        cell->lcInfo.dffEnable = bool_or_default(cell->params, id_DFF_ENABLE);
        cell->lcInfo.carryEnable = bool_or_default(cell->params, id_CARRY_ENABLE);
        cell->lcInfo.negClk = bool_or_default(cell->params, id_NEG_CLK);
        cell->lcInfo.clk = get_net_or_empty(cell, id_CLK);
        cell->lcInfo.cen = get_net_or_empty(cell, id_CEN);
        cell->lcInfo.sr = get_net_or_empty(cell, id_SR);
        cell->lcInfo.inputCount = 0;
        if (get_net_or_empty(cell, id_I1))
            cell->lcInfo.inputCount++;
        if (get_net_or_empty(cell, id_I2))
            cell->lcInfo.inputCount++;
        if (get_net_or_empty(cell, id_I3))
            cell->lcInfo.inputCount++;
        if (get_net_or_empty(cell, id_I4))
            cell->lcInfo.inputCount++;
        if (get_net_or_empty(cell, id_I5))
            cell->lcInfo.inputCount++;
        if (get_net_or_empty(cell, id_I6))
            cell->lcInfo.inputCount++;
    }
}

NEXTPNR_NAMESPACE_END