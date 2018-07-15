/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  Clifford Wolf <clifford@symbioticeda.com>
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
#include "log.h"
#include "nextpnr.h"
#include "util.h"
#include "cells.h"
NEXTPNR_NAMESPACE_BEGIN

// -----------------------------------------------------------------------

IdString Arch::belTypeToId(BelType type) const
{
    if (type == TYPE_ICESTORM_LC)
        return id("ICESTORM_LC");
    if (type == TYPE_ICESTORM_RAM)
        return id("ICESTORM_RAM");
    if (type == TYPE_SB_IO)
        return id("SB_IO");
    if (type == TYPE_SB_GB)
        return id("SB_GB");
    if (type == TYPE_ICESTORM_PLL)
        return id("ICESTORM_PLL");
    if (type == TYPE_SB_WARMBOOT)
        return id("SB_WARMBOOT");
    if (type == TYPE_SB_MAC16)
        return id("SB_MAC16");
    if (type == TYPE_ICESTORM_HFOSC)
        return id("ICESTORM_HFOSC");
    if (type == TYPE_ICESTORM_LFOSC)
        return id("ICESTORM_LFOSC");
    if (type == TYPE_SB_I2C)
        return id("SB_I2C");
    if (type == TYPE_SB_SPI)
        return id("SB_SPI");
    if (type == TYPE_IO_I3C)
        return id("IO_I3C");
    if (type == TYPE_SB_LEDDA_IP)
        return id("SB_LEDDA_IP");
    if (type == TYPE_SB_RGBA_DRV)
        return id("SB_RGBA_DRV");
    if (type == TYPE_ICESTORM_SPRAM)
        return id("ICESTORM_SPRAM");
    return IdString();
}

BelType Arch::belTypeFromId(IdString type) const
{
    if (type == id("ICESTORM_LC"))
        return TYPE_ICESTORM_LC;
    if (type == id("ICESTORM_RAM"))
        return TYPE_ICESTORM_RAM;
    if (type == id("SB_IO"))
        return TYPE_SB_IO;
    if (type == id("SB_GB"))
        return TYPE_SB_GB;
    if (type == id("ICESTORM_PLL"))
        return TYPE_ICESTORM_PLL;
    if (type == id("SB_WARMBOOT"))
        return TYPE_SB_WARMBOOT;
    if (type == id("SB_MAC16"))
        return TYPE_SB_MAC16;
    if (type == id("ICESTORM_HFOSC"))
        return TYPE_ICESTORM_HFOSC;
    if (type == id("ICESTORM_LFOSC"))
        return TYPE_ICESTORM_LFOSC;
    if (type == id("SB_I2C"))
        return TYPE_SB_I2C;
    if (type == id("SB_SPI"))
        return TYPE_SB_SPI;
    if (type == id("IO_I3C"))
        return TYPE_IO_I3C;
    if (type == id("SB_LEDDA_IP"))
        return TYPE_SB_LEDDA_IP;
    if (type == id("SB_RGBA_DRV"))
        return TYPE_SB_RGBA_DRV;
    if (type == id("ICESTORM_SPRAM"))
        return TYPE_ICESTORM_SPRAM;
    return TYPE_NONE;
}

// -----------------------------------------------------------------------

void IdString::initialize_arch(const BaseCtx *ctx)
{
#define X(t) initialize_add(ctx, #t, PIN_##t);
#include "portpins.inc"
#undef X
}

IdString Arch::portPinToId(PortPin type) const
{
    IdString ret;
    if (type > 0 && type < PIN_MAXIDX)
        ret.index = type;
    return ret;
}

PortPin Arch::portPinFromId(IdString type) const
{
    if (type.index > 0 && type.index < PIN_MAXIDX)
        return PortPin(type.index);
    return PIN_NONE;
}

// -----------------------------------------------------------------------

static const ChipInfoPOD *get_chip_info(const RelPtr<ChipInfoPOD> *ptr) { return ptr->get(); }

#if defined(_MSC_VER)
void load_chipdb();
#endif

Arch::Arch(ArchArgs args) : args(args)
{
#if defined(_MSC_VER)
    load_chipdb();
#endif

#ifdef ICE40_HX1K_ONLY
    if (args.type == ArchArgs::HX1K) {
        chip_info = get_chip_info(reinterpret_cast<const RelPtr<ChipInfoPOD> *>(chipdb_blob_1k));
    } else {
        log_error("Unsupported iCE40 chip type.\n");
    }
#else
    if (args.type == ArchArgs::LP384) {
        chip_info = get_chip_info(reinterpret_cast<const RelPtr<ChipInfoPOD> *>(chipdb_blob_384));
    } else if (args.type == ArchArgs::LP1K || args.type == ArchArgs::HX1K) {
        chip_info = get_chip_info(reinterpret_cast<const RelPtr<ChipInfoPOD> *>(chipdb_blob_1k));
    } else if (args.type == ArchArgs::UP5K) {
        chip_info = get_chip_info(reinterpret_cast<const RelPtr<ChipInfoPOD> *>(chipdb_blob_5k));
    } else if (args.type == ArchArgs::LP8K || args.type == ArchArgs::HX8K) {
        chip_info = get_chip_info(reinterpret_cast<const RelPtr<ChipInfoPOD> *>(chipdb_blob_8k));
    } else {
        log_error("Unsupported iCE40 chip type.\n");
    }
#endif

    package_info = nullptr;
    for (int i = 0; i < chip_info->num_packages; i++) {
        if (chip_info->packages_data[i].name.get() == args.package) {
            package_info = &(chip_info->packages_data[i]);
            break;
        }
    }
    if (package_info == nullptr)
        log_error("Unsupported package '%s'.\n", args.package.c_str());

    bel_to_cell.resize(chip_info->num_bels);
    wire_to_net.resize(chip_info->num_wires);
    pip_to_net.resize(chip_info->num_pips);
    switches_locked.resize(chip_info->num_switches);

    // Initialise regularly used IDStrings for performance
    id_glb_buf_out = id("GLOBAL_BUFFER_OUTPUT");
    id_icestorm_lc = id("ICESTORM_LC");
    id_sb_io = id("SB_IO");
    id_sb_gb = id("SB_GB");
    id_cen = id("CEN");
    id_clk = id("CLK");
    id_sr = id("SR");
    id_i0 = id("I0");
    id_i1 = id("I1");
    id_i2 = id("I2");
    id_i3 = id("I3");
    id_dff_en = id("DFF_ENABLE");
    id_neg_clk = id("NEG_CLK");
    id_r = id("R");
    id_s = id("S");
    id_e = id("E");
    id_set_ff = {
            id("SB_DFF"), id("SB_DFFE"), id("SB_DFFSR"),
            id("SB_DFFR"), id("SB_DFFSS"), id("SB_DFFS"),
            id("SB_DFFESR"), id("SB_DFFER"),
            id("SB_DFFESS"), id("SB_DFFES"),
            id("SB_DFFN"), id("SB_DFFNE"),
            id("SB_DFFNSR"), id("SB_DFFNR"),
            id("SB_DFFNSS"), id("SB_DFFNS"),
            id("SB_DFFNESR"), id("SB_DFFNER"),
            id("SB_DFFNESS"), id("SB_DFFNES") };
}

// -----------------------------------------------------------------------

std::string Arch::getChipName()
{
#ifdef ICE40_HX1K_ONLY
    if (args.type == ArchArgs::HX1K) {
        return "Lattice LP1K";
    } else {
        log_error("Unsupported iCE40 chip type.\n");
    }
#else
    if (args.type == ArchArgs::LP384) {
        return "Lattice LP384";
    } else if (args.type == ArchArgs::LP1K) {
        return "Lattice LP1K";
    } else if (args.type == ArchArgs::HX1K) {
        return "Lattice HX1K";
    } else if (args.type == ArchArgs::UP5K) {
        return "Lattice UP5K";
    } else if (args.type == ArchArgs::LP8K) {
        return "Lattice LP8K";
    } else if (args.type == ArchArgs::HX8K) {
        return "Lattice HX8K";
    } else {
        log_error("Unknown chip\n");
    }
#endif
}

// -----------------------------------------------------------------------

IdString Arch::archArgsToId(ArchArgs args) const
{
    if (args.type == ArchArgs::LP384)
        return id("lp384");
    if (args.type == ArchArgs::LP1K)
        return id("lp1k");
    if (args.type == ArchArgs::HX1K)
        return id("hx1k");
    if (args.type == ArchArgs::UP5K)
        return id("up5k");
    if (args.type == ArchArgs::LP8K)
        return id("lp8k");
    if (args.type == ArchArgs::HX8K)
        return id("hx8k");
    return IdString();
}

// -----------------------------------------------------------------------

BelId Arch::getBelByName(IdString name) const
{
    BelId ret;

    if (bel_by_name.empty()) {
        for (int i = 0; i < chip_info->num_bels; i++)
            bel_by_name[id(chip_info->bel_data[i].name.get())] = i;
    }

    auto it = bel_by_name.find(name);
    if (it != bel_by_name.end())
        ret.index = it->second;

    return ret;
}

BelRange Arch::getBelsAtSameTile(BelId bel) const
{
    BelRange br;
    NPNR_ASSERT(bel != BelId());
    // This requires Bels at the same tile are consecutive
    int x = chip_info->bel_data[bel.index].x;
    int y = chip_info->bel_data[bel.index].y;
    int start = bel.index, end = bel.index;
    while (start >= 0 && chip_info->bel_data[start].x == x && chip_info->bel_data[start].y == y)
        start--;
    start++;
    br.b.cursor = start;
    while (end < chip_info->num_bels && chip_info->bel_data[end].x == x && chip_info->bel_data[end].y == y)
        end++;
    br.e.cursor = end;
    return br;
}

WireId Arch::getWireBelPin(BelId bel, PortPin pin) const
{
    WireId ret;

    NPNR_ASSERT(bel != BelId());

    int num_bel_wires = chip_info->bel_data[bel.index].num_bel_wires;
    const BelWirePOD *bel_wires = chip_info->bel_data[bel.index].bel_wires.get();

    for (int i = 0; i < num_bel_wires; i++)
        if (bel_wires[i].port == pin) {
            ret.index = bel_wires[i].wire_index;
            break;
        }

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

    auto it = wire_by_name.find(name);
    if (it != wire_by_name.end())
        ret.index = it->second;

    return ret;
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

    int x = chip_info->pip_data[pip.index].x;
    int y = chip_info->pip_data[pip.index].y;

    std::string src_name = chip_info->wire_data[chip_info->pip_data[pip.index].src].name.get();
    std::replace(src_name.begin(), src_name.end(), '/', '.');

    std::string dst_name = chip_info->wire_data[chip_info->pip_data[pip.index].dst].name.get();
    std::replace(dst_name.begin(), dst_name.end(), '/', '.');

    return id("X" + std::to_string(x) + "/Y" + std::to_string(y) + "/" + src_name + ".->." + dst_name);
}

// -----------------------------------------------------------------------

BelId Arch::getPackagePinBel(const std::string &pin) const
{
    for (int i = 0; i < package_info->num_pins; i++) {
        if (package_info->pins[i].name.get() == pin) {
            BelId id;
            id.index = package_info->pins[i].bel_index;
            return id;
        }
    }
    return BelId();
}

std::string Arch::getBelPackagePin(BelId bel) const
{
    for (int i = 0; i < package_info->num_pins; i++) {
        if (package_info->pins[i].bel_index == bel.index) {
            return std::string(package_info->pins[i].name.get());
        }
    }
    return "";
}
// -----------------------------------------------------------------------

void Arch::estimatePosition(BelId bel, int &x, int &y, bool &gb) const
{
    NPNR_ASSERT(bel != BelId());
    x = chip_info->bel_data[bel.index].x;
    y = chip_info->bel_data[bel.index].y;
    gb = chip_info->bel_data[bel.index].type == TYPE_SB_GB;
}

delay_t Arch::estimateDelay(WireId src, WireId dst) const
{
    NPNR_ASSERT(src != WireId());
    int x1 = chip_info->wire_data[src.index].x;
    int y1 = chip_info->wire_data[src.index].y;

    NPNR_ASSERT(dst != WireId());
    int x2 = chip_info->wire_data[dst.index].x;
    int y2 = chip_info->wire_data[dst.index].y;

    int xd = x2 - x1, yd = y2 - y1;
    int xscale = 120, yscale = 120, offset = 0;

    // if (chip_info->wire_data[src.index].type == WIRE_TYPE_SP4_VERT) {
    //     yd = yd < -4 ? yd + 4 : (yd < 0 ? 0 : yd);
    //     offset = 500;
    // }

    return xscale * abs(xd) + yscale * abs(yd) + offset;
}

// -----------------------------------------------------------------------

std::vector<GraphicElement> Arch::getFrameGraphics() const
{
    std::vector<GraphicElement> ret;

    for (int x = 0; x <= chip_info->width; x++)
        for (int y = 0; y <= chip_info->height; y++) {
            GraphicElement el;
            el.type = GraphicElement::G_LINE;
            el.x1 = x - 0.05, el.x2 = x + 0.05, el.y1 = y, el.y2 = y, el.z = 0;
            ret.push_back(el);
            el.x1 = x, el.x2 = x, el.y1 = y - 0.05, el.y2 = y + 0.05, el.z = 0;
            ret.push_back(el);
        }

    return ret;
}

std::vector<GraphicElement> Arch::getBelGraphics(BelId bel) const
{
    std::vector<GraphicElement> ret;

    auto bel_type = getBelType(bel);

    if (bel_type == TYPE_ICESTORM_LC) {
        GraphicElement el;
        el.type = GraphicElement::G_BOX;
        el.x1 = chip_info->bel_data[bel.index].x + 0.1;
        el.x2 = chip_info->bel_data[bel.index].x + 0.9;
        el.y1 = chip_info->bel_data[bel.index].y + 0.10 + (chip_info->bel_data[bel.index].z) * (0.8 / 8);
        el.y2 = chip_info->bel_data[bel.index].y + 0.18 + (chip_info->bel_data[bel.index].z) * (0.8 / 8);
        el.z = 0;
        ret.push_back(el);
    }

    if (bel_type == TYPE_SB_IO) {
        if (chip_info->bel_data[bel.index].x == 0 || chip_info->bel_data[bel.index].x == chip_info->width - 1) {
            GraphicElement el;
            el.type = GraphicElement::G_BOX;
            el.x1 = chip_info->bel_data[bel.index].x + 0.1;
            el.x2 = chip_info->bel_data[bel.index].x + 0.9;
            if (chip_info->bel_data[bel.index].z == 0) {
                el.y1 = chip_info->bel_data[bel.index].y + 0.10;
                el.y2 = chip_info->bel_data[bel.index].y + 0.45;
            } else {
                el.y1 = chip_info->bel_data[bel.index].y + 0.55;
                el.y2 = chip_info->bel_data[bel.index].y + 0.90;
            }
            el.z = 0;
            ret.push_back(el);
        } else {
            GraphicElement el;
            el.type = GraphicElement::G_BOX;
            if (chip_info->bel_data[bel.index].z == 0) {
                el.x1 = chip_info->bel_data[bel.index].x + 0.10;
                el.x2 = chip_info->bel_data[bel.index].x + 0.45;
            } else {
                el.x1 = chip_info->bel_data[bel.index].x + 0.55;
                el.x2 = chip_info->bel_data[bel.index].x + 0.90;
            }
            el.y1 = chip_info->bel_data[bel.index].y + 0.1;
            el.y2 = chip_info->bel_data[bel.index].y + 0.9;
            el.z = 0;
            ret.push_back(el);
        }
    }

    if (bel_type == TYPE_ICESTORM_RAM) {
        GraphicElement el;
        el.type = GraphicElement::G_BOX;
        el.x1 = chip_info->bel_data[bel.index].x + 0.1;
        el.x2 = chip_info->bel_data[bel.index].x + 0.9;
        el.y1 = chip_info->bel_data[bel.index].y + 0.1;
        el.y2 = chip_info->bel_data[bel.index].y + 1.9;
        el.z = 0;
        ret.push_back(el);
    }

    return ret;
}

std::vector<GraphicElement> Arch::getWireGraphics(WireId wire) const
{
    std::vector<GraphicElement> ret;
    // FIXME
    return ret;
}

std::vector<GraphicElement> Arch::getPipGraphics(PipId pip) const
{
    std::vector<GraphicElement> ret;
    // FIXME
    return ret;
};

// -----------------------------------------------------------------------

bool Arch::getCellDelay(const CellInfo *cell, IdString fromPort, IdString toPort, delay_t &delay) const
{
    if (cell->type == id("ICESTORM_LC")) {
        if ((fromPort == id("I0") || fromPort == id("I1") || fromPort == id("I2") || fromPort == id("I3")) &&
            (toPort == id("O") || toPort == id("LO"))) {
            delay = 450;
            return true;
        } else if (fromPort == id("CIN") && toPort == id("COUT")) {
            delay = 120;
            return true;
        } else if (fromPort == id("I1") && toPort == id("COUT")) {
            delay = 260;
            return true;
        } else if (fromPort == id("I2") && toPort == id("COUT")) {
            delay = 230;
            return true;
        }
    }
    return false;
}

IdString Arch::getPortClock(const CellInfo *cell, IdString port) const
{
    if (cell->type == id("ICESTORM_LC") && bool_or_default(cell->params, id("DFF_ENABLE"))) {
        if (port != id("LO") && port != id("CIN") && port != id("COUT"))
            return id("CLK");
    }
    return IdString();
}

bool Arch::isClockPort(const CellInfo *cell, IdString port) const
{
    if (cell->type == id("ICESTORM_LC") && port == id("CLK"))
        return true;
    return false;
}

bool Arch::isGlobalNet(const NetInfo *net) const
{
    if (net == nullptr)
        return false;
    return net->driver.cell != nullptr && net->driver.port == id_glb_buf_out;
}

bool Arch::isIO(const CellInfo* cell) const
{
    return cell->type == id("SB_IO");
}

bool Arch::isClockPort(const PortRef &port) const
{
    if (port.cell == nullptr)
        return false;
    if (isFF(port.cell))
        return port.port == id("C");
    if (port.cell->type == id("ICESTORM_LC"))
        return port.port == id("CLK");
    if (is_ram(this, port.cell) || port.cell->type == id("ICESTORM_RAM"))
        return port.port == id("RCLK") || port.port == id("WCLK");
    return false;
}

bool Arch::isResetPort(const PortRef &port) const
{
    if (port.cell == nullptr)
        return false;
    if (isFF(port.cell))
        return port.port == id_r || port.port == id_s;
    if (port.cell->type == id_icestorm_lc)
        return port.port == id_sr;
    return false;
}

bool Arch::isEnablePort(const PortRef &port) const
{
    if (port.cell == nullptr)
        return false;
    if (isFF(port.cell))
        return port.port == id_e;
    if (port.cell->type == id_icestorm_lc)
        return port.port == id_cen;
    return false;
}


NEXTPNR_NAMESPACE_END
