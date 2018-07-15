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

#ifndef NEXTPNR_H
#error Include "arch.h" via "nextpnr.h" only.
#endif

NEXTPNR_NAMESPACE_BEGIN

struct ArchArgs
{
};

struct WireInfo;

struct PipInfo
{
    IdString name, bound_net;
    WireId srcWire, dstWire;
    DelayInfo delay;
    std::vector<GraphicElement> graphics;
};

struct WireInfo
{
    IdString name, bound_net;
    std::vector<GraphicElement> graphics;
    std::vector<PipId> downhill, uphill, aliases;
    BelPin uphill_bel_pin;
    std::vector<BelPin> downhill_bel_pins;
    int grid_x, grid_y;
};

struct PinInfo
{
    IdString name;
    WireId wire;
    PortType type;
};

struct BelInfo
{
    IdString name, type, bound_cell;
    std::unordered_map<IdString, PinInfo> pins;
    std::vector<GraphicElement> graphics;
    int grid_x, grid_y;
    bool gb;
};

struct Arch : BaseCtx
{
    std::string chipName;

    std::unordered_map<IdString, WireInfo> wires;
    std::unordered_map<IdString, PipInfo> pips;
    std::unordered_map<IdString, BelInfo> bels;

    std::vector<IdString> bel_ids, wire_ids, pip_ids;
    std::unordered_map<IdString, std::vector<IdString>> bel_ids_by_type;

    std::vector<GraphicElement> frame_graphics;
    float grid_distance_to_delay;

    void addWire(IdString name, int x, int y);
    void addPip(IdString name, IdString srcWire, IdString dstWire, DelayInfo delay);
    void addAlias(IdString name, IdString srcWire, IdString dstWire, DelayInfo delay);

    void addBel(IdString name, IdString type, int x, int y, bool gb);
    void addBelInput(IdString bel, IdString name, IdString wire);
    void addBelOutput(IdString bel, IdString name, IdString wire);
    void addBelInout(IdString bel, IdString name, IdString wire);

    void addFrameGraphic(const GraphicElement &graphic);
    void addWireGraphic(WireId wire, const GraphicElement &graphic);
    void addPipGraphic(PipId pip, const GraphicElement &graphic);
    void addBelGraphic(BelId bel, const GraphicElement &graphic);

    // ---------------------------------------------------------------
    // Common Arch API. Every arch must provide the following methods.

    Arch(ArchArgs args);

    std::string getChipName() { return chipName; }

    IdString archId() const { return id("generic"); }
    IdString archArgsToId(ArchArgs args) const { return id("none"); }

    IdString belTypeToId(BelType type) const { return type; }
    IdString portPinToId(PortPin type) const { return type; }

    BelType belTypeFromId(IdString id) const { return id; }
    PortPin portPinFromId(IdString id) const { return id; }

    BelId getBelByName(IdString name) const;
    IdString getBelName(BelId bel) const;
    uint32_t getBelChecksum(BelId bel) const;
    void bindBel(BelId bel, IdString cell, PlaceStrength strength);
    void unbindBel(BelId bel);
    bool checkBelAvail(BelId bel) const;
    IdString getBoundBelCell(BelId bel) const;
    IdString getConflictingBelCell(BelId bel) const;
    const std::vector<BelId> &getBels() const;
    const std::vector<BelId> &getBelsByType(BelType type) const;
    BelType getBelType(BelId bel) const;
    WireId getWireBelPin(BelId bel, PortPin pin) const;
    BelPin getBelPinUphill(WireId wire) const;
    const std::vector<BelPin> &getBelPinsDownhill(WireId wire) const;

    WireId getWireByName(IdString name) const;
    IdString getWireName(WireId wire) const;
    uint32_t getWireChecksum(WireId wire) const;
    void bindWire(WireId wire, IdString net, PlaceStrength strength);
    void unbindWire(WireId wire);
    bool checkWireAvail(WireId wire) const;
    IdString getBoundWireNet(WireId wire) const;
    IdString getConflictingWireNet(WireId wire) const;
    const std::vector<WireId> &getWires() const;

    PipId getPipByName(IdString name) const;
    IdString getPipName(PipId pip) const;
    uint32_t getPipChecksum(PipId pip) const;
    void bindPip(PipId pip, IdString net, PlaceStrength strength);
    void unbindPip(PipId pip);
    bool checkPipAvail(PipId pip) const;
    IdString getBoundPipNet(PipId pip) const;
    IdString getConflictingPipNet(PipId pip) const;
    const std::vector<PipId> &getPips() const;
    WireId getPipSrcWire(PipId pip) const;
    WireId getPipDstWire(PipId pip) const;
    DelayInfo getPipDelay(PipId pip) const;
    const std::vector<PipId> &getPipsDownhill(WireId wire) const;
    const std::vector<PipId> &getPipsUphill(WireId wire) const;
    const std::vector<PipId> &getWireAliases(WireId wire) const;

    void estimatePosition(BelId bel, int &x, int &y, bool &gb) const;
    delay_t estimateDelay(WireId src, WireId dst) const;
    delay_t getDelayEpsilon() const { return 0.01; }
    delay_t getRipupDelayPenalty() const { return 1.0; }
    float getDelayNS(delay_t v) const { return v; }
    uint32_t getDelayChecksum(delay_t v) const { return 0; }

    const std::vector<GraphicElement> &getFrameGraphics() const;
    const std::vector<GraphicElement> &getBelGraphics(BelId bel) const;
    const std::vector<GraphicElement> &getWireGraphics(WireId wire) const;
    const std::vector<GraphicElement> &getPipGraphics(PipId pip) const;

    bool allGraphicsReload = false;
    bool frameGraphicsReload = false;
    std::unordered_set<BelId> belGraphicsReload;
    std::unordered_set<WireId> wireGraphicsReload;
    std::unordered_set<PipId> pipGraphicsReload;

    bool getCellDelay(const CellInfo *cell, IdString fromPort, IdString toPort, delay_t &delay) const;
    IdString getPortClock(const CellInfo *cell, IdString port) const;
    bool isClockPort(const CellInfo *cell, IdString port) const;
    // Return true if a port is a net
    bool isGlobalNet(const NetInfo *net) const;
    // Return true if cell is a IO
    bool isIO(const CellInfo* cell) const;

    bool isValidBelForCell(CellInfo *cell, BelId bel) const;
    bool isBelLocationValid(BelId bel) const;
};

NEXTPNR_NAMESPACE_END
