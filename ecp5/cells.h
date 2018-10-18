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

#ifndef ECP5_CELLS_H
#define ECP5_CELLS_H

#include "nextpnr.h"

NEXTPNR_NAMESPACE_BEGIN

// Create a standard ECP5 cell and return it
// Name will be automatically assigned if not specified
std::unique_ptr<CellInfo> create_ecp5_cell(Context *ctx, IdString type, std::string name = "");

// Return true if a cell is a LUT
inline bool is_lut(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("LUT4"); }

// Return true if a cell is a flipflop
inline bool is_ff(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("TRELLIS_FF"); }

inline bool is_carry(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("CCU2C"); }

inline bool is_lc(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("TRELLIS_LC"); }

inline bool is_trellis_io(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("TRELLIS_IO"); }

inline bool is_dpram(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("TRELLIS_DPR16X4"); }

inline bool is_pfumx(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("PFUMX"); }

inline bool is_l6mux(const BaseCtx *ctx, const CellInfo *cell) { return cell->type == ctx->id("L6MUX21"); }

void ff_to_slice(Context *ctx, CellInfo *ff, CellInfo *lc, int index, bool driven_by_lut);
void lut_to_slice(Context *ctx, CellInfo *lut, CellInfo *lc, int index);
void ccu2c_to_slice(Context *ctx, CellInfo *ccu, CellInfo *lc);
void dram_to_ramw(Context *ctx, CellInfo *ram, CellInfo *lc);
void dram_to_ram_slice(Context *ctx, CellInfo *ram, CellInfo *lc, CellInfo *ramw, int index);

NEXTPNR_NAMESPACE_END

#endif
