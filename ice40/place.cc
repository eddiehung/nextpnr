/*
 *  nextpnr -- Next Generation Place and Route
 *
 *  Copyright (C) 2018  Clifford Wolf <clifford@symbioticeda.com>
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

#include "design_utils.h"
#include "log.h"
#include "util.h"

NEXTPNR_NAMESPACE_BEGIN

void place_gbs(Context *ctx)
{
    std::vector<BelId> gb_reset;
    std::vector<BelId> gb_cen;

    for (auto bel : ctx->getBels()) {
        BelType type = ctx->getBelType(bel);
        if (type == TYPE_SB_GB) {
            IdString glb_net = ctx->getWireName(ctx->getWireBelPin(bel, PIN_GLOBAL_BUFFER_OUTPUT));
            int glb_id = std::stoi(std::string("") + glb_net.str(ctx).back());
            if (glb_id % 2 == 0)
                gb_reset.push_back(bel);
            else
                gb_cen.push_back(bel);
        }
    }
    for (auto &c : ctx->cells) {
        CellInfo *cell = c.second.get();
        if (cell->type == ctx->id_sb_gb) {
            bool is_reset = false, is_cen = false;
            NPNR_ASSERT(cell->ports.at(ctx->id_glb_buf_out).net != nullptr);
            for (auto user : cell->ports.at(ctx->id_glb_buf_out).net->users) {
                if (ctx->isResetPort(user))
                    is_reset = true;
                if (ctx->isEnablePort(user))
                    is_cen = true;
            }
            NPNR_ASSERT(!is_reset || !is_cen);
            if (is_reset) {
                ctx->bindBel(gb_reset.back(), cell->name, STRENGTH_WEAK);
                gb_reset.pop_back();
            } 
            else if (is_cen) {
                ctx->bindBel(gb_cen.back(), cell->name, STRENGTH_WEAK);
                gb_cen.pop_back();
            }
        }
    }
}

NEXTPNR_NAMESPACE_END
