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

#ifndef TIMING_H
#define TIMING_H

#include "nextpnr.h"

NEXTPNR_NAMESPACE_BEGIN

// Evenly redistribute the total path slack amongst all sinks on each path
void assign_budget(Context *ctx, bool quiet = false);

// Perform timing analysis and print out the fmax, and optionally the
//    critical path
delay_t timing_analysis(Context *ctx, bool print_fmax = false, bool print_histogram = true, bool print_path = false);

NEXTPNR_NAMESPACE_END

#endif
