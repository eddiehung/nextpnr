//#include <cstdio>
//#include <cmath>
//using namespace std;
//
//#include "vtr_util.h"
//#include "vtr_memory.h"
//#include "vtr_log.h"
//
//#include "vpr_types.h"
//#include "vpr_utils.h"
//#include "globals.h"
//#include "path_delay.h"
//#include "path_delay2.h"
//#include "net_delay.h"
//#include "timing_place_lookup.h"
//#include "timing_place.h"
//
//#include "timing_info.h"

static /*vtr*/std::vector</*ClusterNetId,*/ float *> f_timing_place_crit; /* [0..cluster_ctx.clb_nlist.nets().size()-1][1..num_pins-1] */

//static vtr::t_chunk f_timing_place_crit_ch;

/******** prototypes ******************/
static void alloc_crit(/*vtr::t_chunk *chunk_list_ptr*/);

static void free_crit(/*vtr::t_chunk *chunk_list_ptr*/);

/**************************************/

/* Allocates space for the f_timing_place_crit data structure *
* I chunk the data to save space on large problems.           */
static void alloc_crit(/*vtr::t_chunk *chunk_list_ptr*/) {
	auto& cluster_ctx = g_vpr_ctx.clustering();
	float *tmp_ptr;

	f_timing_place_crit.resize(cluster_ctx.clb_nlist.nets().size());

	for (const auto& net : cluster_ctx.clb_nlist.nets()) {
        auto net_id = net.second.get();
		//tmp_ptr = (float *) vtr::chunk_malloc(
		//		(cluster_ctx.clb_nlist.net_sinks(net_id).size()) * sizeof(float), chunk_list_ptr);
		tmp_ptr = (float *) malloc((cluster_ctx.clb_nlist.net_sinks(net_id).size()) * sizeof(float));
		f_timing_place_crit[net_id->udata] = tmp_ptr - 1; /* [1..num_sinks] */
	}
}

///**************************************/
static void free_crit(/*vtr::t_chunk *chunk_list_ptr*/){
//    vtr::free_chunk_memory(chunk_list_ptr);
    for (auto i : f_timing_place_crit)
        free(i+1);
}

/**************************************/
void load_criticalities(/*SetupTimingInfo& timing_info,*/ float crit_exponent /*, const ClusteredPinAtomPinsLookup& pin_lookup*/) {
	/* Performs a 1-to-1 mapping from criticality to f_timing_place_crit.
	  For every pin on every net (or, equivalently, for every tedge ending
	  in that pin), f_timing_place_crit = criticality^(criticality exponent) */

    auto& cluster_ctx = g_vpr_ctx.clustering();
	for (const auto& net : cluster_ctx.clb_nlist.nets()) {
		auto net_id = net.second.get();
		if (cluster_ctx.clb_nlist.net_is_global(net_id))
			continue;

		int ipin = 1;
        for (const auto& clb_pin : cluster_ctx.clb_nlist.net_sinks(net_id)) {
            //int ipin = cluster_ctx.clb_nlist.pin_net_index(clb_pin);

            float clb_pin_crit = calculate_clb_net_pin_criticality(/*timing_info, pin_lookup,*/ clb_pin, net_id);

            /* The placer likes a great deal of contrast between criticalities.
            Since path criticality varies much more than timing, we "sharpen" timing
            criticality by taking it to some power, crit_exponent (between 1 and 8 by default). */
            f_timing_place_crit[net_id->udata][ipin] = pow(clb_pin_crit, crit_exponent);

			++ipin;
        }
	}
}


float get_timing_place_crit(/*ClusterNetId*/ NetInfo* net_id, int ipin) {
    return f_timing_place_crit[net_id->udata][ipin];
}

void set_timing_place_crit(/*ClusterNetId*/ NetInfo* net_id, int ipin, float val) {
    f_timing_place_crit[net_id->udata][ipin] = val;
}

/**************************************/
void alloc_lookups_and_criticalities(/*t_chan_width_dist chan_width_dist,
		t_router_opts router_opts,
		t_det_routing_arch *det_routing_arch, t_segment_inf * segment_inf,
		const t_direct_inf *directs,
		const int num_directs*/) {

//	compute_delay_lookup_tables(router_opts, det_routing_arch, segment_inf,
//			chan_width_dist, directs, num_directs);

	alloc_crit(/*&f_timing_place_crit_ch*/);
}

/**************************************/
void free_lookups_and_criticalities() {
//	//TODO: May need to free f_timing_place_crit ?
	free_crit(/*&f_timing_place_crit_ch*/);

//	free_place_lookup_structs();
}

/**************************************/