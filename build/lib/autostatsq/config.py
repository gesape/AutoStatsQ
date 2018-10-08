from pyrocko.guts import Object, Float, Int, String, Bool, List, Tuple, Dict
import os
from pyrocko.gf import TPDef


class GeneralSettings(Object):
	data_dir = String.T(
		help='Upper data directory')
	list_station_lists = List.T(
		help='List of station list names, either txt or xml files')

class CatalogConfig(Object):
	search_events = Bool.T()
	use_local_catalog = Bool.T()
	subset_of_local_catalog = Bool.T()
	use_local_subsets = Bool.T()
	subset_fns = Dict.T()
	# catalog
	catalog_fn = String.T(
		help='Name of local catalog to save to or to read from')
	dist = Float.T(
		help='Max. distance to show on catalog plot (map) [deg]')
	min_mag = Float.T()
	max_mag=Float.T()
	tmin_str = String.T()
	tmax_str = String.T()
	wedges_width = Int.T(
		help='Find one event in backazimuthal wedge of this width')
	mid_point = List.T(
		help='Enter some estimation on centre point of array, needed for\
             subset based on distances and plotting only.')
	median_ev_in_bin = Bool.T(default=True)
	#weighted_magn_baz_ev = Bool.T()

	min_dist_km = Float.T()
	max_dist_km = Float.T()
	depth_options = Dict.T()  

	plot_catalog_all = Bool.T()
	plot_hist_wedges = Bool.T()
	plot_wedges_vs_dist = Bool.T()
	plot_wedges_vs_magn = Bool.T()
	plot_dist_vs_magn = Bool.T()
	plot_catalog_subset = Bool.T()
	#pl_opt = List.T(
#		help='[map center lat, centre lon, map radius [m]]')



class ArrTConfig(Object):
	calc_first_arr_t = Bool.T(
		help='calculate first arrival times using cake?')
	phase_select = String.T()
	calc_est_R = Bool.T(default=False)


class CakeTTTGenerator(Object):
	calc_ttt=Bool.T(
		help='calculate first arrival times using ttt interpolation?')
	dir_ttt=String.T()
	earthmodel_id = String.T()
	tabulated_phases = List.T(
		TPDef.T(),
		help='list of tabulated phase definitions usable shifters')
	dist_min = Float.T()
	dist_max = Float.T()
	dist_acc = Float.T()
	s_depth_min = Float.T()
	s_depth_max = Float.T()
	s_acc = Float.T()
	r_depth_min = Float.T(default=0.)
	r_depth_max = Float.T(default=0.)
	r_acc = Float.T()
	t_acc = Float.T()


class MetaDataDownloadConfig(Object):
	# Data and Metadata download
	download_data = Bool.T()
	download_metadata = Bool.T()
	local_metadata = List.T(default=[])
	use_downmeta = Bool.T(default=True)
	components_download = String.T()
	#components = List.T()
	token = Dict.T(default={})
	sites = List.T()
	# start and end time with respect to origin time
	# start time before origin time!
	dt_start = Float.T()
	dt_end = Float.T()


class RestDownRotConfig(Object):
	# Restitution
	rest_data = Bool.T()
	freqlim = Tuple.T(4, Float.T())

	# Downsampling & Rotate NE --> RT
	rotate_data = Bool.T()
	deltat_down = Int.T()


class SynthDataConfig(Object):
	# syn. data
	make_syn_data = Bool.T()
	engine_path = String.T()
	store_id = String.T()


class GainfactorsConfig(Object):
	# Gainfactors
	calc_gainfactors = Bool.T()
	gain_factor_method = List.T()
	''' gain relative to options:
	   * 'scale_one' for reference amplitude = 1.
	   *  tuple ('reference_nsl', refernce_id) - gain
	      relative to one specific station,
	      reference_id e.g. 'BFO'
	   * 'median_all_avail': first assesses all abs. amplitudes,
	      in the end chooses a reference station of those that recorded
	      most events based on median gain of first event
	'''
	fband = Dict.T()
	taper_xfrac = Float.T()
	wdw_st_arr = Int.T()
	wdw_sp_arr = Int.T()
	phase_select = String.T()
	components = List.T()

	plot_median_gain_on_map = Bool.T()
	#map_gain_size = Tuple.T(2, Float.T())
	plot_allgains = Bool.T()


class PSDConfig(Object):
	# PSD section
	calc_psd = Bool.T()
	tinc = Int.T()
	tpad = Int.T()
	dt_start = Int.T()
	dt_end = Int.T()
	n_poly = Int.T(default=25)
	norm_factor = Int.T(default=10)
	f_ign = Float.T(default=0.02)
	only_first = Bool.T(default=True)
	plot_psds = Bool.T()
	plot_ratio_extra = Bool.T(default=False)
	plot_m_rat = Bool.T(default=False)
	plot_flat_ranges = Bool.T()
	plt_neigh_ranges = Bool.T(default=False)


class OrientConfig(Object):
	# orientations
	orient_rayl = Bool.T()
	bandpass = Tuple.T(3, Float.T())
	start_before_ev = Float.T()
	stop_after_ev = Float.T()
	plot_heatmap = Bool.T()
	plot_distr = Bool.T()
	plot_orient_map_fromfile = Bool.T()
	orient_map_label = List.T(optional=True)


class maps(Object):
    map_size = List.T()
    pl_opt = List.T()
    pl_topo = Bool.T()


class AutoStatsQConfig(Object):
	Settings = List.T()
