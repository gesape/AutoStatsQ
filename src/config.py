from pyrocko.guts import Object, Float, Int, String, Bool, List, Tuple, Dict
from pyrocko.gf import TPDef


class GeneralSettings(Object):
    work_dir = String.T(
        help='Upper data directory')
    list_station_lists = List.T(
        help='List of station list names, either txt or xml files')
    st_white_list = List.T(default=[])

class CatalogConfig(Object):
    search_events = Bool.T()
    use_local_catalog = Bool.T()
    subset_of_local_catalog = Bool.T()
    use_local_subsets = Bool.T()
    subset_fns = Dict.T()
    # catalog
    catalog_fn = String.T(optional=True,
        help='Name of local catalog to save to or to read from')
    dist = Float.T(default=165.0,
        help='Max. distance to show on catalog plot (map) [deg]')
    min_mag = Float.T()
    max_mag = Float.T()
    tmin_str = String.T()
    tmax_str = String.T()
    wedges_width = Int.T(
        help='Find one event in backazimuthal wedge of this width')
    mid_point = List.T(optional=True,
        help='Enter some estimation on centre point of array, needed for\
             subset based on distances and plotting only.')
    median_ev_in_bin = Bool.T(default=True)
    # weighted_magn_baz_ev = Bool.T()

    min_dist_km = Float.T()
    max_dist_km = Float.T()
    depth_options = Dict.T()

    plot_catalog_all = Bool.T(default=False)
    plot_hist_wedges = Bool.T(default=False)
    plot_wedges_vs_dist = Bool.T(default=False)
    plot_wedges_vs_magn = Bool.T(default=False)
    plot_dist_vs_magn = Bool.T(default=False)
    plot_catalog_subset = Bool.T(default=False)
    # pl_opt = List.T(
    # help='[map center lat, centre lon, map radius [m]]')


class ArrTConfig(Object):
    calc_first_arr_t = Bool.T(
        help='calculate first arrival times using cake?')
    phase_select = String.T()
    calc_est_R = Bool.T(default=False)
    v_rayleigh = Float.T(default=4.0)


#class CakeTTTGenerator(Object):
#    calc_ttt = Bool.T(default=False)
#    dir_ttt = String.T(optional=True)
#    earthmodel_id = String.T(optional=True)
#    tabulated_phases = List.T(
#        TPDef.T(), optional=True)
#    dist_min = Float.T(optional=True)
#    dist_max = Float.T(optional=True)
#    dist_acc = Float.T(optional=True)
#    s_depth_min = Float.T(optional=True)
#    s_depth_max = Float.T(optional=True)
#    s_acc = Float.T(optional=True)
#    r_depth_min = Float.T(optional=True)
#    r_depth_max = Float.T(optional=True)
#    r_acc = Float.T(optional=True)
#    t_acc = Float.T(optional=True)


class MetaDataDownloadConfig(Object):
    # Data and Metadata download
    download_data = Bool.T()
    download_metadata = Bool.T()
    local_metadata = List.T(default=[])
    local_data = List.T(default=[])
    local_waveforms_only = Bool.T(default=False)
    sds_structure = Bool.T(default=False)
    use_downmeta = Bool.T(default=True)
    channels_download = String.T()
    all_channels = Bool.T(default=False)
    # components = List.T()
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
    deltat_down = Float.T(default=2.0)


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
    snr_thresh = Float.T(default=1.)
    debug_mode = Bool.T(default=False)
    phase_select = String.T()
    components = List.T()

    plot_median_gain_on_map = Bool.T()
    plot_allgains = Bool.T()


class PSDConfig(Object):
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
    ccmin = Float.T(optional=True)
    plot_distr = Bool.T()
    plot_orient_map_fromfile = Bool.T()
    plot_angles_vs_events = Bool.T(default=False)
    orient_map_label = List.T(optional=True)
    debug_mode = Bool.T(default=False)

class TimingConfig(Object):
    timing_test = Bool.T()
    bandpass = Tuple.T(3, Float.T())
    time_wdw = List.T()
    cc_thresh = Float.T()
    search_locations = Bool.T()
    debug_mode = Bool.T(default=False)


class TeleCheckConfig(Object):
    tele_check = Bool.T()

class maps(Object):
    map_size = List.T()
    pl_opt = List.T()
    pl_topo = Bool.T(default=False)


class AutoStatsQConfig(Object):
    Settings = List.T()
