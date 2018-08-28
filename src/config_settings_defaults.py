from .config import GeneralSettings, CatalogConfig, ArrTConfig,\
CakeTTTGenerator,\
MetaDataDownloadConfig, RestDownRotConfig, SynthDataConfig,\
GainfactorsConfig, PSDConfig, OrientConfig, AutoStatsQConfig
from pyrocko.gf import TPDef


def generate_default_config():
    gensettings = GeneralSettings(
        data_dir='/media/gesap/TOSHIBA EXT/data/teleseismics9/',
        list_station_lists=['/home/gesap/Documents/AlpArray/station_net/\
AlpArray_permanent_stations_01_Sept2016.csv',
'/home/gesap/Documents/AlpArray/station_net/SwathD.xml'])
        # '/home/gesap/Documents/AlpArray/station_net/\
# AlpArray_temporary_stations_01_Sept2016.csv'

    catalogconf = CatalogConfig(
                                # Filename of catalog (local or name for saving catalog)
                                catalog_fn='catalog_Mgr6.5.txt',
                                # catalog plot (backazimuthal view) max. distance to show
                                dist=165.,
                                min_mag=6.5,
                                max_mag=8.5,
                                tmin_str='2015-01-01 00:00:00',
                                tmax_str='2018-03-02 00:00:00',
                                # backazimuthal distribution
                                wedges_width=15,  # deg
                                # enter a midpoint of array for global event maps
                                mid_point=[46.98, 10.74],
                                # enter lat, lon, radius for output maps showing gains,
                                # polarisation, etc.
                                pl_opt = [46, 12, 700000],
                                search_events=False,
                                use_local_catalog=True,
                                subset_of_local_catalog=True,
                                min_dist_km=1000.,
                                max_dist_km=9999999.9,
                                depth_options=[['deep', 60, 1000], ['shallow', 1, 60]],  
                                plot_catalog_all=False,
                                plot_hist_wedges=False,
                                plot_wedges_vs_dist=False,
                                plot_wedges_vs_magn=False,
                                plot_dist_vs_magn=False,
                                plot_catalog_subset=False,
                                median_ev_in_bin=False)#,
                                #weighted_magn_baz_ev=True)

    arrTconf = ArrTConfig(calc_first_arr_t=False,
                          # The first arriving phase will
                          # be used for the time window.
                          # Type 'cake list-phase-map' into terminal to find more
                          # phase names for cake.
                          phase_select = 'P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                        'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p')

    cake_ttt_gen = CakeTTTGenerator(calc_ttt=False,
                          earthmodel_id='prem-no-ocean.f',
                          dir_ttt='/media/gesap/TOSHIBA EXT/data/teleseismics9/ttt/',
                          tabulated_phases = [TPDef(id='p', definition='P,p')],
                          dist_min = 1000.,
                          dist_max = 9999.,
                          dist_acc = 2.,
                          s_depth_min = 10.,
                          s_depth_max = 50.,
                          s_acc = 2.,
                          r_depth_min = 0.,
                          r_depth_max = 0.,
                          r_acc = 2.,
                          t_acc = 3.,
                          )

    metaDataconf = MetaDataDownloadConfig(
        download_data=False,
        download_metadata=False,
        components_download='HH*',
        components=['HHZ', 'HHN', 'HHE'],
        token={'geofon': '/home/gesap/Documents/AlpArray/download-waveforms/swathD/token.asc'},
        sites=['geofon', 'orfeus', 'iris', 'bgr', 'ingv', 'lmu'],
        dt_start=0.1,
        dt_end=1.5)

    RestDownconf = RestDownRotConfig(rest_data=False,
                                     freqlim=(0.005, 0.01, 0.2, 0.25),
                                     rotate_data=False,
                                     deltat_down=2)

    synthsconf = SynthDataConfig(make_syn_data=False,
                                 engine_path='/media/gesap/TOSHIBA EXT/gf_stores',
                                 store_id='global_2s')

    gainfconf = GainfactorsConfig(calc_gainfactors=False,
                                  gain_factor_method=['reference_nsl', ('GE', 'MATE')],#['median_all_avail'],
                                  fband={'order': 4, 'corner_hp': 0.01, 'corner_lp': 0.2},
                                  taper_xfrac=0.25,  # richtig?
                                  wdw_st_arr=5,
                                  wdw_sp_arr=60,
                                  phase_select = 'first(P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                                 'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p)',
                                  components=['Z','R','T'],
                                  plot_mean_gain_on_map=False,
                                  map_gain_size=(30.,30.),
                                  plot_allgains=False)

    psdsconf = PSDConfig(calc_psd=False,
    					 tinc=600,
    					 tpad=200,
    					 dt_start=60,
    					 dt_end=1800,
    					 n_poly = 25,
	                     norm_factor = 50,
	                     f_ign = 0.02,
                         plot_psds=False,
                         plot_ratio_extra=False,
                         plot_m_rat=False,
                         plot_flat_ranges=False,
                         plt_neigh_ranges=False)

    orientconf = OrientConfig(orient_rayl=False,
    	                      bandpass=(3, 0.01, 0.05),
    	                      start_before_ev=30,
    	                      stop_after_ev=480,
                              plot_heatmap=False,
                              plot_distr=False)

    config = AutoStatsQConfig(
      Settings=[gensettings, catalogconf, arrTconf, cake_ttt_gen, metaDataconf,
                RestDownconf,
                synthsconf, gainfconf, psdsconf, orientconf])

    return config