import os
from .config import GeneralSettings, CatalogConfig, ArrTConfig,\
MetaDataDownloadConfig, RestDownRotConfig, SynthDataConfig,\
GainfactorsConfig, PSDConfig, OrientConfig, AutoStatsQConfig,\
TimingConfig, maps, TeleCheckConfig


def generate_default_config():
    gensettings = GeneralSettings(
        work_dir='/some/data/directory/',
        list_station_lists=['/path/to/station-file/file.csv',
                            '/path/to/station-file/file.xml'])

    catalogconf = CatalogConfig(
                                search_events=True,
                                use_local_catalog=False,
                                subset_of_local_catalog=True,
                                use_local_subsets=False,
                                subset_fns={},
                                # Filename of catalog (local or name for saving catalog)
                                # catalog_fn='catalog.txt',
                                # catalog plot (backazimuthal view) max. distance to show
                                # dist=165.,
                                min_mag=6.5,
                                max_mag=8.5,
                                tmin_str='2000-01-01 00:00:00',
                                tmax_str='2019-10-01 00:00:00',
                                # backazimuthal distribution
                                wedges_width=15,  # deg
                                # enter a midpoint of array for global event maps
                                # mid_point=[46.98, 10.74],
                                # enter lat, lon, radius for output maps showing gains,
                                # polarisation, etc.
                                min_dist_km=1000.,
                                max_dist_km=9999999.9,
                                depth_options={'deep': [25000, 1000000],
                                               'shallow': [100, 40000]},
                                plot_catalog_all=False,
                                plot_hist_wedges=False,
                                plot_catalog_subset=False)
                                #weighted_magn_baz_ev=True)

    arrTconf = ArrTConfig(calc_first_arr_t=False,
                          # The first arriving phase will
                          # be used for the time window.
                          # Type 'cake list-phase-map' into terminal to find more
                          # phase names for cake.
                          phase_select='P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                       'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p',
                          calc_est_R=False)

    #cake_ttt_gen = CakeTTTGenerator(calc_ttt=False)

    metaDataconf = MetaDataDownloadConfig(
        download_data=False,
        download_metadata=False,
        channels_download='HH*',
        #local_metadata=[],
        use_downmeta=False,
        token={'geofon': os.path.join(os.path.expanduser('~'), '.eidatoken')},
        sites=['geofon', 'orfeus', 'iris'],
        dt_start=0.1,
        dt_end=1.5)

    RestDownconf = RestDownRotConfig(rest_data=False,
                                     freqlim=(0.005, 0.01, 0.2, 0.25),
                                     rotate_data=False,
                                     deltat_down=2.)

    synthsconf = SynthDataConfig(make_syn_data=False,
                                 engine_path='/path/to/GF_stores',
                                 store_id='global_2s')

    gainfconf = GainfactorsConfig(calc_gainfactors=False,
                                  gain_factor_method=['reference_nsl', ('NET', 'STAT')],#['median_all_avail'],
                                  fband={'order': 4, 'corner_hp': 0.01, 'corner_lp': 0.2},
                                  taper_xfrac=0.25,
                                  wdw_st_arr=5,
                                  wdw_sp_arr=120,
                                  # phase_select='first(P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                  #              'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p)',
                                  components=['Z'],
                                  plot_median_gain_on_map=False,
                                  plot_allgains=False,
                                  debug_mode=False)

    psdsconf = PSDConfig(calc_psd=False,
                         tinc=600,
                         tpad=200,
                         dt_start=60,
                         dt_end=120,
                         n_poly=25,
                         norm_factor=50,
                         f_ign=0.02,
                         plot_psds=False,
                         plot_ratio_extra=False,
                         plot_m_rat=False,
                         plot_flat_ranges=False)

    orientconf = OrientConfig(orient_rayl=False,
                              bandpass=(3, 0.01, 0.05),
                              start_before_ev=30,
                              stop_after_ev=480,
                              ccmin=0.80,
                              plot_heatmap=False,
                              plot_distr=False,
                              plot_orient_map_fromfile=False,
                              plot_angles_vs_events=False,
                              debug_mode=False)

    timingconf = TimingConfig(timing_test=False,
                              bandpass=(3, 0.01, 0.1),
                              time_wdw=['firstP', 1200],
                              cc_thresh=0.6,
                              # search_avail_stats=True,
                              debug_mode=False)

    _tele_check = TeleCheckConfig(tele_check=False)

    
    _maps = maps(map_size=[30.0, 30.0],
                 pl_opt=['lat', 'lon', 'radius', 'split'],
                 pl_topo=False)


    config = AutoStatsQConfig(
      Settings=[gensettings, catalogconf, arrTconf, metaDataconf,
                RestDownconf,
                synthsconf, gainfconf, psdsconf, orientconf, timingconf, 
                _tele_check, _maps])

    return config
