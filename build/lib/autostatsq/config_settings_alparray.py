from config import GeneralSettings, CatalogConfig, ArrTConfig,\
MetaDataDownloadConfig, RestDownRotConfig, SynthDataConfig,\
GainfactorsConfig, PSDConfig, OrientConfig, AutoStatsQConfig

def generate_default_config():
    gensettings = GeneralSettings(
        data_dir='/media/gesap/TOSHIBA EXT/data/teleseismics7/',
        list_station_lists=['/home/gesap/Documents/AlpArray/station_net/SwathD.xml', 
        '/home/gesap/Documents/AlpArray/station_net/\
AlpArray_permanent_stations_01_Sept2016.csv'])
        # '/home/gesap/Documents/AlpArray/station_net/\
# AlpArray_temporary_stations_01_Sept2016.csv'

    catalogconf = CatalogConfig(
                                # Filename of catalog (local or name for saving catalog)
                                catalog_fn='global_2015-2018_Mgr6,5.txt',
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
                                use_local_catalog=False,
                                subset_of_local_catalog=False,
                                use_local_subsets=True,
                                subset_fns={'deep': '/media/gesap/TOSHIBA EXT/data/teleseismics10/results/catalog/catalog_Mgr6.5_deep_geofon.txt',
                                             'shallow': '/media/gesap/TOSHIBA EXT/data/teleseismics10/results/catalog/catalog_Mgr6.5_shallow_geofon.txt'},
                                min_dist_km=1000.,
                                max_dist_km=9999999.9,
                                plot_catalog_all=False,
                                plot_hist_wedges=False,
                                plot_wedges_vs_dist=False,
                                plot_wedges_vs_magn=False,
                                plot_dist_vs_magn=False,
                                plot_catalog_subset=False,
                                median_ev_in_bin=False,
                                weighted_magn_baz_ev=True)

    arrTconf = ArrTConfig(calc_first_arr_t=False,
                          # The first arriving phase will
                          # be used for the time window.
                          # Type 'cake list-phase-map' into terminal to find more
                          # phase names for cake.
                          phase_select = 'P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                        'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p')

    metaDataconf = MetaDataDownloadConfig(
        download_data=False,
        download_metadata=False,
        components_download='HH*',
        components=['HHZ', 'HHN', 'HHE'],
        token='/home/gesap/Documents/AlpArray/download-waveforms/swathD/token.asc',
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
                                  gain_factor_method='median_all_avail',
                                  fband={'order': 4, 'corner_hp': 0.01, 'corner_lp': 0.2},
                                  taper_xfrac=0.25,  # richtig?
                                  wdw_static_length=60.,
                                  wdw_phase_position=0.5,
                                  phase_select = 'first(P|p|P(cmb)P(icb)P(icb)p(cmb)p|' +\
                                                 'P(cmb)Pv(icb)p(cmb)p|P(cmb)P<(icb)(cmb)p)',
                                  plot_mean_gain_on_map=False,
                                  map_gain_size=(30.,30.),
                                  gain_rel_syn=False)

    psdsconf = PSDConfig(calc_psd=False,
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
      Settings=[gensettings, catalogconf, arrTconf, metaDataconf, RestDownconf,
                synthsconf, gainfconf, psdsconf, orientconf])

    return config