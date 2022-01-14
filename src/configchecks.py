import sys, logging, os
from .config import GeneralSettings, CatalogConfig, ArrTConfig,\
MetaDataDownloadConfig, RestDownRotConfig, SynthDataConfig,\
GainfactorsConfig, PSDConfig, OrientConfig, TimingConfig, TeleCheckConfig,\
maps, AutoStatsQConfig


logs = logging.getLogger('CONFIG')
logger = logs

def check_general_settings(gensettings):

    error = False
    if os.path.exists(gensettings.work_dir):
        logger.debug(' Check passed: Working directory exists.')
    else:
        logger.error('ERROR: Working directory not found: %s.'
                     % gensettings.work_dir)
        error = True

    if len(gensettings.list_station_lists) < 1:
        logger.error('ERROR: List of station lists is empty. Please enter filenames.')
        error = True
    else:
        for f in gensettings.list_station_lists:
            if os.path.exists(f):
                logger.debug(' Check passed: Found station file %s.' % f)
            else:
                logger.error('ERROR: Station file not found: %s.' % f)
                error = True
    
    if len(gensettings.st_use_list) > 0:                      
        logger.debug('Using only stations: %s' % st_use_list)

    return error


def check_catalog_settings(catalogconf):
    cf = catalogconf

    error = False

    if not cf.search_events:
        if not (cf.use_local_catalog or cf.use_local_subsets):
            logger.error(' ERROR: Either use option search_events or provide a local catalog (use_local_catalog and catalog_fn) '
                         + ' or option use_local_subsets.')
            error = True

        if cf.use_local_catalog and cf.use_local_subsets:
            logger.error(' ERROR: Either use option use_local_catalog  '
                         + ' or option use_local_subsets.')
            error = True

        if cf.subset_of_local_catalog and cf.use_local_subsets:
            logger.error(' ERROR: Either use option use_local_catalog and subset_of_local_catalog  '
                         + ' or option use_local_subsets.')
            error = True
    else:
        if cf.use_local_catalog or cf.use_local_subsets:
            logger.error(' ERROR: Either use option search_events or provide a local catalog (use_local_catalog and catalog_fn) '
                         + ' or option use_local_subsets and subset_fns.')
            error = True

    if cf.subset_fns:
        for key, value in cf.subset_fns.items():
            if os.path.exists(value):
                logger.debug(' Check passed: %s catalog file exists.' % key)
            else:
                logger.error(' ERROR: %s catalog file not found: %s.' % (key, value))
                error = True
    
    
    if cf.catalog_fn:
        if os.path.exists(cf.catalog_fn):
            logger.debug(' Check passed: Catalog file exists.')
        else:
            logger.error(' ERROR: Catalog file not found: %s.'
                         % cf.catalog_fn)
            error = True

    # select events
    if cf.min_mag < 6.0:
        logger.warning(' Events with smaller magnitudes than Mw 6.0 are unlikely'
                       + ' to have distinct P and Rayleigh waves needed for the test.'
                       + ' We recommend using Mw 6.5-7.5.')
    if cf.max_mag > 8.0:
        logger.warning(' Events with larger magnitudes than Mw 8.0 may be very complex'
                       + ' and therefore may not be forward modeled well for comparison.'
                       + ' We recommend using Mw 6.5-7.5.')

    logger.debug(' t min for catalog search: %s' % cf.tmin_str)
    logger.debug(' t max for catalog search: %s' % cf.tmax_str)

    if cf.min_dist_km < 1000:
        logger.info(' If possible, the events should be at a distance >> than the inter-station distances.')
    
    if cf.min_dist_km > cf.max_dist_km:
        logger.error(' min distance larger max distance.')
        error = True

    if cf.depth_options['shallow']:
        logger.debug(' Catalog subset depth, shallow events: %s-%s km.'
                      % (cf.depth_options['shallow'][0]/1000, 
                        cf.depth_options['shallow'][1]/1000))

        if cf.depth_options['shallow'][0] > cf.depth_options['shallow'][1]:
            logger.error(' ERROR: depth_options: shallow: min depth larger max depth.')
            error = True

    if cf.depth_options['deep']:
        logger.debug(' Catalog subset depth, deep events: %s-%s km.'
                      % (cf.depth_options['deep'][0]/1000, 
                        cf.depth_options['deep'][1]/1000))

        if cf.depth_options['deep'][0] > cf.depth_options['deep'][1]:
            logger.error(' ERROR: depth_options: deep: min depth larger max depth.')
            error = True

    if cf.depth_options['shallow'] and cf.depth_options['deep']:
        if cf.depth_options['shallow'][1] > cf.depth_options['deep'][1]:
            logger.error(' ERROR: depth_options: deep settings shallower than shallow settings.')
            error = True

    if cf.wedges_width > 15:
        logger.debug(' The wedges width defines in which azimuthal bins events are found.'
                     + ' A too large bin results in less events. Set to smaller values to increase number of events')

    if cf.median_ev_in_bin:
        logger.debug(' Using median event per azimuthal bin.')

    if not cf.mid_point:
        logger.debug(' Mid point of station array/network for subsets and of maps for catalog plots determined autmatically.')
    else:
        if not len(cf.mid_point) == 2:
            logger.error(' Mid point of station array/network for subsets and of maps for catalog plots must have length=2 (lat, lon) or'
                         + ' can be left empty for internal computation.')
            error = True

    if cf.dist > 165.0 and (cf.plot_catalog_all or cf.plot_catalog_subset):
        logger.error(' Plotting distance dist [deg] for azimuthal plots cannot be larger than 165 degrees')
        error = True
    
    logger.debug(' Plotting distance dist [deg] for azimuthal plots set to %s degrees.' % cf.dist)

    return error


def check_arrivaltime_settings(arrTconf):
    pass

def check_metadata_settings(metaDataconf):
    cf = metaDataconf
    error = False

    if not cf.use_downmeta:
        logger.warning(' Set MetaDataDownloadConfig use_downmeta to true to use the downloaded metadata.')
    
    if len(cf.channels_download) !=3:
        logger.error(' MetaDataDownloadConfig channels_download takes a string with three characters as input. Use * as wildcard, e.g. HH*.')
        error = True

    if cf.dt_start > 0.5 or cf.dt_start < 0 or cf.dt_end < 0.5:
        logger.warning( 'MetaDataDownloadConfig dt_start is in hours before and dt_end in hours after origin time. Make sure time windows are long enough to include body waves, surface waves and sufficient time for filtering etc.')
    
    return error


def check_restitution_settings(RestDownconf):
    cf = RestDownconf

    if not cf.freqlim[0] < cf.freqlim[1] < cf.freqlim[2] < cf.freqlim[3]:
        logging.error(' RestDownRotConfig frelim is defining the frequency range (corner frequencies) for which data is restituted. First entry must be smallest, last largest.')
    pass

def check_synthetics_settings(synthsconf):
    cf = synthsconf
    error = False

    if cf.make_syn_data and cf.engine_path and cf.store_id:
        if os.path.exists(os.path.join(cf.engine_path, cf.store_id)):
            logger.debug(' Check passed: synthsconf Engine path and gfdb exists.')
        else:
            logger.error( 'ERROR:  SynthDataConfig Engine path and/or gfdb not found.')
            error = True

    if cf.make_syn_data and not (cf.engine_path and cf.store_id):
        logger.error(' ERROR: SynthDataConfig To compute synthetic data, the engine path amd gfdb must be defined.')
        error = True

    return error


def check_gain_settings(gainfconf, RestDownconf):
    cf = gainfconf
    error = False

    # check methods
    if len(cf.gain_factor_method) == 1:
        if cf.gain_factor_method[0] not in ['scale_one', 'median_all_avail', 'syn_comp']:
            logging.info(cf.gain_factor_method)
            logging.error(' Gain factor method not known. Choose between reference_nsl, syn_comp, median_all_avail, scale_one.')
            error = True

    elif len(cf.gain_factor_method) == 2:
        if cf.gain_factor_method[0] != 'reference_nsl':
            logging.error(' Gain factor method not known. Choose between reference_nsl, syn_comp, median_all_avail, scale_one.')
            error = True
    else:
        logging.error(' Gain factor method not known. Choose between reference_nsl, syn_comp, median_all_avail, scale_one.')
        error = True

    # check len twd wrt frequency
    if cf.fband['corner_hp'] > cf.fband['corner_lp']:
        logging.error(' GainfactorsConfig HP corner of BP filter must be smaller than LP corner.')
        error = True

    if cf.fband['corner_lp'] > RestDownconf.freqlim[2]:
        logging.error(' GainfactorsConfig LP corner of BP filter can be at maximum be same as limit of restitution frequency limit (3rd entry of RestDownRotConfig.freqlim.')
        error = True

    # component Z
    for x in cf.components:
        if x not in ['Z', 'N', 'E']:
            logging.error(' GainfactorsConfig Components must Z, N and/or E, e.g. [Z] or [Z,N,E].')
            error = True

    return error


def check_psd_settings(psdsconf):
    # twd long enough?
    # explain tinc and dt_start etc

    # n_poly large enough, not too large
    # norm factor?
    # f_ign 
    # --> all a simple waring if far from default values, but okay to test

    pass

def check_orient_settings(orientconf):
    cf = orientconf
    error = False

    if isinstance(cf.bandpass[0], int):
        logging.error(' First part of OrientConfig bandpass is steepness of filter (integer).')
        error = True

    if cf.bandpass[1] > cf.bandpass[2]:
        logging.error(' OrientConfig bandpass has form [filter_steepness, fmin, fmax]. fmin must be smaller fmax.')
        error = True

    if cf.bandpass[2] > 0.06:
        logging.warning(' OrientConfig bandpass should be chosen to emphasize surface waves. Usually fmax is best chosen smaller than 0.05.')
    
    if cf.stop_after_ev < 400:
        logging.warning(' OrientConfig stop_after_ev must enable a long time window to do the test with deep frequencies without filer effects. Please choose 600 s if possible.')
    
    if cf.ccmin < 0.7:
        logging.warning(' OrientConfig Make sure to set a sufficiently large cross-correlation threshold to get meaningful results. Test e.g. using 0.7, 0.75, 0.8 or even higher.')
    
    return error


def check_timing_settings(timingconf):
    error = False
    cf = timingconf

    # bandpass
    # twd
    # cc thresh large enough

    if isinstance(cf.bandpass[0], int):
        logging.error(' First part of TimingConfig bandpass is steepness of filter (integer).')
        error = True

    if cf.bandpass[1] > cf.bandpass[2]:
        logging.error(' TimingConfig bandpass has form [filter_steepness, fmin, fmax]. fmin must be smaller fmax.')
        error = True

    return error


def check_telecheck_settings(tc):
    pass


def check_maps_settings(maps):
    error = False

    # copy test from main script
    if maps.pl_opt == ['automatic']:
        logging.debug(' Map plotting settings determined automatically.')

    elif len(maps.pl_opt) == 4:
        for o in maps.pl_opt[0:2]:
            if isinstance(o,str):
                logs.error(' Please make sure the %s in the map plotting option in the confic files is set.' % o)
                
    else:
        logs.error(' Set pl_opt to *automatic* or provide [lat, lon, radius, cscale].')
        error = True

    if maps.outformat not in ['png', 'pdf']:
        logs.error(' Outformat for maps must be pdf or png.')
        error = True

    return error



def check_config(configfile):

    gensettings, catalogconf, arrTconf, metaDataconf, RestDownconf,\
        synthsconf, gainfconf, psdsconf, orientconf, timingconf, tc, maps =\
        AutoStatsQConfig.load(filename=configfile).Settings

    error_gen = check_general_settings(gensettings)
    error_cat = check_catalog_settings(catalogconf)
    error_arr = check_arrivaltime_settings(arrTconf)
    error_md = check_metadata_settings(metaDataconf)
    error_rest = check_restitution_settings(RestDownconf)
    error_syn = check_synthetics_settings(synthsconf)
    error_gain = check_gain_settings(gainfconf, RestDownconf)
    error_psd = check_psd_settings(psdsconf)
    error_ori = check_orient_settings(orientconf)
    error_tim = check_timing_settings(timingconf)
    error_tele = check_telecheck_settings(tc)
    error_map = check_maps_settings(maps)


    if error_gen or error_cat or error_arr or error_md or error_rest or error_syn\
      or error_gain or error_psd or error_ori or error_tim or error_tele or error_map:
        
        logs.info('\n')
        logs.error('--> --> One more more errors in config file found. Please check error messages above.')
        sys.exit()