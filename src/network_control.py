import os, sys, glob
import numpy as num
import datetime
import logging
# import linecache
import gc
import math
import argparse

# import matplotlib.pyplot as plt
from pyrocko import util, model, orthodrome, pile, trace, io
from pyrocko import cake, gf
from pyrocko.client import catalog, fdsn
from pyrocko.client.fdsn import EmptyResult
from pyrocko.io import stationxml
from pyrocko import orthodrome as od
# from pyrocko.fdsn import station as fs

# from .gainfactors import *
from . import gainfactors as gainf 
from .catalog import subset_events_dist_cat, subset_events_dist_evlist 
from .catalogplots import *
from .gainplots import plot_median_gain_map_from_file
from .tele_check import TeleCheck
from . import freq_psd as fp
from . import orient
from . import timing as tt
from . import call_tele_check as tele
from.configchecks import check_config
from .config_settings_defaults import generate_default_config
from .config import GeneralSettings, CatalogConfig, ArrTConfig,\
MetaDataDownloadConfig, RestDownRotConfig, SynthDataConfig,\
GainfactorsConfig, PSDConfig, OrientConfig, TimingConfig, TeleCheckConfig,\
maps, AutoStatsQConfig
from .calc_ttt import *
from .make_report import gen_report


fdsn.g_timeout = 120.

'''
Quality control of array stations

0.  Read station lists
    - needs comma-spread list of stations with information on
      net, stat, lat, lom, elev, depth

1.  Catalog search for teleseismic events
    - queries the geofon catalog
    - opt. plots of catalog statistics + map

2.  Subset of events for quality control
    - search for subset of events that will be used for quality control
    - opt. plot: map

3.  Download data and metadata
    - query sites defined in settings, fdsn-download
    - both optional, local data can be used.

4.  Data preparation: restitution of data

5.  Rotation NE --> RT

6.  Synthetic data

7.  Calc. Gain factors (relative)
    - returns list of mean gain factor of each station rel. to ref. station +
      list of gain factors of all stations/events rel. to ref. station
    - opt. plot: x-y plot of all gain factors

8. PSDs
    - comp. synth and obs. PSDs (>30 min twd)
    - providing ranges of flat ratio for MT inv

9. Rayleigh wave polarization analysis for orientation
'''


def get_pl_opt(lats_all, lons_all):
    ''' compute automated map dimensions'''

    lat_m = num.mean(lats_all)
    lon_m = num.mean(lons_all)

    corners = [od.Loc(num.min(lats_all), num.min(lons_all)),
               od.Loc(num.max(lats_all), num.max(lons_all))]

    dist1 = od.distance_accurate50m(od.Loc(lat_m, lon_m), corners[0])
    dist2 = od.distance_accurate50m(od.Loc(lat_m, lon_m), corners[1])
    radius = max(dist1, dist2)*1.2

    print(lat_m, lon_m, radius)

    return [lat_m, lon_m, radius, 'split']


def main(): 
  
    # Any event with "bad" data to be excluded?
    exclude_event = ['2017-11-04 09:00:19.000']

    # Command line input handling
    desc = 'AutoStatsQ - Automated station quality control for MT inversion'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--config')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--generate_config', action='store_true')
    parser.add_argument('--report', action='store_true')
    parser.add_argument('-l', '--loglevel',
                        help='Verbosity in the output.', default='INFO',
                        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO',
                                 'DEBUG'])
    parser.add_argument('--logoutput', '-o', default=None,
                        help='File to save the log')
    args = parser.parse_args()

    # Set verbose level for logging
    verbo = getattr(logging, args.loglevel)
    if args.logoutput is None:
        logging.basicConfig(level=verbo)
    else:
        logging.basicConfig(filename=args.logoutput, level=verbo)

    # Generate a (template) config file:
    if args.generate_config:
        logging.info('Welcome to AutoStatsQ - a station quality control checking tool.\n')
        # Set Logger name and verbosity
        logs = logging.getLogger('Generate config')
        logs.setLevel(verbo)

        fn_config = 'AutoStatsQ_settings.config'
        if os.path.exists('AutoStatsQ_settings.config'):
            logs.error('file exists: %s' % fn_config)

        config = generate_default_config()

        config.dump(filename=fn_config)
        logs.info('created a fresh config file %s' % fn_config)

    if not args.generate_config and not args.config:
        logging.info('Welcome to AutoStatsQ - a station quality control checking tool.\n')
        logging.error('AutoStatsQ needs a config file.')
        print(parser.print_help())

    if args.config:
        logging.info('Checking configuration file.')
        check_config(args.config)

    if args.report:
        logging.info('Generating html report from results.')
        if not args.config:
            logging.error('AutoStatsQ needs a config file.')
        gen_report(args.config)


    # run AutoStatsQ
    if args.run:
        # Set Logger name and verbosity
        logs = logging.getLogger('Run')
        logs.setLevel(verbo)

        logging.info('Welcome to AutoStatsQ - a station quality control checking tool.\n')

        # read existing config file:

        gensettings, catalogconf, arrTconf, metaDataconf, RestDownconf,\
        synthsconf, gainfconf, psdsconf, orientconf, timingconf, tc, maps =\
        AutoStatsQConfig.load(filename=args.config).Settings

        data_dir = gensettings.work_dir
        # os.path.join takes properly into account the case of trailing slash
        os.makedirs(os.path.join(data_dir, 'results'), exist_ok=True)

        sites = metaDataconf.sites

        ''' 0. Read station lists '''
        st_lats = []
        st_lons = []
        ns = []
        all_stations = []

        for stat_list in gensettings.list_station_lists:
            if stat_list.endswith('.csv'):
                with open(stat_list, 'r') as f:
                    for line in f.readlines():
                        if len(line.strip().split(',')) == 6:
                            n, s, lat, lon, elev, d = line.strip().split(',')
                            # Consider elevation and depth as 0 if empty
                            if not len(elev):
                                elev = '0'
                            if not len(d):
                                d = '0'
                            n_s = '%s.%s' % (n, s)
                            if n_s in gensettings.st_use_list or gensettings.st_use_list == []:
                                all_stations.append(model.Station(network=n, station=s,
                                                                  lat=float(lat), lon=float(lon),
                                                                  elevation=float(elev), depth=d))
                                st_lats.append(float(lat))
                                st_lons.append(float(lon))
                                ns.append((n, s))

                        elif len(line.strip().split(',')) == 5:
                            n, s, lat, lon, elev = line.strip().split(',')
                            # Consider elevation and depth as 0 if empty
                            if not len(elev):
                                elev = '0'
                            n_s = '%s.%s' % (n, s)
                            if n_s in gensettings.st_use_list or gensettings.st_use_list == []:
                                all_stations.append(model.Station(network=n, station=s,
                                                                  lat=float(lat), lon=float(lon),
                                                                  elevation=float(elev)))                            
                                st_lats.append(float(lat))
                                st_lons.append(float(lon))
                                ns.append((n, s))

            elif stat_list.endswith('.xml'):
                zs = stationxml.load_xml(filename=stat_list)
                for net in zs.network_list:
                    for stat in net.station_list:
                        n_s = '%s.%s' % (net.code, stat.code)
                        if n_s in gensettings.st_use_list or gensettings.st_use_list == []:

                            st_lats.append(float(stat.latitude.value))
                            st_lons.append(float(stat.longitude.value))
                            ns.append((net.code, stat.code))
                            all_stations.append(model.Station(network=net.code,
                                                station=stat.code,
                                                lat=float(stat.latitude.value),
                                                lon=float(stat.longitude.value),
                                                elevation=float(stat.elevation.value)))

            elif stat_list.endswith('.yaml') or stat_list.endswith('.pf'):
                zs = model.station.load_stations(filename=stat_list)
                for stat in zs:
                    n_s = '%s.%s' % (stat.network, stat.station)
                    if n_s in gensettings.st_use_list or gensettings.st_use_list == []:
                        st_lats.append(float(stat.lat))
                        st_lons.append(float(stat.lon))
                        ns.append((stat.network, stat.station))
                        all_stations.append(model.Station(network=stat.network,
                                            station=stat.station,
                                            lat=float(stat.lat),
                                            lon=float(stat.lon),
                                            elevation=float(stat.elevation)))
            else:
                msg = 'Station file extension not known: %s. Please use .xml, .csv or .yaml.' % stat_list
                logs.error(msg)
                raise Exception('Unknown file extension')


        logs.info(' Number of stations: %d' % len(ns))
        if len(ns) == 0:
            logs.error('No stations found.')
            sys.exit()

        ##### FOR SHORT TESTING
        '''
        all_stations = all_stations[161:175]
        ns = ns[161:175]
        st_lats = st_lats[161:175]
        st_lons = st_lons[161:175]

        print(ns)
        print('stations:', len(ns))
        '''
        #####

        ''' 1. Catalog search for teleseismic events '''
        # Set Logger name and verbosity
        logs = logging.getLogger('Catalog search')
        logs.setLevel(verbo)

        tmin = util.ctimegm(catalogconf.tmin_str)
        tmax = util.ctimegm(catalogconf.tmax_str)
        # os.path.join takes properly into account the case of trailing slash
        os.makedirs(os.path.join(data_dir, 'results', 'catalog'), exist_ok=True)
        ev_catalog = []

        if catalogconf.search_events is True:

            geofon = catalog.GlobalCMT()
            event_names = geofon.get_event_names(
                time_range=(tmin, tmax),
                magmin=catalogconf.min_mag)

            for ev_name in event_names:
                ev_catalog.append(geofon.get_event(ev_name))

            # logs.info('%d events found.' % (len(ev_catalog)))
            catfilename = os.path.join(data_dir, 'results/catalog', 'catalog_Mgr%s.txt' % (catalogconf.min_mag))
            model.dump_events(ev_catalog, catfilename)
            logs.info('length catalog: %d' % len(ev_catalog))

        if catalogconf.use_local_catalog is True:

            if not hasattr(catalogconf, 'catalog_fn') or catalogconf.catalog_fn is None:
                catalogconf.catalog_fn = os.path.join(data_dir, 'results/catalog',
                                                      'catalog_Mgr%s.txt' % (catalogconf.min_mag))

            if catalogconf.subset_of_local_catalog is False:
                ev_catalog = model.load_events(catalogconf.catalog_fn)

            if catalogconf.subset_of_local_catalog is True:
                if not catalogconf.mid_point:

                    mid_point = orthodrome.geographic_midpoint(num.asarray(st_lats),
                                                               num.asarray(st_lons))

                else: 
                    mid_point = catalogconf.mid_point

                ev_catalog = subset_events_dist_cat(catalogconf.catalog_fn,
                                             catalogconf.min_mag,
                                             catalogconf.max_mag,
                                             catalogconf.tmin_str,
                                             catalogconf.tmax_str,
                                             mid_point[0],
                                             mid_point[1],
                                             catalogconf.min_dist_km,
                                             catalogconf.max_dist_km)

            logs.info('length catalog: %d' % len(ev_catalog))

        if not ev_catalog and not catalogconf.use_local_subsets:
            logs.error('A catalog is needed to continue.')
            raise Exception('No catalog!')

        ''' 2. Subset of events for quality control'''
        # Set Logger name and verbosity
        logs = logging.getLogger('Subset events')
        logs.setLevel(verbo)

        subsets_events = {}
        no_bins = int(360/catalogconf.wedges_width)

        if not catalogconf.mid_point:
            mid_point = orthodrome.geographic_midpoint(num.asarray(st_lats),
                                                       num.asarray(st_lons))    
        else: 
            mid_point = catalogconf.mid_point

        for d, val in catalogconf.depth_options.items():
            if not catalogconf.use_local_subsets:
                ev_cat = subset_events_dist_evlist(ev_catalog,
                                                   catalogconf.min_mag,
                                                   catalogconf.max_mag,
                                                   catalogconf.tmin_str,
                                                   catalogconf.tmax_str,
                                                   mid_point[0],
                                                   mid_point[1],
                                                   val[0],
                                                   val[1],
                                                   catalogconf.min_dist_km,
                                                   catalogconf.max_dist_km)


                dist_array = num.empty((len(ev_cat), len(ns)))
                bazi_array = num.empty((len(ev_cat), len(ns)))
                bazi_mp_array = num.empty((len(ev_cat)))

                for i_ev, ev in enumerate(ev_cat):
                    dist_array[i_ev, :] = [float(orthodrome.distance_accurate50m_numpy(
                                          ev.lat, ev.lon, lat, lon))
                                          for (lat, lon) in zip(st_lats, st_lons)]

                    bazi_array[i_ev, :] = [orthodrome.azibazi(ev.lat, ev.lon, lat, lon)[1]
                                           for (lat, lon) in zip(st_lats, st_lons)]

                    bazi_mp_array[i_ev] = orthodrome.azibazi(ev.lat, ev.lon,
                                                             mid_point[0],
                                                             mid_point[1])[1]

                if catalogconf.plot_catalog_all is True:
                    auxdir = os.path.join(data_dir, 'results/catalog')
                    os.makedirs(auxdir, exist_ok=True)
                    pltfilename = 'catalog_global_Mgr%s_%s-%s_%s.%s' % \
                                  (catalogconf.min_mag,
                                   catalogconf.tmin_str[0:10],
                                   catalogconf.tmax_str[0:10], d,
                                   maps.outformat)
                    fn = os.path.join(auxdir, pltfilename)

                    logs.info(' Plotting catalog azimuthal for full catalog: %s.' % d)
                    gmtplot_catalog_azimuthal(ev_cat, mid_point,
                                              catalogconf.dist, fn,
                                              catalogconf.wedges_width)
                
                wedges_array = num.floor(bazi_array / catalogconf.wedges_width)

                # hist based on chosen mid_point of array
                wedges_array_mp = num.floor(bazi_mp_array / catalogconf.wedges_width)
                mean_wedges_mp = num.where(wedges_array_mp < 0,
                                           no_bins+wedges_array_mp,
                                           wedges_array_mp)
                bins_hist = [a for a in range(no_bins+1)]
                hist, bin_edges = num.histogram(mean_wedges_mp, bins=bins_hist)


                if catalogconf.plot_hist_wedges is True:
                    logs.info(' Plotting catalog histogram for full catalog: %s.' % d)
                    plot_catalog_hist(ev_cat, dist_array, mean_wedges_mp,
                                      bins_hist, data_dir, catalogconf.min_mag, d,
                                      catalogconf.wedges_width, no_bins)

                # if catalogconf.plot_dist_vs_magn is True:
                #     logs.info(' Plotting catalog distance vs. magnitude for full catalog: %s.' % d)
                #     plot_distmagn(dist_array, ev_cat, data_dir, d)


                # find 'best' subset of catalog events
                subset_catalog = []

                for bin_nr in range(no_bins):
                    # get indices of all events in current bin:
                    bin_ev_ind = num.argwhere(mean_wedges_mp == bin_nr)

                    if len(bin_ev_ind) == 0:
                        logs.warning('no event for %d - %d deg' % \
                                     (bin_nr*catalogconf.wedges_width,
                                     (bin_nr+1)*catalogconf.wedges_width))

                    if len(bin_ev_ind) == 1:
                        subset_catalog.append(ev_cat[int(bin_ev_ind[0])])

                    if len(bin_ev_ind) > 1:
                        # choose event
                        # if around it bins with no event choose more to that side,
                        # if on both sides
                        # no events if possible two events at both bin margins

                        if bin_nr != 0 and bin_nr != no_bins-1 and\
                          hist[bin_nr-1] == 0 and hist[bin_nr+1] == 0:
                            # choose min und max bazi for better azimuthal coverage
                            min_bazi_ev_ind = bin_ev_ind[num.argmin(bazi_mp_array[bin_ev_ind])][0]
                            max_bazi_ev_ind = bin_ev_ind[num.argmax(bazi_mp_array[bin_ev_ind])][0]
                            subset_catalog.append(ev_cat[min_bazi_ev_ind])
                            subset_catalog.append(ev_cat[max_bazi_ev_ind])

                        elif bin_nr != 0 and hist[bin_nr-1] == 0:
                            # choose one which is more to that side
                            ev_ind_next = bin_ev_ind[
                                                     num.argsort(
                                                                bazi_mp_array[bin_ev_ind],
                                                                axis=0)[0]][0][0]
                            subset_catalog.append(ev_cat[ev_ind_next])

                        elif bin_nr != no_bins-1 and hist[bin_nr+1] == 0:
                            # choose one more to that side
                            ev_ind_next = bin_ev_ind[num.argsort(
                                                                bazi_mp_array[bin_ev_ind],
                                                                axis=0)
                                                     [len(bin_ev_ind)-1]][0][0]
                            subset_catalog.append(ev_cat[ev_ind_next])

                        else:
                            if catalogconf.median_ev_in_bin is True:
                                # only in middle if uniformly distributed eqs in bin!
                                median_ev_ind = bin_ev_ind[num.argsort(
                                                                  bazi_mp_array[bin_ev_ind],
                                                                  axis=0)
                                                           [len(bin_ev_ind)//2]][0][0]
                                subset_catalog.append(ev_cat[median_ev_ind])

                            # if catalogconf.weighted_magn_baz_ev is True:
                            #     # alternatively weighting with magnitude:
                            #     # weighting of magnitude
                            #     mags = [ev.magnitude for ev in [ev_catalog[i_ev[0]]
                            #             for i_ev in bin_ev_ind]]
                            #     max_mag = max(mags)
                            #     min_mag = min(mags)
                            #     W_mag = (num.max(bazi_mp_array[bin_ev_ind]) -
                            #              num.min(bazi_mp_array[bin_ev_ind])) /\
                            #             (max_mag - min_mag) + 0.2
                            #     d_to_opt = []
                            #     baz_op = (bin_nr+1) * catalogconf.wedges_width - catalogconf.wedges_width/2.
                            #     if baz_op > 180:
                            #         baz_op = baz_op - 360
                            #     for ii in bin_ev_ind:
                            #         d_to_opt.append(num.sqrt(
                            #                         W_mag *
                            #                         num.square(max(mags) -
                            #                         ev_catalog[ii[0]].magnitude) +
                            #                         num.square(abs(baz_op -
                            #                         bazi_mp_array[ii]))))

                            #     i_best_ev = bin_ev_ind[num.argmin(d_to_opt)]
                            #     subset_catalog.append(ev_catalog[i_best_ev[0]])
                
                            '''
                            print('best and median event index',
                                  i_best_ev, median_ev_ind,
                                  baz_op, bazi_mp_array[i_best_ev],
                                  bazi_mp_array[median_ev_ind],
                                  catalog[i_best_ev].magnitude,
                                  catalog[median_ev_ind].magnitude)
                            '''
              
                logs.info('Subset of %d events was generated for %s.' % \
                          (len(subset_catalog), d))

                # sort subset catalog by time:
                subset_catalog.sort(key=lambda x: x.time)

                # append to subsets-dict
                subsets_events[d] = subset_catalog

                catfilename = 'catalog_Mgr%s_%s.txt' % (catalogconf.min_mag, d)
                catfullpath = os.path.join(data_dir, 'results/catalog', catfilename)
                model.dump_events(subset_catalog, catfullpath)
                # print([(util.time_to_str(ev.time), ev.magnitude, ev.depth)

            else:
                try:
                    if not hasattr(catalogconf, 'subset_fns') or (d not in catalogconf.subset_fns):
                        catalogconf.subset_fns[d] = os.path.join(data_dir,
                                                                 'results/catalog',
                                                                 'catalog_Mgr%s_%s.txt' % (catalogconf.min_mag, d))

                    subset_catalog = model.load_events(catalogconf.subset_fns[d])
                except Exception:
                    subset_catalog = []

                if subset_catalog == []:
                    logs.error('Catalog empty, %s' % d)
                    sys.exit()

                logs.info(' Subset of catalog for %s: %s events.\n' % (d, len(subset_catalog)))

                subsets_events[d] = subset_catalog

            if catalogconf.plot_hist_wedges is True:
                dist_array_subset = num.empty((len(subset_catalog), len(ns)))
                bazi_array_subset = num.empty((len(subset_catalog), len(ns)))
                bazi_mp_array_subset = num.empty((len(subset_catalog)))

                for i_ev, ev in enumerate(subset_catalog):
                    dist_array_subset[i_ev, :] = [float(orthodrome.distance_accurate50m_numpy(
                                          ev.lat, ev.lon, lat, lon))
                                          for (lat, lon) in zip(st_lats, st_lons)]

                    bazi_array_subset[i_ev, :] = [orthodrome.azibazi(ev.lat, ev.lon, lat, lon)[1]
                                           for (lat, lon) in zip(st_lats, st_lons)]

                    bazi_mp_array_subset[i_ev] = orthodrome.azibazi(ev.lat, ev.lon,
                                                             mid_point[0],
                                                             mid_point[1])[1]

                wedges_array_subset = num.floor(bazi_array_subset / catalogconf.wedges_width)

                # hist based on chosen mid_point of array
                wedges_array_mp_subset = num.floor(bazi_mp_array_subset / catalogconf.wedges_width)
                mean_wedges_mp_subset = num.where(wedges_array_mp_subset < 0,
                                           no_bins+wedges_array_mp_subset,
                                           wedges_array_mp_subset)
                bins_hist_subset = [a for a in range(no_bins+1)]
                hist_subset, bin_edges_subset = num.histogram(mean_wedges_mp_subset, bins=bins_hist_subset)

                logs.info(' Plotting catalog histogram for subset: %s.' % d)
                plot_catalog_hist(subset_catalog, dist_array_subset, mean_wedges_mp_subset,
                                  bins_hist_subset, data_dir, catalogconf.min_mag, d,
                                  catalogconf.wedges_width, no_bins, full_cat=False)

            #if catalogconf.plot_dist_vs_magn is True:
            #    dist_array_subset = num.empty((len(subset_catalog), len(ns)))

            #    for i_ev, ev in enumerate(subset_catalog):
            #        dist_array_subset[i_ev, :] = [float(orthodrome.distance_accurate50m_numpy(
            #                              ev.lat, ev.lon, lat, lon))
            #                              for (lat, lon) in zip(st_lats, st_lons)]

            #    logs.info(' Plotting catalog distance vs. magnitude for subset: %s.' % d)
            #    plot_distmagn(dist_array_subset, subset_catalog, data_dir, d, full_cat=False)

            if catalogconf.plot_catalog_subset is True:
                '''
                plot all events of catalog
                '''
                os.makedirs(os.path.join(data_dir, 'results/catalog'), exist_ok=True)
                _tmin = util.time_to_str(min([ev.time for ev in subset_catalog]))
                _tmax = util.time_to_str(max([ev.time for ev in subset_catalog]))                
                catfilename = 'catalog_global_Mgr%s_%s-%s_%s_subset.%s' % \
                              (catalogconf.min_mag, _tmin[0:10], _tmax[0:10], d, maps.outformat)
                fn = os.path.join(data_dir, 'results/catalog', catfilename)

                logs.info(' Plotting catalog azimuthal for subset: %s.' % d)
                gmtplot_catalog_azimuthal(subset_catalog, mid_point, 
                                               catalogconf.dist, fn, catalogconf.wedges_width)

            if exclude_event != []:
                new_subset_catalog = []
                for ev in subset_catalog:
                    if util.time_to_str(ev.time) not in exclude_event:
                        new_subset_catalog.append(ev)
                subset_catalog = new_subset_catalog


            ''' 2.1 Calculate arrival times for all event/station pairs '''
            # Method a) using cake
            arrT_array = None
            arrT_R_array = None
            if arrTconf.calc_first_arr_t is True:
                logs.info(' Calculating arrival times of first phase, %s subset.'
                          % d)
                data_dir = gensettings.work_dir
                os.makedirs(os.path.join(data_dir, 'ttt'), exist_ok=True)
                dist_array_sub = num.empty((len(subset_catalog), len(ns)))

                for i_ev, ev in enumerate(subset_catalog):
                    dist_array_sub[i_ev, :] = [float(orthodrome.distance_accurate50m_numpy(
                                          ev.lat, ev.lon, lat, lon))
                                          for (lat, lon) in zip(st_lats, st_lons)]

                arrT_array = num.empty((dist_array_sub.shape))
                depths = [ev.depth for ev in subset_catalog]
                vmodel = cake.load_model('prem-no-ocean.f')
                phases = [cake.PhaseDef(pid) for pid in arrTconf.phase_select.split('|')]
                nev = len(subset_catalog)
                for i_ev, ev in enumerate(subset_catalog):
                    ds = depths[i_ev]
                    print('Event: %5d/%s' % (i_ev,nev), end='\r')
                    for i_st in range(len(ns)):
                        dist = dist_array_sub[i_ev, i_st]
                        arrivals = vmodel.arrivals(distances=[dist*cake.m2d],
                                                   phases=phases,
                                                   zstart=ds)

                        min_t = min(arrivals, key=lambda x: x.t).t

                        arrT_array[i_ev, i_st] = ev.time + min_t

                num.save(os.path.join(data_dir, 'ttt', 'ArrivalTimes_%s' % (d)), arrT_array)

                logs.info(' Arrival times of first phase ready, %s subset.\n'
                          % d)

            if arrTconf.calc_est_R is True:
                logs.info(' Calculating arrival times of Rayleigh waves, %s subset.'
                          % d)
                data_dir = gensettings.work_dir
                os.makedirs(os.path.join(data_dir, 'ttt'), exist_ok=True)
                dist_array_sub = num.empty((len(subset_catalog), len(ns)))
                nev = len(subset_catalog)

                for i_ev, ev in enumerate(subset_catalog):
                    print('Event: %5d/%s' % (i_ev,nev), end='\r')
                    dist_array_sub[i_ev, :] = [float(orthodrome.distance_accurate50m_numpy(
                                          ev.lat, ev.lon, lat, lon))
                                          for (lat, lon) in zip(st_lats, st_lons)]

                arrT_R_array = num.empty((dist_array_sub.shape))

                for i_ev, ev in enumerate(subset_catalog):
                    for i_st in range(len(ns)):
                        dist = dist_array_sub[i_ev, i_st]  # m
                        #print('------')
                        #print('d', dist)
                        #print('dt', dist/4000.)
                        #print('t', ev.time)
                        #print('new t', ev.time + dist/4000.)
                        arrT_R_array[i_ev, i_st] = ev.time + dist/(arrTconf.v_rayleigh*1000.)

                num.save(os.path.join(data_dir, 'ttt', 'ArrivalTimes_estR_%s' % (d)), arrT_R_array)

                logs.info(' Arrival times of Rayleigh waves ready, %s subset.\n'
                          % d)

            # Method b) interpolating from fomosto travel time tables
            '''
            arrT_array = None
            if ct.calc_ttt is True:
                data_dir = gensettings.work_dir
                os.makedirs(data_dir+'ttt', exist_ok=True)
                arrT_array_ttt = num.empty((len(subset_catalog), len(ns)))

                n = len(subset_catalog)
                # array n*3 (r_depth, s_depth, distance) in m
                coords = num.zeros((n,3))
                # coords[:,0] receiver_depths are zero
                coords[:,1] = [ev.depth for ev in subset_catalog]
                print('depth_min', num.min(coords[:,1]))
                print('depth_max', num.max(coords[:,1]))

                origin_times = num.asarray([ev.time for ev in subset_catalog])

                for i_st, (lat, lon) in enumerate(zip(st_lats, st_lons)):
                    coords[:,2] = [float(orthodrome.distance_accurate50m_numpy(
                                          ev.lat, ev.lon, lat, lon))
                                          for ev in subset_catalog]

                    intp_times = num.asarray(get_ttt(ct, coords, val))

                    arrT_array_ttt[:, i_st] = origin_times + intp_times
                    if i_st ==2:
                        print(coords[:])
                        print(arrT_array_ttt[:, i_st])

                num.save('%sttt/ArrivalTimes_%s' % (data_dir, d), arrT_array_ttt)
                # print(arrT_array_ttt[0:3,0:50])
            '''    

        ''' 3. Download data and metadata '''
        # Set Logger name and verbosity
        logs = logging.getLogger('Download data')
        logs.setLevel(verbo)

        data_pile = None
        # token = open(metaDataconf.token, 'rb').read()

        if metaDataconf.download_data is True:   ### clean up!
            
            logs.info(' Starting data downloading section.')

            for subset_catalog in subsets_events.values():
                for ev in subset_catalog:
                    ev_t_str = util.time_to_str(ev.time)
                    ev_t = datetime.datetime.strptime(ev_t_str, "%Y-%m-%d %H:%M:%S.%f")
                    dt_s = datetime.timedelta(hours=metaDataconf.dt_start)
                    dt_e = datetime.timedelta(hours=metaDataconf.dt_end)
                    t_start = util.str_to_time(str(ev_t - dt_s),
                                               format='%Y-%m-%d %H:%M:%S.OPTFRAC')
                    t_end = util.str_to_time(str(ev_t + dt_e),
                                             format='%Y-%m-%d %H:%M:%S.OPTFRAC')
                    ev_t_str = ev_t_str.replace(' ', '_')
                    ev_dir = ev_t_str + '/'
                    dir_make = os.path.join(data_dir, ev_dir)
                    os.makedirs(dir_make, exist_ok=True)

                    for ns_now in ns:
                        mseed_fn_st = os.path.join(dir_make, '%s_%s_%s' %
                                                   (ns_now[0], ns_now[1], ev_t_str))

                        if glob.glob(mseed_fn_st + '*'):
                            continue
                        #else:
                        #    print(glob.glob(mseed_fn_st+'*'))
                        
                        selection = [(ns_now[0], ns_now[1], '*',
                                      metaDataconf.channels_download,
                                      t_start, t_end)]
                        logs.debug(selection)

                        for site in sites:
                            mseed_fn = mseed_fn_st + site + 'tr.mseed'

                            if site in metaDataconf.token:
                                token = open(metaDataconf.token[site], 'rb').read()
                            else:
                                token = None
                            try:
                                if token:
                                    request_waveform = fdsn.dataselect(site=site,
                                                                       selection=selection,
                                                                       token=token)
                                else:
                                    request_waveform = fdsn.dataselect(site=site,
                                                                       selection=selection)

                                with open(mseed_fn, 'wb') as wffile:
                                    wffile.write(request_waveform.read())

                            except fdsn.EmptyResult:
                                logs.warning('%s no data from %s' % (ns_now, site))

                            except:
                                logs.error('exception unknown %s' % (ns_now,))

                            else:
                                logs.debug('%s data downloaded from %s' % (ns_now, site))
                                break

                    if not os.listdir(dir_make):
                        os.rmdir(dir_make)

            logs.info(' Finished data downloading section.\n')


        if metaDataconf.download_metadata is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Download metadata')
            logs.setLevel(verbo)

            logs.info(' Starting to download station metadata.')

            cat_tmin = min([ev.time for ev in subset_catalog for k, subset_catalog in subsets_events.items()])
            cat_tmax = max([ev.time for ev in subset_catalog for k, subset_catalog in subsets_events.items()])

            # one file for all events
            selection = []
            for ns_now in ns:
                selection.append((ns_now[0], ns_now[1], '*',
                                  metaDataconf.channels_download,
                                  cat_tmin, cat_tmax))

            meta_fn = os.path.join(data_dir, 'Resp_all')

            for site in sites:
                # This sometimes does not work properly, why? Further testing?...
                logs.info(site)
                try:
                    request_response = fdsn.station(
                            site=site, selection=selection, level='response')
                except EmptyResult:
                    logs.warning('no metadata for %s from %s' % (selection[1], site))
                    continue
                request_response.dump_xml(filename='%s_%s.xml' % (meta_fn, site))
                #except:
                #    print('no metadata at all', site, selection[1])

            logs.info(' Finished to download station metadata.\n')


        ''' 4. Data preparation: restitution of data '''
        if RestDownconf.rest_data is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Restitution')
            logs.setLevel(verbo)

            logs.info(' Starting restitution of data.')
            responses = []

            if metaDataconf.local_metadata != []:
                logs.debug(' Looking for local metadata.')
                for file in metaDataconf.local_metadata:
                    responses.append(stationxml.load_xml(filename=file))
            
            if metaDataconf.use_downmeta is True:
                logs.debug(' Looking for downloaded metadata.')
                for site in sites:
                    stations_fn = os.path.join(data_dir, 'Resp_all_' + str(site) + '.xml')
                    responses.append(stationxml.load_xml(filename=stations_fn))

            if not metaDataconf.use_downmeta and not metaDataconf.local_metadata:
                logs.error(' No response files found. Set *use_downmeta* to True '
                           + 'or provide other local metadata using option'
                           + ' *local_metadata*.')
                sys.exit()

            i_resp = len(responses)
            logs.info(' Number of responses found: %s.' % i_resp)

            if metaDataconf.local_data and not metaDataconf.sds_structure:
                logs.info('Accessing local data.')
                p_local = pile.make_pile(paths=metaDataconf.local_data,
                                         show_progress=True)

            nst = len(all_stations)
            for key, subset_catalog in subsets_events.items(): 
                logs.info(' Starting restitution of data, %s subset.' % d)
                nev = len(subset_catalog)

                for i_ev, ev in enumerate(subset_catalog):
                    logs.debug(util.time_to_str(ev.time))

                    ev_t_str = util.time_to_str(ev.time).replace(' ', '_')
                    tmin = ev.time+metaDataconf.dt_start*3600
                    tmax = ev.time+metaDataconf.dt_end*3600

                    # stations_fn = data_dir + ev_t_str + '_resp_geofon.xml'
                    # response = stationxml.load_xml(filename=stations_fn)
                    # print(data_dir+ev_t_str)
                    if metaDataconf.local_waveforms_only is False:
                        p = pile.make_pile(paths=os.path.join(data_dir, ev_t_str), 
                                           show_progress=False)
                    
                    dir_make = os.path.join(data_dir, 'rest', ev_t_str)
                    os.makedirs(dir_make, exist_ok=True)
                    transf_taper = 1/min(RestDownconf.freqlim)

                    for i_st, st in enumerate(all_stations):
                        print('Event: %5d/%s, Station: %5d/%s' % (i_ev, nev, i_st, nst),
                              end='\r')
                        #if st.network not in net_check:
                        #    continue                        
                        nsl = st.nsl()
                        trs = []
                        if not metaDataconf.local_waveforms_only:
                            #trs = p.all(
                            #    trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2])
                            trs.extend(p.all(tmin=tmin, tmax=tmax,
                                             trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2]))

                        if metaDataconf.local_data:
                            if metaDataconf.sds_structure is False:
                                trs.extend(p_local.all(tmin=tmin, tmax=tmax,
                                                 trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2]))

                            else:
                                year = util.time_to_str(ev.time)[0:4]
                                jul_day = util.julian_day_of_year(ev.time)
                                local_data_dirs = metaDataconf.local_data
                                for i_ldd, ldd in enumerate(local_data_dirs):
                                    path_str = os.path.join(ldd, year, st.network, st.station)
                                    p = pile.make_pile(paths=path_str, regex='.%s' % jul_day, show_progress=False)
                                    trs.extend(p.all(tmin=tmin, tmax=tmax,
                                                     trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2]))
                            # trace.snuffle(trs)

                        if trs:
                            comps = [tr.channel for tr in trs]

                            if metaDataconf.all_channels == False:
                                
                                if 'HHZ' in comps and 'HHN' in comps and 'HHE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HHZ', 'HHN', 'HHE']]
                                elif 'HHZ' in comps and 'HH2' in comps and 'HH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HHZ', 'HH2', 'HH3']]
                                elif 'HHZ' in comps and 'HH1' in comps and 'HH1' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HHZ', 'HH1', 'HH2']]

                                elif 'BHZ' in comps and 'BHN' in comps and 'BHE' in comps:
                                    trs = [tr for tr in trs  if tr.channel in ['BHZ', 'BHN', 'BHE']] 
                                elif 'BHZ' in comps and 'BH2' in comps and 'BH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['BHZ', 'BH2', 'BH3']]
                                elif 'BHZ' in comps and 'BH1' in comps and 'BH2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['BHZ', 'BH1', 'BH2']] 

                                elif 'HNZ' in comps and 'HNN' in comps and 'HNE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HNZ', 'HNN', 'HNE']]
                                elif 'HNZ' in comps and 'HN2' in comps and 'HN3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HNZ', 'HN2', 'HN3']]
                                elif 'HNZ' in comps and 'HN1' in comps and 'HN2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['HNZ', 'HN1', 'HN2']]

                                elif 'EHZ' in comps and 'EHN' in comps and 'EHE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['EHZ', 'EHN', 'EHE']]
                                elif 'EHZ' in comps and 'EH2' in comps and 'EH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['EHZ', 'EH2', 'EH3']]                            
                                elif 'EHZ' in comps and 'EH1' in comps and 'EH2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['EHZ', 'EH1', 'EH2']] 

                                elif 'LHZ' in comps and 'LHN' in comps and 'LHE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['LHZ', 'LHN', 'LHE']]
                                elif 'LHZ' in comps and 'LH2' in comps and 'LH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['LHZ', 'LH2', 'LH3']]                            
                                elif 'LHZ' in comps and 'LH1' in comps and 'LH2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['LHZ', 'LH1', 'LH2']]

                                elif 'CNZ' in comps and 'CNN' in comps and 'CNE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CNZ', 'CNN', 'CNE']]
                                elif 'CNZ' in comps and 'CN2' in comps and 'HH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CNZ', 'CN2', 'CN3']]
                                elif 'CNZ' in comps and 'CN1' in comps and 'CN2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CNZ', 'CN1', 'CN2']]

                                elif 'CHZ' in comps and 'CHN' in comps and 'CHE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CHZ', 'CHN', 'CHE']]
                                elif 'CHZ' in comps and 'CH2' in comps and 'CH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CHZ', 'CH2', 'CH3']]
                                elif 'CHZ' in comps and 'CH1' in comps and 'CH2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['CHZ', 'CH1', 'CH2']]

                                elif 'DNZ' in comps and 'DNN' in comps and 'DNE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DNZ', 'DNN', 'DNE']]
                                elif 'DNZ' in comps and 'DN2' in comps and 'DN3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DNZ', 'DN2', 'DN3']]
                                elif 'DNZ' in comps and 'DN1' in comps and 'DN2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DNZ', 'DN1', 'DN2']]
                                elif 'DN1' in comps and 'DN2' in comps and 'DN3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DN1', 'DN2', 'DN3']]

                                elif 'DHZ' in comps and 'DHN' in comps and 'DHE' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DHZ', 'DHN', 'DHE']]
                                elif 'DHZ' in comps and 'DH2' in comps and 'DH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DHZ', 'DH2', 'DH3']]
                                elif 'DHZ' in comps and 'DH1' in comps and 'DH2' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DHZ', 'DH1', 'DH2']]
                                elif 'DH1' in comps and 'DH2' in comps and 'DH3' in comps:
                                    trs = [tr for tr in trs if tr.channel in ['DH1', 'DH2', 'DH3']]

                                else:
                                    # print('no BH* or HH* data for station %s found' % (str(nsl)))
                                    # print('found these:', [tr.channel for tr in trs])
                                    continue

                                for tr in trs:
                                    cnt_resp = 0
                                    if RestDownconf.deltat_down > 0.0 and not RestDownconf.deltat_down == tr.deltat:
                                        tr.downsample_to(0.1)
                                    for resp_now in responses:
                                        try:
                                            polezero_resp = resp_now.get_pyrocko_response(
                                                nslc=tr.nslc_id,
                                                timespan=(tr.tmin, tr.tmax),
                                                fake_input_units='M'
                                                )

                                            restituted = tr.transfer(
                                                tfade=transf_taper,
                                                freqlimits=RestDownconf.freqlim,
                                                transfer_function=polezero_resp,
                                                invert=True)

                                            fname = '%s_%s_%s__%s_%srest2.mseed' \
                                                    % (tr.nslc_id[0], tr.nslc_id[1],
                                                       tr.nslc_id[2], tr.nslc_id[3],
                                                       ev_t_str)
                                            rest_fn = os.path.join(dir_make, fname)
                                            # rest_fn = dir_make + '/' + str(tr.nslc_id[0]) + '_' + \
                                            #                     str(tr.nslc_id[1])\
                                            #           + '_' + str(tr.nslc_id[2]) + '_' +\
                                            #           '_' + str(tr.nslc_id[3]) + '_' +\
                                            #           ev_t_str + 'rest2.mseed'
                                            io.save(restituted, rest_fn)

                                        except stationxml.NoResponseInformation:
                                            cnt_resp += 1
                                            if cnt_resp == i_resp:
                                                logs.warning('no resp found: %s' % (tr.nslc_id,))

                                        except trace.TraceTooShort:
                                            logs.error('trace too short: %s' % (tr.nslc_id,))

                                        except ValueError:
                                            logs.error('downsampling does not work: %s' % (tr.nslc_id,))

                                        else:
                                            break

                            if metaDataconf.all_channels == True:
                                ch_list = set([ch[0:2] for ch in comps])
                                for ch_ in ch_list:
                                    if ch_+'N' in comps and ch_+'E' in comps and ch_+'Z' in comps:
                                        trs = [tr for tr in trs if tr.channel in [ch_+'N', ch_+'E', ch_+'Z']]
                                    else:
                                      continue

                                    for tr in trs:
                                        if RestDownconf.deltat_down > 0.0 and not RestDownconf.deltat_down == tr.deltat:
                                            tr.downsample_to(0.1)
                                        cnt_resp = 0
                                        for resp_now in responses:
                                            try:
                                                polezero_resp = resp_now.get_pyrocko_response(
                                                    nslc=tr.nslc_id,
                                                    timespan=(tr.tmin, tr.tmax),
                                                    fake_input_units='M'
                                                    )

                                                restituted = tr.transfer(
                                                    tfade=transf_taper,
                                                    freqlimits=RestDownconf.freqlim,
                                                    transfer_function=polezero_resp,
                                                    invert=True)

                                                fname = '%s_%s_%s__%s_%srest2.mseed' % \
                                                        (tr.nslc_id[0], tr.nslc_id[1],
                                                         tr.nslc_id[2], tr.nslc_id[3],
                                                         ev_t_str)
                                                rest_fn = os.path.join(dir_make, fname)
                                                # rest_fn = dir_make + '/' + str(tr.nslc_id[0]) + '_' +\
                                                #           str(tr.nslc_id[1])\
                                                #           + '_' + str(tr.nslc_id[2]) + '_' +\
                                                #           '_' + str(tr.nslc_id[3]) + '_' +\
                                                #           ev_t_str + 'rest2.mseed'
                                                io.save(restituted, rest_fn)

                                            except stationxml.NoResponseInformation:
                                                cnt_resp += 1
                                                if cnt_resp == i_resp:
                                                    logs.error('no resp found: %s' % tr.nslc_id)

                                            except trace.TraceTooShort:
                                                logs.error('trace too short: %s' % tr.nslc_id)

                                            except ValueError:
                                                logs.error('downsampling does not work: %s' % tr.nslc_id)

                                            else:
                                                break
                    if not os.listdir(dir_make):
                        os.rmdir(dir_make)

            logs.info(' Finished restitution of data, %s subset.' % d)
            logs.info(' Saved restituted data in directory %s.\n' % os.path.join(data_dir, 'rest'))


        ''' 5. Rotation NE --> RT '''
        if RestDownconf.rotate_data is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Rotation')
            logs.setLevel(verbo)

            def save_rot_down_tr(tr, dir_rot, ev_t_str):
                fname = '%s_%s_%s__%s_%srrd2.mseed' % \
                        (tr.nslc_id[0], tr.nslc_id[1],
                         tr.nslc_id[2], tr.nslc_id[3],
                         ev_t_str)
                rot_fn = os.path.join(dir_rot, fname)
                # rot_fn = dir_rot + '/' + str(tr.nslc_id[0]) + '_' +\
                #          str(tr.nslc_id[1])\
                #          + '_' + str(tr.nslc_id[2]) + '_' +\
                #          '_' + str(tr.nslc_id[3]) + '_' +\
                #          ev_t_str + 'rrd2.mseed'
                io.save(tr, rot_fn)

            def downsample_rotate(dir_rest, dir_rot, all_stations, st_xml, deltat_down):

                ev_data_pile = pile.make_pile(dir_rest, regex='rest2',
                                              show_progress=False)
                for st in all_stations:

                    #if st.network not in net_check:
                    #    continue 
                    nsl = st.nsl()
                    trs = ev_data_pile.all(
                            trace_selector=lambda tr: tr.nslc_id[:2] == nsl[:2],
                            want_incomplete=True)
                    if not trs:
                        continue
                    chan_avail = [tr.channel for tr in trs]
                    loc_avail = [tr.location for tr in trs]
                    # lens_trs = [len(tr.ydata) for tr in trs]

                    ch_list = set([ch[0:2] for ch in chan_avail])
                    
                    for l in loc_avail:
                        for ch_ in ch_list:
                            lens_trs = [len(tr.ydata) for tr in trs if tr.channel[0:2] == ch_ and tr.location == l]
                            trs_ch = [tr for tr in trs if tr.channel[0:2] == ch_ and tr.location == l]
                            if not lens_trs or 0 in lens_trs or len(lens_trs) != 3:
                                continue

                            else:
                                az1 = None
                                tr1 = None
                                tr2 = None

                                nslcs = [tr.nslc_id for tr in trs_ch]
                                tmin = max([tr.tmin for tr in trs_ch])
                                tmax = min([tr.tmax for tr in trs_ch])
                                trmin = math.ceil(tmin)
                                trmax = int(tmax)
                                for st_now in st_xml:
                                    stat = st_now.get_pyrocko_stations(nslcs=nslcs,
                                           timespan=(trmin, trmax))
                                    if len(stat) == 1:
                                        break

                                if len(stat) != 0:
                                    stat = stat[0]
                                    tr1_ch = '0'
                                    tr2_ch = '0'

                                    test = []
                                    naming = ''
                                    for tr in trs_ch:
                                        if tr.channel.endswith('2'):
                                            test.append('2')
                                        elif tr.channel.endswith('3'):
                                            test.append('3')
                                        elif tr.channel.endswith('1'):
                                            test.append('1')
                                    
                                    if '1' in test and '2' in test and '3' not in test:
                                        naming = '1,2'
                                    elif  '2' in test and '3' in test and '1' not in test:
                                        naming = '2,3'
                                    elif '1' in test and '2' in test and '3' in test:
                                        naming = '1,2,3'

                                    for tr in trs_ch:
                                        if tr.channel.endswith('N') or\
                                         (tr.channel.endswith('2') and naming == '2,3')\
                                          or (tr.channel.endswith('1') and naming == '1,2')\
                                          or (tr.channel.endswith('2') and naming == '1,2,3'):

                                            for ch in stat.channels:
                                                if tr.channel == ch.name:
                                                    az1 = ch.azimuth
                                                    tr1 = tr.copy()

                                                    try:
                                                        tr1.chop(trmin, trmax)
                                                        if deltat_down > 0.0:
                                                            tr1.downsample_to(deltat=deltat_down)
                                                        save_rot_down_tr(tr1, dir_rot, ev_t_str)

                                                    except trace.NoData:
                                                        logs.error('N/2 comp no data in twd %s' % nsl)
                                                        tr1 = None
                                                    except ValueError:
                                                        tr1 = None
                                                        logs.error('N/2 downsampling not successful')
                                                    except util.UnavailableDecimation:
                                                        logs.error('unavailable decimation %s' % tr1.station)
                                                        tr1 = None

                                        if tr.channel.endswith('E') or\
                                         (tr.channel.endswith('3') and naming == '2,3')\
                                          or (tr.channel.endswith('2') and naming == '1,2')\
                                          or (tr.channel.endswith('1') and naming == '1,2,3'):

                                            for ch in stat.channels:
                                                if tr.channel == ch.name:
                                                    tr2 = tr.copy()

                                                    try:
                                                        tr2.chop(trmin, trmax)
                                                        if deltat_down > 0.0:                                                
                                                            tr2.downsample_to(deltat=deltat_down)
                                                        save_rot_down_tr(tr2, dir_rot, ev_t_str)

                                                    except trace.NoData:
                                                        logs.error('E/3 comp no data in twd %s' % nsl)
                                                        tr2 = None
                                                    except ValueError:
                                                        tr2 = None
                                                        logs.error('E/3 downsampling not successful')
                                                    except util.UnavailableDecimation:
                                                        logs.error('unavailable decimation %s' % tr2.station)
                                                        tr2 = None

                                        if tr.channel.endswith('Z')\
                                           or tr.channel.endswith('3') and naming == '1,2,3':
                                            trZ = tr.copy()

                                            try:
                                                trZ.chop(trmin, trmax)
                                                if deltat_down > 0.0:
                                                    trZ.downsample_to(deltat=deltat_down)
                                                trZ.set_channel('Z')
                                                save_rot_down_tr(trZ, dir_rot, ev_t_str)

                                            except trace.NoData:
                                                logs.error('E/3 comp no data in twd %s' % nsl)
                                            except ValueError:
                                                trZ = None
                                                logs.error('Z downsampling not successful')
                                            except util.UnavailableDecimation:
                                                logs.error('unavailable decimation %s' % trZ.station)
                                                trZ = None

                                if az1 is not None and tr1 is not None\
                                  and tr2 is not None:

                                    stat.set_event_relative_data(ev)
                                    baz = stat.backazimuth
                                    az_r = baz + 180 - az1
                                   
                                    if str(tr1.channel).endswith('N') is True\
                                       or str(tr1.channel).endswith('2') is True and naming == '2,3'\
                                       or str(tr1.channel).endswith('1') is True and naming == '1,2'\
                                       or str(tr1.channel).endswith('2') is True and naming == '1,2,3':
                                        tr1_ch = tr1.channel

                                    if str(tr1.channel).endswith('E') is True\
                                      or str(tr1.channel).endswith('3') is True and naming == '2,3'\
                                      or str(tr1.channel).endswith('2') is True and naming == '1,2'\
                                      or str(tr1.channel).endswith('1') is True and naming == '1,2,3':
                                        tr2_ch = tr1.channel

                                    if str(tr2.channel).endswith('N') is True\
                                     or str(tr2.channel).endswith('2') is True and naming == '2,3'\
                                     or str(tr2.channel).endswith('1') is True and naming == '1,2'\
                                     or str(tr2.channel).endswith('2') is True and naming == '1,2,3':
                                        tr1_ch = tr2.channel

                                    if str(tr2.channel).endswith('E') is True\
                                     or str(tr2.channel).endswith('3') is True and naming == '2,3'\
                                     or str(tr2.channel).endswith('2') is True and naming == '1,2'\
                                     or str(tr2.channel).endswith('1') is True and naming == '1,2,3':
                                        tr2_ch = tr2.channel

                                    if tr1_ch != '0' and tr2_ch != '0':
                                        tracescopy=[tr1,tr2]
                                        tracestouse=[]
                                        tracestouse=[tr for tr in tracescopy if num.any(num.diff(tr.ydata))]
                                        lens_trsuse = [len(tr.ydata) for tr in tracestouse ]
                                        if not lens_trsuse or 0 in lens_trsuse or len(lens_trsuse) != 2:
                                            continue
                                        if len(lens_trsuse)==2:
                                            rots = trace.rotate(tracestouse, azimuth=az_r, 
                                                            in_channels=[tr1_ch, tr2_ch],
                                                            out_channels=['R', 'T'])
                                            for tr in rots:
                                                save_rot_down_tr(tr, dir_rot, ev_t_str)

            st_xml = []
            if metaDataconf.local_metadata != []:
                for file in metaDataconf.local_metadata:
                    st_xml.append(stationxml.load_xml(filename=file))
            
            if metaDataconf.use_downmeta is True:            
                for site in sites:
                    stations_fn = os.path.join(data_dir, 'Resp_all_' + str(site) 
                                               + '.xml')
                    st_xml.append(stationxml.load_xml(filename=stations_fn))

            i_st_xml = len(st_xml)
            for key, subset_catalog in subsets_events.items():
                logs.info(' Starting downsampling and rotation, subset %s' % key)
                nev = len(subset_catalog)

                for i_ev, ev in enumerate(subset_catalog):
                    logs.debug(' Event %s' % util.time_to_str(ev.time))
                    print('Event: %5d/%s' % (i_ev,nev), end='\r')
                    gc.collect()
                    ev_t_str = util.time_to_str(ev.time).replace(' ', '_')
                    os.makedirs(os.path.join(data_dir, 'rrd'), exist_ok=True)
                    dir_rot = os.path.join(data_dir, 'rrd', ev_t_str)
                    dir_rest = os.path.join(data_dir, 'rest', ev_t_str)
                    downsample_rotate(dir_rest, dir_rot, all_stations, st_xml, 
                                      RestDownconf.deltat_down)
                    logs.debug(' Saved ev %s' % util.time_to_str(ev.time))

                    if not os.listdir(os.path.join(data_dir, 'rrd')):
                        os.rmdir(os.path.join(data_dir, 'rrd'))

            logs.info(' Done with rotation to ZRT. Saved rotated and'
                      + ' downsampled data in directory %s.\n' % 
                      os.path.join(data_dir, 'rrd'))


        ''' 6. Synthetic data '''
        if synthsconf.make_syn_data is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Synthetics')
            logs.setLevel(verbo)

            logs.info(' Starting synthetic data generation section.')

            if not os.path.exists(synthsconf.engine_path):
                logs.error('ERROR: Engine path not found: %s.' % synthsconf.engine_path
                           + ' Please add a valid path to *engine_path* to compute synthetic data.' )
                sys.exit()

            if not os.path.exists(os.path.join(synthsconf.engine_path, synthsconf.store_id)):
                logs.error('ERROR: GF store not found: %s.' % os.path.join(synthsconf.engine_path, synthsconf.store_id))
                sys.exit()

            logs.info(' GF store found: %s.' % os.path.join(synthsconf.engine_path, synthsconf.store_id))

            freqlim = RestDownconf.freqlim
            transf_taper = 1/min(freqlim)
            store_id = synthsconf.store_id
            engine = gf.LocalEngine(store_superdirs=
                                   [synthsconf.engine_path])
            os.makedirs(os.path.join(data_dir, 'synthetics'), exist_ok=True)
            loc = '0'
            nst = len(all_stations)

            for key, subset_catalog in subsets_events.items():
                logs.info(' Generating synthetic data for %s subset.' % key)
                nev = len(subset_catalog)

                for i_ev, ev in enumerate(subset_catalog):
                    
                    ev_t_str = util.time_to_str(ev.time).replace(' ', '_')
                    #logs.info(ev_t_str)

                    dir_syn_ev = os.path.join(data_dir, 'synthetics', ev_t_str)
                    os.makedirs(dir_syn_ev, exist_ok=True)
                    # ev.duration = ev.
                    # ev.duration
                    # oder
                    # source.stf = gf.BoxcarSTF(duration=)
                    # scaling mit magnitude
                    if ev.duration < 1:
                        logs.warning(' Warning ev.duration: %s' % ev.duration)
                    #ev.duration = None
                    ev.time = ev.time + ev.duration/2
                    source = gf.MTSource.from_pyrocko_event(ev)

                    for ist, st in enumerate(all_stations):
                        print('Event: %5d/%s, Station %5d/%s' 
                              % (i_ev,nev, ist,nst), end='\r')
                        targets = []
                        sta = st.station
                        net = st.network

                        for cha in ['Z', 'R', 'T']:
                            target = gf.Target(
                                    codes=(net, sta, loc, cha),
                                    quantity='displacement',
                                    lat=st.lat,
                                    lon=st.lon,
                                    store_id=store_id)
                            azi, _ = target.azibazi_to(ev)

                            if cha == 'R':
                                target.azimuth = azi - 180.
                                target.dip = 0.

                            elif cha == 'T':
                                target.azimuth = azi - 90.
                                target.dip = 0.

                            if cha == 'Z':
                                target.azimuth = 0.
                                target.dip = -90.

                            targets.append(target)

                            try:
                                response = engine.process(source, targets)
                                trs_syn = response.pyrocko_traces()

                            except Exception:
                                continue

                            else:

                                for tr in trs_syn:
                                    # adding zeros at beginning of trace for taper in
                                    # tr.transfer
                                    tr.extend(tmin=tr.tmin-transf_taper)
                                    tr.transfer(
                                                tfade=transf_taper,
                                                freqlimits=freqlim,
                                                invert=False)
                                    net, sta, loc, cha = tr.nslc_id

                                    filename = '%s/syn_%s_%s_%s_%s.mseed'\
                                             % (dir_syn_ev, net, sta, cha, ev_t_str)
                                    io.save(tr, filename, format='mseed')

            logs.info(' Done with computation of synth. waveforms.'
                      + ' Saved in directory %s. \n' % 
                      os.path.join(data_dir, 'synthetics'))


        ''' 7. Gain factors '''
        if gainfconf.calc_gainfactors is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Gain factors')
            logs.setLevel(verbo)

            logs.info(' Starting evaluation of gain factors.')
            logs.info(' Method: %s' % gainfconf.gain_factor_method)

            dir_gains = os.path.join(data_dir, 'results', 'gains')
            os.makedirs(dir_gains, exist_ok=True)
            twd = (gainfconf.wdw_st_arr, gainfconf.wdw_sp_arr)

            # arrival times
            if arrT_array is None:
                try:
                    data_dir = gensettings.work_dir
                    atfile = os.path.join(data_dir, 'ttt', 'ArrivalTimes_deep.npy')
                    arrT_array = num.load(atfile)
                except Exception:
                    logs.error('Please calculate arrival times first!')
                    raise Exception('Arrival times not calculated!')
                    sys.exit()

            def run_autogain(data_dir, all_stations, subset_catalog,
                             gain_factor_method, dir_gains, 
                             twd, arrT_array, comp):
                
                datapath = os.path.join(data_dir, 'rrd')

                if len(gain_factor_method) == 1:
                    gain_factor_method = gain_factor_method[0]
                data_pile = pile.make_pile(datapath, show_progress=False)

                if data_pile.is_empty():
                    logs.error(' No data found in %s. Gain test failed.' % datapath)
                    sys.exit()
                
                syn_data_pile = None
                if gain_factor_method == 'syn_comp':
                    syn_data_pile = pile.make_pile(os.path.join(data_dir, 
                                                   'synthetics'))

                    if syn_data_pile.is_empty():
                        logs.error(' No synthetic data found in %s. Gain test failed.' 
                                    % os.path.join(data_dir, 'synthetics'))
                        sys.exit()
                
                fband = gainfconf.fband
                taper = trace.CosFader(xfrac=gainfconf.taper_xfrac)

                ag = gainf.AutoGain(data_pile, stations=all_stations,
                                          events=subset_catalog,
                                          arrT=arrT_array,
                                          snr_thresh=gainfconf.snr_thresh,
                                          component=comp,
                                          gain_rel_to=gain_factor_method,
                                          syn_data_pile=syn_data_pile)

                ag.process(fband, taper, twd, gainfconf.debug_mode)

                # Store mean results in YAML format:
                logs.debug(' Saving mean gains: gains_median_and_mean%s.txt %s' % (c, dir_gains))
                ag.save_median_and_mean_and_stdev('gains_median_and_mean%s.txt' % c, directory=dir_gains)
                # Store all results in comma-spread text file:
                ag.save_single_events('gains_all_events%s.txt' % c,
                                      directory=dir_gains, plot=gainfconf.plot_allgains)
                ag = None
            
            for c in gainfconf.components:
                logs.info(' Gain test component %s.' % c)
                run_autogain(data_dir, all_stations, subsets_events['deep'],
                             gainfconf.gain_factor_method,
                             dir_gains, twd, arrT_array, c)

            gc.collect()
            logs.info(' Finished evaluation of gain factors.'
                      + ' Results saved in directory %s.\n'
                      % dir_gains)

        if gainfconf.plot_median_gain_on_map is True:
            logs = logging.getLogger('Gain plot')
            logs.setLevel(verbo)

            dir_gains = os.path.join(data_dir, 'results', 'gains')

            skip_plot = False

            if maps.pl_opt == ['automatic']:
                pl_opt = get_pl_opt(st_lats, st_lons)

            elif len(maps.pl_opt) == 4:
                for o in maps.pl_opt[0:2]:
                    if isinstance(o,str):
                        logs.warning(' Please make sure the %s in the map plotting option in the confic files is set.' % o)
                        skip_plot = True
                if not skip_plot:
                    pl_opt = maps.pl_opt

            else:
                logs.warning(' Set pl_opt to *automatic* or provide [lat, lon, radius, cscale].')
                skip_plot = True

            if not skip_plot:
                for c in gainfconf.components:
                    logs.info(' Plotting gain factors on map, component %s.' % c)
                    plot_median_gain_map_from_file(ns, st_lats, st_lons, pl_opt, maps.pl_topo,
                                                   'gains_median_and_mean%s.txt' % c, dir_gains, c,
                                                   maps.map_size, maps.outformat)

                    logs.info(' Map plot(s) saved in directory %s.\n' % dir_gains)


        ''' 8. Frequency spectra'''
        # psd plot for each station, syn + obs
        # plot ratio of syn and obs for each station, single events
        # plot mean of psd ratio for each station
        # plot flat ranges of psd ratio
        # output flat-ratio-ranges as yaml file

        if psdsconf.calc_psd is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('PSD')
            logs.setLevel(verbo)

            logs.info(' Starting PSD test section.')
            dir_f = os.path.join(data_dir, 'results', 'freq')
            os.makedirs(dir_f, exist_ok=True)

            datapath = os.path.join(data_dir, 'rrd')
            syndatapath = os.path.join(data_dir, 'synthetics')

            logs.info('Data path: %s\nSynthetic data path: %s' % (datapath, syndatapath))

            nst = len(all_stations)

            if arrT_array is None:
                try:
                    data_dir = gensettings.work_dir
                    atfile = os.path.join(data_dir, 'ttt', 'ArrivalTimes_deep.npy')
                    arrT_array = num.load(atfile)
                except:
                    logs.error('Please calculate arrival times first!')
                    raise Exception('Arrival times not calculated!')

            if arrT_R_array is None:
                try:
                    data_dir = gensettings.work_dir
                    atrfile = os.path.join(data_dir, 'ttt', 'ArrivalTimes_estR_deep.npy')
                    arrT_R_array = num.load(atrfile)
                except Exception:
                    logs.error('Please calculate R arrival times first!')
                    raise Exception('R arrival times not calculated!')

            st_numbers = [i_st for i_st in range(len(all_stations))]
            flat_f_ranges_ll = []
            freq_rat_list_y_ll = []
            # flat_by_next_ll = []
            freq_neigh_list_y_ll = []

            nslc_list = []
            # nslc_list2 = []
            # freq_neigh_list_y_ll = []

            for i_st, st in zip(st_numbers, all_stations):
                logs.info('%s %s' % (i_st, st.station))
            # for i_st, st in enumerate(all_stations):
                st_data_pile = pile.make_pile(datapath,
                                              regex='%s_%s_' % (st.network, st.station),
                                              show_progress=False)
                locs = list(set(list(st_data_pile.locations.keys())))
                for l in locs:
                    freq_rat_list_st, freq_rat_list_y, nslc_list_st = fp.prep_psd_fct(
                      i_st, st, nst, l, subsets_events['deep'],
                      dir_f,
                      arrT_array, arrT_R_array,
                      datapath, syndatapath,
                      psdsconf.tinc, psdsconf.tpad,
                      psdsconf.dt_start, psdsconf.dt_end,
                      psdsconf.n_poly,
                      psdsconf.norm_factor,
                      psdsconf.f_ign,
                      plot_psds=psdsconf.plot_psds,
                      plot_ratio_extra=psdsconf.plot_ratio_extra,
                      plot_m_rat=psdsconf.plot_m_rat,
                      plot_flat_ranges=psdsconf.plot_flat_ranges)  # ,
                      # plot_neighb_ranges=psdsconf.plt_neigh_ranges)

                    if freq_rat_list_st != [] and nslc_list_st != []:
                        flat_f_ranges_ll.extend(freq_rat_list_st)
                        freq_rat_list_y_ll.extend(freq_rat_list_y)
                        nslc_list.extend(nslc_list_st)
                '''
                if flat_by_next != [] and nslc_list_st != []:
                    flat_by_next_ll.extend(flat_by_next)
                    freq_neigh_list_y_ll.extend(flat_by_next_y)
                    nslc_list2.extend(nslc_list_st)
                '''

            fp.dump_flat_ranges(flat_f_ranges_ll, freq_rat_list_y_ll,
                                nslc_list, dir_f,
                                fname_ext='linefit2',
                                only_first=psdsconf.only_first)
            '''
            dump_flat_ranges(flat_by_next_ll, freq_neigh_list_y_ll,
                                nslc_list2, dir_f,
                                fname_ext='neighbour2')
            '''
            logs.info(' Finished PSD test section.')
            logs.info(' Results saved in directoriy %s.' % dir_f)


        # 9. Rayleigh wave polarization analysis for orientation
        if orientconf.orient_rayl is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Orientation')
            logs.setLevel(verbo)

            logs.info('Starting Rayleigh wave orientation section.')

            dir_ro = os.path.join(data_dir, 'results', 'orient')
            os.makedirs(dir_ro, exist_ok=True)
            datapath = os.path.join(data_dir, 'rrd')
            list_median_a = []
            list_mean_a = []
            list_stdd_a = []
            list_switched = []
            used_stats = []
            list_all_angles = []
            n_ev = []
            nst = len(all_stations)

            st_numbers = [i_st for i_st in range(nst)]
            for i_st, st in zip(st_numbers, all_stations):
                logs.debug(st.station)
                st_data_pile = pile.make_pile(datapath,
                                              regex='%s_%s_' % (st.network, st.station),
                                              show_progress=False)
                if not st_data_pile:
                    logs.debug('No data found for station %s' % st.station)
                    continue
                locs = list(set(list(st_data_pile.locations.keys())))
                for loc in locs:
                    out = orient.prep_orient(
                                            datapath,
                                            st, i_st, nst,
                                            loc,
                                            subsets_events['shallow'],
                                            dir_ro,
                                            arrTconf.v_rayleigh,
                                            orientconf.bandpass,
                                            orientconf.start_before_ev,
                                            orientconf.stop_after_ev,
                                            orientconf.ccmin,
                                            orientconf.plot_heatmap,
                                            orientconf.plot_distr,
                                            orientconf.debug_mode)

                    if out:
                        list_median_a.append(out[0])
                        list_mean_a.append(out[1])
                        list_stdd_a.append(out[2])
                        list_switched.append(out[3])
                        list_all_angles.append(out[4])
                        n_ev.append(out[5])
                        used_stats.append((st.network, st.station, loc))

            orient.write_output(list_median_a, list_mean_a, list_stdd_a,
                                list_switched,
                                n_ev, used_stats, dir_ro, orientconf.ccmin)

            orient.write_all_output_csv(list_all_angles, used_stats, dir_ro)
            logs.info(' Saved output of orient test in directory %s.' % dir_ro)

        if orientconf.plot_orient_map_fromfile is True:
            logs = logging.getLogger('Orientation')
            logs.setLevel(verbo)
            logs.info(' Plotting output of orient test: Map plot.')
            dir_ro = os.path.join(data_dir, 'results', 'orient')

            skip_plot = False
            if maps.pl_opt == ['automatic']:
                pl_opt = get_pl_opt(st_lats, st_lons)

            elif len(maps.pl_opt) == 4:
                for o in maps.pl_opt[0:2]:
                    if isinstance(o,str):
                        logs.warning(' Please make sure the %s in the map plotting option in the confic files is set.' % o)
                        skip_plot = True

                if not skip_plot:
                    pl_opt = maps.pl_opt

            else:
                logs.warning(' Set pl_opt to *automatic* or provide [lat, lon, radius, cscale].')
                skip_plot = True

            if not skip_plot:
                orient.plot_corr_angles(ns, st_lats, st_lons,
                                        'CorrectionAngles.yaml', dir_ro,
                                        pl_opt, maps.pl_topo,
                                        maps.map_size, maps.outformat,
                                        orientconf.orient_map_label)
            logs.info(' Saved map plot of orient test in directory %s.' % dir_ro)

        if orientconf.plot_angles_vs_events is True:
            logs = logging.getLogger('Orientation')
            logs.setLevel(verbo)
            logs.info(' Plotting output of orient test: Correction angles vs. events.')
            dir_ro = os.path.join(data_dir, 'results', 'orient')
            orient.plot_corr_time(ns, 'AllCorrectionAngles.yaml', dir_ro)
            logs.info(' Saved plot in directory %s.' % dir_ro)

        if orientconf.plot_angles_vs_baz is True:
            dir_ro = os.path.join(data_dir, 'results', 'orient')
            orient.plot_corr_baz(ns, 'AllCorrectionAngles.yaml', 
                                 'CorrectionAngles.yaml', dir_ro, 
                                 subsets_events['shallow'], all_stations)


        if timingconf.timing_test is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Timing')
            logs.setLevel(verbo)

            logs.info(' Starting timing test.')
            if arrT_array is None:
                try:
                    data_dir = gensettings.work_dir
                    arrT_array = num.load(os.path.join(data_dir, 'ttt', 'ArrivalTimes_deep.npy'))
                except:
                    logs.error('Please calculate arrival times first!')
                    raise Exception('Arrival times must be calculated first!')

            subset_catalog = subsets_events['deep']
            nev = len(subset_catalog)
            datapath = os.path.join(gensettings.work_dir, 'rrd')
            syndatapath = os.path.join(gensettings.work_dir, 'synthetics')
            dir_time = os.path.join(gensettings.work_dir, 'results/timing')
            os.makedirs(dir_time, exist_ok=True)
            p_obs = pile.make_pile(datapath, show_progress=False)
            p_syn = pile.make_pile(syndatapath, show_progress=False)

            # print(p_obs, p_syn)

            #if timingconf.search_avail_stats is True:
            #    nslc_list = []
            #    for t in p_obs.iter_traces(trace_selector=lambda tr: tr.channel=='R'):
            #        nslc_list.append(t.nslc_id)
            #    nslc_list = list(set(nslc_list))
            #    stations = nslc_list

            # else:
            
            stations = all_stations

            tshifts = num.empty((len(stations), len(subset_catalog)))
            tshifts.fill(num.nan)

            for i_ev, ev in enumerate(subset_catalog):
                tshifts[:, i_ev] = tt.ccs_allstats_one_event(i_ev, nev, ev, stations, all_stations,
                                                             p_obs, p_syn,
                                                             dir_time, timingconf.bandpass,
                                                             arrT_array, timingconf.cc_thresh,
                                                             timingconf.time_wdw,
                                                             debug_mode=timingconf.debug_mode)
            tshifts_cor = tt.correct_for_med_tshifts(tshifts)
            tt.plot_matrix(tshifts, tshifts_cor, stations, dir_time)

            # get mean and stdev
            medians = num.nanmedian(tshifts_cor, axis=1)
            means = num.nanmean(tshifts_cor, axis=1)
            stdevs = num.nanstd(tshifts_cor, axis=1)

            n_evs = [tshifts_cor.shape[1] - num.isnan(tshifts_cor[i_st,:]).sum()
                     for i_st in range(tshifts_cor.shape[0])]

            # plot
            outfile = os.path.join(dir_time, 'timing_errors_allStats.%s' % maps.outformat)
            tt.plot_tshifts(tshifts_cor, means, stdevs, outfile, stations)
            tt.save_mms(medians, means, stdevs, stations, dir_time, n_evs)

            logs.info(' Finished timing error test.' + 
                      ' Results saved in directory %s.' % dir_time)

        if tc.tele_check is True:
            # Set Logger name and verbosity
            logs = logging.getLogger('Teleseismic')
            logs.setLevel(verbo)

            logs.info('Starting interactive tele-check')
            
            subset_catalog = subsets_events['deep']
            datapath = os.path.join(gensettings.work_dir, 'rrd')
            dir_tc = os.path.join(gensettings.work_dir, 'results/tele_check')
            os.makedirs(dir_tc, exist_ok=True)

            def load_snuffling(win):
                s = TeleCheck()
                s.setup()
                win.pile_viewer.viewer.add_snuffling(s, reloaded=True)

            for ev in subset_catalog:
                ev.name = util.time_to_str(ev.time).replace(' ','_')
                p_obs = pile.make_pile(os.path.join(datapath, ev.name),
                                       show_progress=False)
                if p_obs.is_empty():
                    continue
                p_obs.snuffle(stations=all_stations, events=[ev],
                              launch_hook=load_snuffling)
            
            filename_list = glob.glob(os.path.join(dir_tc, '*.cor'))
            tele.get_correction_statistcs(all_stations, filename_list)

    logging.info('AutoStatsQ run finished.')

if __name__ == '__main__':
    main()
