# https://github.com/gesape/AutoStatsQ - GPLv3
'''
P-wave orientation using the maximum windowed radial/transverse energy ratio.

Configuration is a section within Settings, following OrientConfig::

    - !autostatsq.config.OrientPEnergyConfig
      orient_p_energy: false
      time_wdw: [-5.0, 30.0]
      bandpass: [3.0, 0.03, 0.2]
      energy_ratio_min: 3.0
      snr_min: false  # A float, e.g. 3.0, enables SNR.
      snr_noise_wdw: [-60.0, -10.0]
      resolve_polarity: false
      rz_cc_min: false  # A float enables an inclusive minimum on |CC(R,Z)|.
      plot_angles_vs_events: false
      plot_angles_vs_baz: false
      plot_angles_rose: false

The optional fields may be omitted; SNR filtering and polarity selection
default to disabled. Enable orient_p_energy explicitly to calculate results.

Set orient_p_energy to false and change energy_ratio_min to regenerate
filtered output and plots from saved measurements, without reading waveforms.
Both deep and shallow event subsets are used. 
Energy cannot distinguish angles separated by 180 degrees: energy-only 
corrections are in [-90, 90).

Optional polarity selection uses demeaned zero-lag R-Z correlation and reports
corrections in [-180, 180).
With resolve_polarity enabled, optional rz_cc_min filters saved measurements
by abs(CC(R,Z)) >= rz_cc_min. Changing this threshold requires no waveforms.
Resolved-result filenames include rzcc followed by the threshold or -off.

Optional SNR is vertical RMS in time_wdw divided by RMS in snr_noise_wdw.
All measured events are saved before applying the energy and SNR thresholds.
Changing either threshold needs no waveform processing. Enabling a feature
that was disabled during measurement, or changing polarity/noise-window
settings, requires recalculation. 
'''

import logging
import os

import numpy as num
from pyrocko import orthodrome, pile, trace, util
from pyrocko.guts import Any, Bool, Dict, Float, String, Tuple, load
from . import orient


RAW_FILENAME = 'AllCorrectionAngles_erall.yaml'
ROT_ANGLES = range(-90, 90)
FILTER_TPAD = 300.0  # Seconds before and after the analysis window.


class P_energy_by_stat(orient.rota_ev_by_stat):
    ev_energy_ratio = Dict.T(String.T(), Float.T())
    ev_snr = Dict.T(String.T(), Float.T(), default={})
    ev_rz_corr = Dict.T(String.T(), Float.T(), default={})


class P_energy_all(orient.dict_stats_all_rota):
    time_wdw = Tuple.T(2, Float.T())
    bandpass = Tuple.T(3, Float.T())
    arrival_tables = Dict.T(String.T(), String.T())
    snr_min = Any.T(default=False)
    snr_noise_wdw = Tuple.T(2, Float.T(), default=(-60.0, -10.0))
    resolve_polarity = Bool.T(default=False)
    rz_cc_min = Any.T(default=False)


def check_settings(conf):
    start, stop = conf.time_wdw
    order, fmin, fmax = conf.bandpass
    if not all(num.isfinite([start, stop, order, fmin, fmax])):
        raise ValueError('P-energy window and bandpass must be finite.')
    if start >= stop:
        raise ValueError('P-energy window start must precede its end.')
    if order < 1 or order != int(order) or not 0 < fmin < fmax:
        raise ValueError('P-energy bandpass: [positive integer order, fmin, fmax].')
    if not num.isfinite(conf.energy_ratio_min) or conf.energy_ratio_min < 1:
        raise ValueError('P-energy ratio threshold must be finite and >= 1.')

    if isinstance(conf.snr_min, float):
        if not num.isfinite(conf.snr_min) or conf.snr_min <= 0:
            raise ValueError('P-energy SNR threshold must be positive and finite.')
        noise_start, noise_stop = conf.snr_noise_wdw
        if (not all(num.isfinite([noise_start, noise_stop]))
                or not noise_start < noise_stop <= 0):
            raise ValueError('SNR noise window must have start < end <= P arrival.')
    if conf.resolve_polarity is True and isinstance(conf.rz_cc_min, float):
        if not num.isfinite(conf.rz_cc_min) or not 0.0 <= conf.rz_cc_min <= 1.0:
            raise ValueError('R-Z correlation minimum must be finite and between 0.0 and 1.0.')


def get_vertical_windows(st_data_pile, st, loc, arrT, conf):
    '''Filter padded Z before chopping the P signal and optional noise windows.'''
    start_twd, end_twd = [arrT + dt for dt in conf.time_wdw]
    filter_start, filter_stop = start_twd, end_twd
    
    if isinstance(conf.snr_min, float):
        noise_start, noise_stop = [arrT + dt for dt in conf.snr_noise_wdw]
        filter_start = min(filter_start, noise_start)
        filter_stop = max(filter_stop, noise_stop)
    traces = st_data_pile.all(
        tmin=filter_start-FILTER_TPAD, tmax=filter_stop+FILTER_TPAD,
        trace_selector=lambda tr: tr.nslc_id == (st.network, st.station, loc, 'Z'),
        want_incomplete=False)
    if len(traces) != 1:
        raise ValueError('Missing, incomplete or multiple Z traces.')
    trZ = traces[0].copy()
    order, fmin, fmax = conf.bandpass
    if fmax >= 0.5 / trZ.deltat:
        raise ValueError('Vertical bandpass maximum must be below Nyquist.')
    if not num.all(num.isfinite(trZ.ydata)):
        raise ValueError('Nonfinite vertical waveform samples.')
    trZ.bandpass(int(order), fmin, fmax)
    snr = None

    if isinstance(conf.snr_min, float):
        noise = trZ.chop(tmin=noise_start, tmax=noise_stop, inplace=False)
        noise_rms = num.sqrt(num.mean(noise.ydata**2))
    trZ.chop(tmin=start_twd, tmax=end_twd)

    if isinstance(conf.snr_min, float):
        signal_rms = num.sqrt(num.mean(trZ.ydata**2))
        snr = float(signal_rms / noise_rms) if noise_rms > 0 else num.inf
    return trZ, snr


def get_polarity_angle(trR, trT, trZ, angle):
    '''Choose positive zero-lag R-Z correlation for prepared Z-up/R-away data.'''
    if (not num.isclose(trR.deltat, trZ.deltat)
            or abs(trR.tmin-trZ.tmin) > trR.deltat * 1.e-4
            or trR.ydata.size != trZ.ydata.size):
        raise ValueError('R/Z sample times or lengths do not match.')
    rot_R, rot_T = trace.rotate(
        [trR, trT], azimuth=angle,
        in_channels=['R', 'T'], out_channels=['2', '3'])
    r = rot_R.ydata.astype(float)
    z = trZ.ydata.astype(float)

    rscale, zscale = num.max(num.abs(r)), num.max(num.abs(z))
    if not num.isfinite(rscale) or not num.isfinite(zscale) or min(rscale, zscale) == 0:
        raise ValueError('No finite R/Z energy for polarity selection.')
    r, z = r / rscale, z / zscale
    r, z = r - num.mean(r), z - num.mean(z)
    denom = num.sqrt(num.dot(r, r) * num.dot(z, z))
    
    if not num.isfinite(denom) or denom == 0:
        raise ValueError('Undefined R-Z correlation for polarity selection.')
    
    corr = float(num.clip(num.dot(r, z) / denom, -1.0, 1.0))
    if corr == 0:
        raise ValueError('Zero R-Z correlation cannot select a polarity.')
    # if corr * expected_rz_sign < 0:
    if corr < 0:
        angle = (angle + 180) % 360
        if angle >= 180:
            angle -= 360
        corr = -corr
    return int(angle), corr


def get_max_energy_ratio(trR, trT):
    '''Find the best integer correction angle, as in the Rayleigh test.'''

    scale = max(num.max(num.abs(trR.ydata)), num.max(num.abs(trT.ydata)))

    if not num.isfinite(scale) or scale == 0:
        raise ValueError('No finite horizontal energy.')

    trR, trT = trR.copy(), trT.copy()
    trR.set_ydata(trR.ydata.astype(float) / scale)
    trT.set_ydata(trT.ydata.astype(float) / scale)
    ratios = []

    for angle in ROT_ANGLES:
        rot_R, rot_T = trace.rotate(
            [trR, trT], azimuth=angle,
            in_channels=['R', 'T'], out_channels=['2', '3'])

        er = float(num.dot(rot_R.ydata, rot_R.ydata))
        et = float(num.dot(rot_T.ydata, rot_T.ydata))
        ratios.append(er / et if et > 0 else num.inf)

    # Prefer the smallest absolute correction for tied maxima.
    best = max(range(len(ratios)),
               key=lambda i: (ratios[i], -abs(ROT_ANGLES[i])))

    return int(ROT_ANGLES[best]), float(ratios[best])


def prep_orient(st_data_pile, st, i_st, loc, catalog, arrT_array, conf):
    '''Filter padded R/T data, chop about the saved P arrival, maximize E_R/E_T.'''
    
    logs = logging.getLogger('prep_orient_P_energy')
    dict_ev_angle, dict_ev_ratio = {}, {}
    dict_ev_snr, dict_ev_rz_corr = {}, {}
    order, fmin, fmax = conf.bandpass
    tpad = FILTER_TPAD

    for i_ev, ev in enumerate(catalog):
        arrT = arrT_array[i_ev, i_st]
        
        if not num.isfinite(arrT):
            logs.warning('No P arrival for %s.%s, %s.',
                         st.network, st.station, util.tts(ev.time))
            continue

        start_twd, end_twd = [arrT + dt for dt in conf.time_wdw]

        try:
            horizontals = []
            for cha in ['R', 'T']:
                traces = st_data_pile.all(
                    tmin=start_twd-tpad, tmax=end_twd+tpad,
                    trace_selector=lambda tr: (
                        tr.nslc_id == (st.network, st.station, loc, cha)),
                    want_incomplete=False)

                if len(traces) != 1:
                    raise ValueError('Missing, incomplete or multiple %s traces.' % cha)

                tr = traces[0].copy()

                if fmax >= 0.5 / tr.deltat:
                    raise ValueError('Bandpass maximum must be below Nyquist.')

                tr.bandpass(int(order), fmin, fmax)
                tr.chop(tmin=start_twd, tmax=end_twd)
                horizontals.append(tr)

            trR, trT = horizontals
            if (not num.isclose(trR.deltat, trT.deltat)
                    or abs(trR.tmin-trT.tmin) > trR.deltat * 1.e-4
                    or trR.ydata.size != trT.ydata.size
                    or trR.ydata.size < 2):
                raise ValueError('R/T sample times or lengths do not match.')

            #trace.snuffle([trR,trT])

            snr, rz_corr = None, None
            # if conf.use_snr is True or conf.resolve_polarity is True:
            if isinstance(conf.snr_min, float) or conf.resolve_polarity is True:
                trZ, snr = get_vertical_windows(st_data_pile, st, loc, arrT, conf)
            angle, ratio = get_max_energy_ratio(trR, trT)
            if conf.resolve_polarity is True:
                angle, rz_corr = get_polarity_angle(
                    trR, trT, trZ, angle)

        except (ValueError, trace.NoData, trace.TraceTooShort) as exc:
            logs.warning('Skipping %s.%s.%s, %s: %s',
                         st.network, st.station, loc, util.tts(ev.time), exc)
            continue
        
        ev_time = util.time_to_str(ev.time)
        dict_ev_angle[ev_time] = angle
        dict_ev_ratio[ev_time] = ratio
        if isinstance(conf.snr_min, float):
            dict_ev_snr[ev_time] = snr
        if conf.resolve_polarity is True:
            dict_ev_rz_corr[ev_time] = rz_corr

    return dict_ev_angle, dict_ev_ratio, dict_ev_snr, dict_ev_rz_corr



def get_result_tag(conf):
    '''Filename suffix for ER, SNR and optional polarity-correlation filtering.'''
    snr = conf.snr_min if isinstance(conf.snr_min, float) else '-off'
    tag = 'er%s_snr%s' % (conf.energy_ratio_min, snr)
    if conf.resolve_polarity is True:
        cc = conf.rz_cc_min if isinstance(conf.rz_cc_min, float) else '-off'
        tag += '_rzcc%s' % cc
    return tag


def write_all_output_yaml(list_all_angles, list_all_ratios, used_stats,
                          dir_ro, threshold, conf, arrival_tables,
                          list_all_snr, list_all_rz_corr, selection_conf=None):
    '''Reuse the existing writer, extending each station with energy ratios.'''
    orient.write_all_output_yaml(
        list_all_angles, used_stats, dir_ro, threshold)
    selection = conf if selection_conf is None else selection_conf
    if threshold == 'all':
        filename = os.path.join(dir_ro, RAW_FILENAME)
    else:
        filename = os.path.join(
            dir_ro, 'AllCorrectionAngles_%s.yaml' % get_result_tag(selection))
    os.replace(os.path.join(dir_ro, 'AllCorrectionAngles_cc%s.yaml' % threshold),
               filename)
    saved = load(filename=filename)

    results = P_energy_all(
        dict_stats_all=[P_energy_by_stat(
            station=item.station, ev_rota=item.ev_rota, ev_energy_ratio=ratios,
            ev_snr=snr, ev_rz_corr=rz_corr)
            for item, ratios, snr, rz_corr in zip(
                saved.dict_stats_all, list_all_ratios, list_all_snr, list_all_rz_corr)],
        time_wdw=conf.time_wdw, bandpass=conf.bandpass,
        arrival_tables=arrival_tables,
        snr_min=selection.snr_min,
        snr_noise_wdw=conf.snr_noise_wdw,
        resolve_polarity=conf.resolve_polarity,
        rz_cc_min=selection.rz_cc_min if conf.resolve_polarity is True else False)
    results.regularize()
    results.validate()
    results.dump(filename=filename)


def get_m_angle(angles, period=180.0):
    '''Mean, median and RMS angular deviation for axial or resolved angles.'''
    if not angles:
        return num.nan, num.nan, num.nan, 0
    angles = num.asarray(angles, dtype=float)
    factor = 360.0 / period
    half_period = period / 2.0
    vector = num.mean(num.exp(1j * factor * num.deg2rad(angles)))
    if abs(vector) < 1.e-12:
        return num.nan, num.nan, num.nan, len(angles)
    mean_a = (num.rad2deg(num.angle(vector)) / factor + half_period) % period - half_period
    deviations = (angles - mean_a + half_period) % period - half_period
    median_a = (mean_a + num.median(deviations) + half_period) % period - half_period
    std_a = num.sqrt(num.mean(deviations**2)) if len(angles) > 1 else num.nan
    return float(median_a), float(mean_a), float(std_a), len(angles)


def reanalyse_saved_results(dir_ro, conf):
    '''Apply current ER/SNR/R-Z thresholds to saved measurements and save results.'''
    saved = load(filename=os.path.join(dir_ro, RAW_FILENAME))
    if isinstance(conf.snr_min, float):
        if not isinstance(saved.snr_min, float):
            raise ValueError('Saved results have no SNR measurements; rerun with a float snr_min.')
        if not num.isfinite(conf.snr_min) or conf.snr_min <= 0:
            raise ValueError('P-energy SNR threshold must be positive and finite.')
        if conf.snr_noise_wdw != saved.snr_noise_wdw:
            raise ValueError('Noise window has changed; rerun to measure SNR for that window.')

    if conf.resolve_polarity != saved.resolve_polarity:
        raise ValueError('Polarity settings have changed; rerun the P-energy calculation.')
    if conf.resolve_polarity is True and isinstance(conf.rz_cc_min, float):
        if not num.isfinite(conf.rz_cc_min) or not 0.0 <= conf.rz_cc_min <= 1.0:
            raise ValueError('R-Z correlation minimum must be finite and between 0.0 and 1.0.')
    
    period = 360.0 if saved.resolve_polarity else 180.0
    list_snr, list_rz_corr = [], []
    used_stats, list_angles, list_ratios = [], [], []
    list_median_a, list_mean_a, list_stdd_a, n_ev = [], [], [], []
    
    for item in saved.dict_stats_all:
        ratios = {ev: ratio for ev, ratio in item.ev_energy_ratio.items()
                  if ratio >= conf.energy_ratio_min}
    
        if isinstance(conf.snr_min, float):
            if any(ev not in item.ev_snr for ev in ratios):
                raise ValueError('Missing saved SNR values; rerun the P-energy calculation.')
            ratios = {ev: ratio for ev, ratio in ratios.items()
                      if item.ev_snr[ev] >= conf.snr_min}
    
        if conf.resolve_polarity is True and isinstance(conf.rz_cc_min, float):
            if any(ev not in item.ev_rz_corr for ev in ratios):
                raise ValueError('Missing saved R-Z correlations; rerun with resolve_polarity enabled.')
            ratios = {ev: ratio for ev, ratio in ratios.items()
                      if num.isfinite(item.ev_rz_corr[ev])
                      and abs(item.ev_rz_corr[ev]) >= conf.rz_cc_min}
    
        angles = {ev: angle for ev, angle in item.ev_rota.items() if ev in ratios}
        median_a, mean_a, std_a, nev = get_m_angle(list(angles.values()), period)
        list_snr.append({ev: value for ev, value in item.ev_snr.items() if ev in angles})
        list_rz_corr.append({ev: value for ev, value in item.ev_rz_corr.items() if ev in angles})
        used_stats.append(item.station)
        list_angles.append(angles)
        list_ratios.append(ratios)
        list_median_a.append(median_a)
        list_mean_a.append(mean_a)
        list_stdd_a.append(std_a)
        n_ev.append(nev)
    write_all_output_yaml(list_angles, list_ratios, used_stats, dir_ro,
                          conf.energy_ratio_min, saved, saved.arrival_tables,
                          list_snr, list_rz_corr, selection_conf=conf)
    orient.write_output(list_median_a, list_mean_a, list_stdd_a, [],
                        n_ev, used_stats, dir_ro, conf.energy_ratio_min)
    if used_stats:
        os.replace(
            os.path.join(dir_ro, 'CorrectionAngles_cc%s.yaml' % conf.energy_ratio_min),
            os.path.join(dir_ro, 'CorrectionAngles_%s.yaml' % get_result_tag(conf)))
    return used_stats


def plot_corr_rose(nsl, filename_all, filename_stats, dir_ro, threshold,
                   bins_deg=1.0, dpi=100, snr_min=False,
                   resolve_polarity=False, rz_cc_min=False):
    '''Reuse the shared rose plotter with an energy-ratio caption.'''
    original_subplots = orient.plt.subplots
    selection_text = 'ER ≥ %s' % threshold
    if isinstance(snr_min, float):
        selection_text += ', SNR ≥ %s' % snr_min
    else:
        selection_text += ', SNR off'
    if resolve_polarity is True:
        if isinstance(rz_cc_min, float):
            selection_text += ', |CC(R,Z)| ≥ %s' % rz_cc_min
        else:
            selection_text += ', R-Z CC filter off'

    def subplots_er(*args, **kwargs):
        fig, ax = original_subplots(*args, **kwargs)
        original_savefig = fig.savefig

        def savefig_er(*save_args, **save_kwargs):
            ax.set_title(ax.get_title().replace(
                'CC > %s' % threshold, selection_text), pad=28)
            return original_savefig(*save_args, **save_kwargs)

        fig.savefig = savefig_er
        return fig, ax

    try:
        orient.plt.subplots = subplots_er
        orient.plot_corr_rose(
            nsl, filename_all, filename_stats, dir_ro, threshold,
            bins_deg=bins_deg, dpi=dpi)
    finally:
        orient.plt.subplots = original_subplots


def plot_output(dir_ro, conf, events, stations, maps, nsl):
    '''Only saved YAML files are used for angles and summary statistics.'''
    logs = logging.getLogger('Orientation_P_energy')
    threshold = conf.energy_ratio_min
    # fn_stats = 'CorrectionAngles_er%s.yaml' % threshold
    result_tag = get_result_tag(conf)
    fn_all = 'AllCorrectionAngles_%s.yaml' % result_tag
    fn_stats = 'CorrectionAngles_%s.yaml' % result_tag
    jobs = []
    
    if conf.plot_angles_vs_events:
        # jobs.append((orient.plot_corr_time, (nsl, fn_all, dir_ro, threshold), {}))
        jobs.append((orient.plot_corr_time, (nsl, fn_all, dir_ro, result_tag), {}))
    
    if conf.plot_angles_vs_baz and events:
        jobs.append((orient.plot_corr_baz,
                     (nsl, fn_all, fn_stats, dir_ro, events, stations, result_tag), {}))
    
    if conf.plot_angles_rose:
        jobs.append((plot_corr_rose,
                     (nsl, fn_all, fn_stats, dir_ro, threshold),
                     dict(bins_deg=conf.rose_bins_deg, dpi=conf.rose_dpi,
                          snr_min=conf.snr_min,
                          resolve_polarity=conf.resolve_polarity,
                          rz_cc_min=conf.rz_cc_min)))
    
    if conf.plot_orient_map_fromfile:
        pl_opt = maps.pl_opt
        if pl_opt == ['automatic'] and stations:
            lats = [s.lat for s in stations]
            lons = [s.lon for s in stations]
            lat_m, lon_m = float(num.mean(lats)), float(num.mean(lons))
            corners = [orthodrome.Loc(min(lats), min(lons)),
                       orthodrome.Loc(max(lats), max(lons))]
            radius = 1.2 * max(orthodrome.distance_accurate50m(
                orthodrome.Loc(lat_m, lon_m), corner) for corner in corners)
            pl_opt = [lat_m, lon_m, max(radius, 1000.0), 'split']
    
        jobs.append((orient.plot_corr_angles,
                     ([(s.network, s.station) for s in stations],
                      [s.lat for s in stations], [s.lon for s in stations],
                      fn_stats, dir_ro, pl_opt, maps.pl_topo,
                      maps.map_size, maps.outformat, result_tag, conf.orient_map_label), {}))
    
    for function, args, kwargs in jobs:
        try:
            function(*args, **kwargs)
            if function == plot_corr_rose:
                old_suffix = '_rose_cc%s.png' % threshold
                new_suffix = '_rose_%s.png' % result_tag
                for filename in os.listdir(dir_ro):
                    if filename.endswith(old_suffix):
                        new_filename = filename[:-len(old_suffix)] + new_suffix
                        os.replace(os.path.join(dir_ro, filename),
                                   os.path.join(dir_ro, new_filename))
        except Exception:
            logs.exception('P-energy plot failed (%s); saved results are available.',
                           function.__name__)


def run_orient(data_dir, stations, subsets_events, arrT_by_subset, conf):
    '''Calculate and save P-energy results using existing P-arrival tables.'''
    
    check_settings(conf)
    
    logs = logging.getLogger('Orientation_P_energy')
    dir_ro = os.path.join(data_dir, 'results', 'orient_P-energy')

    if conf.orient_p_energy is True:
        os.makedirs(dir_ro, exist_ok=True)
        
        datapath = os.path.join(data_dir, 'rrd')
        
        arrival_tables = {}
        for subset in ['deep', 'shallow']:
            expected_shape = (len(subsets_events[subset]), len(stations))
            
            if arrT_by_subset[subset].shape != expected_shape:
                raise ValueError(
                    'P-arrival table for %s has shape %s; expected %s. '
                    'Use the same event subsets and station order as for '
                    'arrival-table generation.'
                    % (subset, arrT_by_subset[subset].shape, expected_shape))

            arrival_tables[subset] = os.path.join(
                data_dir, 'ttt', 'ArrivalTimes_%s.npy' % subset)

        used_stats, list_all_angles, list_all_ratios = [], [], []
        list_all_snr, list_all_rz_corr = [], []
        for i_st, st in enumerate(stations):
            st_data_pile = pile.make_pile(
                datapath, regex='%s_%s_' % (st.network, st.station),
                show_progress=False)
            
            for loc in sorted(st_data_pile.locations):
                angles, ratios = {}, {}
                snr, rz_corr = {}, {}
                for subset in ['deep', 'shallow']:
                    ev_angles, ev_ratios, ev_snr, ev_rz_corr = prep_orient(
                        st_data_pile, st, i_st, loc, subsets_events[subset],
                        arrT_by_subset[subset], conf)
                    snr.update(ev_snr)
                    rz_corr.update(ev_rz_corr)
                    angles.update(ev_angles)
                    ratios.update(ev_ratios)

                used_stats.append((st.network, st.station, loc))
                list_all_angles.append(angles)
                list_all_ratios.append(ratios)
                list_all_snr.append(snr)
                list_all_rz_corr.append(rz_corr)

        write_all_output_yaml(list_all_angles, list_all_ratios, used_stats,
                              dir_ro, 'all', conf, arrival_tables,
                              list_all_snr, list_all_rz_corr)
        
        reanalyse_saved_results(dir_ro, conf)
        
        logs.info('Saved P-energy results in %s.', dir_ro)
