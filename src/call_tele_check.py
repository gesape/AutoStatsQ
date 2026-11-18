# https://github.com/gesape/AutoStatsQ - GPLv3
#
# AutoStasQ is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
#
# The AutoStatsQ Developers, 21st Century
# ---|P------/S----------~Lg----------

'''
File containing functionality to run interactive tele check.
'''

import grond
import logging
import numpy as num
from pyrocko.guts import Object, List, Int, Tuple, Float, String, dump_all
import datetime
import os
import tempfile
from pyrocko import guts, util

# Parameters for saved TeleCheck orientation results and plots.
polarization_filename = 'PolarizationAngles.yaml'
relative_filename = 'RelativeAlignmentAngles.yaml'
gain_filename = 'GainFactors.yaml'
rotation_failure_threshold = 3.0
rotation_step_deg = 1.0
max_alignment_iterations = 10000
rose_bin_deg = 5.0
plot_dpi = 100
angle_decimals = 1
linearity_decimals = 2
yaml_decimals = 2
# Fixed marker areas in points squared for linearity 0 and 1.
linearity_symbol_area_min = 8.0
linearity_symbol_area_max = 80.0
linearity_legend_values = (0.0, 0.5, 1.0)
linearity_symbol_exponent = 3.0


def _rounded_yaml_output(value):
    """Copy YAML data with floats rounded to the output precision."""

    if isinstance(value, (float, num.floating)):
        return round(float(value), yaml_decimals)

    if isinstance(value, Object):
        return type(value)(**{
            prop.name: _rounded_yaml_output(getattr(value, prop.name))
            for prop in value.T.properties})

    if isinstance(value, list):
        return [_rounded_yaml_output(item) for item in value]

    if isinstance(value, tuple):
        return tuple(_rounded_yaml_output(item) for item in value)

    if isinstance(value, dict):
        return {key: _rounded_yaml_output(item) for key, item in value.items()}

    return value


class OrientationEvent(Object):
    event_name = String.T()
    angle_deg = Float.T(optional=True,
                        help='Absent for unavailable or rejected estimates.')
    linearity = Float.T(optional=True)
    relative_waveform_misfit = Float.T(optional=True)


class OrientationSummary(Object):
    median_angle_deg = Float.T(optional=True)
    mean_angle_deg = Float.T(optional=True)
    scatter_deg = Float.T(optional=True,
        help='RMS wrapped angular deviation about the mean, in degrees.')
    n_events = Int.T()
    linearity_weighted_mean_angle_deg = Float.T(optional=True)
    linearity_weighted_median_angle_deg = Float.T(optional=True)
    linearity_weighted_stddev_deg = Float.T(optional=True,
        help='Linearity-weighted RMS wrapped deviation about the weighted mean.')


class OrientationStation(Object):
    station = Tuple.T(3, String.T())
    summary = OrientationSummary.T()
    events = List.T(OrientationEvent.T())


class OrientationResults(Object):
    method = String.T()
    angle_units = String.T(default='degrees')
    angular_period_deg = Float.T()
    dict_stats_all = List.T(OrientationStation.T())


class GainEvent(Object):
    event_name = String.T()
    channel = String.T()
    gain_factor = Float.T(optional=True,
                          help='Relative amplitude factor; absent when unavailable.')


class GainSummary(Object):
    channel = String.T(help='Analyzed component, e.g. Z.')
    median_factor = Float.T(optional=True)
    mean_factor = Float.T(optional=True)
    stdev_factor = Float.T(optional=True)
    n_events = Int.T()


class GainStation(Object):
    station = Tuple.T(3, String.T())
    summary = List.T(GainSummary.T())
    events = List.T(GainEvent.T())


class GainResults(Object):
    method = String.T(default='Cross-correlation relative gain')
    factor_definition = String.T(
        default='Relative amplitude factor as in .cor; the displayed correction is its reciprocal.')
    dict_stats_all = List.T(GainStation.T())


def save_gain_results(directory, event_time, stations, factors, components):
    """Save station/event gain factors and component-wise summaries."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, gain_filename)
    result = guts.load(filename=path) if os.path.exists(path) else GainResults(dict_stats_all=[])

    if not isinstance(result, GainResults):
        raise ValueError('Not a TeleCheck gain result document.')

    saved = {tuple(item.station): item for item in result.dict_stats_all}
    event_name = event_time.replace(' ', '_', 1)

    for item in saved.values():
        item.events = [event for event in item.events if event.event_name != event_name]

    for station in stations:
        station = tuple(station)
        item = saved.setdefault(station, GainStation(station=station, summary=[], events=[]))
        measured_components = set()

        for nslc, factor in sorted(factors.items()):
            if tuple(nslc[:3]) != station:
                continue

            measured_components.add(nslc[3][-1])
            values = dict(event_name=event_name, channel=nslc[3])

            if num.isfinite(factor) and factor > 0.:
                values['gain_factor'] = float(factor)

            item.events.append(GainEvent(**values))

        for component in components:
            if component not in measured_components:
                item.events.append(GainEvent(event_name=event_name, channel=component))

    for item in saved.values():
        item.events.sort(key=lambda event: (event.event_name, event.channel))
        item.summary = []

        for component in sorted({event.channel[-1] for event in item.events}):
            values = num.asarray([
                event.gain_factor for event in item.events
                if event.channel.endswith(component) and event.gain_factor is not None
                and num.isfinite(event.gain_factor) and event.gain_factor > 0.], dtype=float)
            summary = GainSummary(channel=component, n_events=int(values.size))

            if values.size:
                summary.median_factor = float(num.median(values))
                summary.mean_factor = float(num.mean(values))
                summary.stdev_factor = float(num.std(values))

            item.summary.append(summary)

    result.dict_stats_all = [saved[station] for station in sorted(saved)]
    result.regularize()
    result.validate()
    descriptor, temporary = tempfile.mkstemp(prefix='.' + gain_filename, dir=directory)
    os.close(descriptor)

    try:
        _rounded_yaml_output(result).dump(filename=temporary)
        os.replace(temporary, path)
    finally:

        if os.path.exists(temporary):
            os.unlink(temporary)

    logging.info('Saved TeleCheck event gain results: %s', path)


def _orientation_mapping(saved):
    """Convert Guts objects or a station mapping to internal measurements."""

    if isinstance(saved, dict):
        return saved

    if not isinstance(saved, OrientationResults):
        raise ValueError('Not a TeleCheck orientation result document.')

    stations = {}

    for item in saved.dict_stats_all:
        events = {}
        for event in item.events:
            usable = event.angle_deg is not None and num.isfinite(event.angle_deg)
            events[event.event_name.replace('_', ' ', 1)] = dict(
                angle_deg=event.angle_deg, usable=bool(usable),
                linearity=event.linearity,
                relative_waveform_misfit=event.relative_waveform_misfit)

        stations['.'.join(item.station)] = dict(
            events=events,
            summary={name: getattr(item.summary, name) for name in (
                'median_angle_deg', 'mean_angle_deg', 'scatter_deg', 'n_events',
                'linearity_weighted_mean_angle_deg',
                'linearity_weighted_median_angle_deg', 'linearity_weighted_stddev_deg')})

    return stations


def _orientation_document(stations, period):
    # """Make a compact, typed station/event document; diagnostics stay internal."""
    """Build typed station/event results with rounded measurements."""
    station_results = []

    for key, item in sorted(stations.items()):
        events = []

        for time, measurement in sorted(item['events'].items()):
            angle = measurement.get('angle_deg')
            usable = (measurement.get('usable', False) and angle is not None
                      and num.isfinite(angle))
            values = dict(event_name=time.replace(' ', '_', 1))

            if usable:
                values['angle_deg'] = round(float(angle), angle_decimals)

            linearity = measurement.get('linearity')

            if period == 180. and linearity is not None and num.isfinite(linearity):
                values['linearity'] = round(float(linearity), linearity_decimals)

            misfit = measurement.get('relative_waveform_misfit')

            if period == 360. and misfit is not None and num.isfinite(misfit):
                values['relative_waveform_misfit'] = float(misfit)

            events.append(OrientationEvent(**values))

        station_results.append(OrientationStation(
            station=tuple(key.split('.')), events=events,
            summary=OrientationSummary(**item['summary'])))

    result = OrientationResults(
        method='P polarization' if period == 180. else 'Relative network alignment',
        angular_period_deg=float(period), dict_stats_all=station_results)
    result.regularize()
    result.validate()

    return result


def polarization_measurement(x, y):
    """Calculate polarization angle and linearity for x=T, y=R."""
    result = dict(angle_deg=None, linearity=None, usable=False,
                  status='invalid_signal')
    data = num.asarray([x, y], dtype=float)

    if data.shape[1] < 2 or not num.isfinite(data).all():
        return result

    eigenvalues, eigenvectors = num.linalg.eigh(num.cov(data))

    if not num.isfinite(eigenvalues).all() or eigenvalues[-1] <= 0.0:
        return result

    linearity = float(num.clip(1.0 - eigenvalues[0] / eigenvalues[-1], 0., 1.))
    result['linearity'] = linearity

    if num.isclose(eigenvalues[0], eigenvalues[-1], rtol=1.e-12, atol=0.):
        result['status'] = 'undefined_axis'
        return result

    pc = eigenvectors[:, -1]
    result.update(
        angle_deg=float((-num.degrees(num.arctan2(pc[1], pc[0]))) % 180. - 90.),
        usable=True, status='ok')

    return result


def relative_orientation(nsls, xs, ys):
    # """Existing pairwise waveform-rotation method, with final accumulated angles."""
    """Calculate relative network rotations and waveform misfits."""
    results = {nsl: dict(angle_deg=None, relative_waveform_misfit=None,
                         usable=False, status='insufficient_stations')
               for nsl in nsls}
    n = len(nsls)

    if n < 2:
        return results

    amplitudes = num.max(num.hypot(xs, ys), axis=1)
    valid = (num.isfinite(xs).all(axis=1) & num.isfinite(ys).all(axis=1)
             & num.isfinite(amplitudes) & (amplitudes > 0.))

    if not valid.all():
        valid_indices = num.flatnonzero(valid)
        subset = [nsls[i] for i in valid_indices]

        for nsl in nsls:
            results[nsl]['status'] = 'invalid_signal'

        results.update(relative_orientation(subset, xs[valid], ys[valid]))

        return results

    xs = xs / amplitudes[:, None]
    ys = ys / amplitudes[:, None]
    phis = num.arange(-180., 180. + rotation_step_deg / 2., rotation_step_deg)
    d = num.empty((n, n, phis.size))

    for ia in range(n):

        for iphi, phi in enumerate(num.deg2rad(phis)):
            xrot = num.cos(phi) * xs[ia] + num.sin(phi) * ys[ia]
            yrot = -num.sin(phi) * xs[ia] + num.cos(phi) * ys[ia]
            d[ia, :, iphi] = num.sqrt(num.sum(
                (xrot[None, :] - xs)**2 + (yrot[None, :] - ys)**2, axis=1))

    dmins = num.min(d, axis=2)
    denominator = float(num.median(dmins))

    if denominator > 0. and num.isfinite(denominator):
        scores = num.mean(dmins / denominator, axis=1)
    elif num.all(dmins == 0.):
        scores = num.zeros(n)
    else:
        for result in results.values():
            result['status'] = 'undefined_misfit'

        return results

    # failed = scores > ROTATION_FAILURE_THRESHOLD
    failed = scores > rotation_failure_threshold
    phimins = phis[num.argmin(d, axis=2)]
    rotations = num.zeros(n)
    converged = False

    for _ in range(max_alignment_iterations):
        worst = num.argmax(num.mean(num.abs(phimins), axis=1))
        phimod = (phis[None, :] + phimins[worst, :, None] + 180.) % 360. - 180.
        update = phis[num.argmin(num.mean(num.abs(phimod), axis=0))]

        if abs(update) < rotation_step_deg:
            converged = True
            break

        phimins[worst, :] = (phimins[worst, :] + update + 180.) % 360. - 180.
        phimins[:, worst] = (phimins[:, worst] - update + 180.) % 360. - 180.

        if not failed[worst]:
            rotations[worst] += update

    for i, nsl in enumerate(nsls):
        usable = bool(converged and not failed[i] and num.isfinite(scores[i]))
        results[nsl].update(
            angle_deg=float((rotations[i] + 180.) % 360. - 180.) if usable else None,
            relative_waveform_misfit=float(scores[i]), usable=usable,
            status=('not_converged' if not converged else
                    'waveform_mismatch' if failed[i] else 'ok'))

    return results


def _linearity_weighted_summary(events, period):
    """Calculate linearity-weighted angular mean, median and scatter."""
    summary = dict(linearity_weighted_mean_angle_deg=None,
                   linearity_weighted_median_angle_deg=None,
                   linearity_weighted_stddev_deg=None)
    pairs = []

    for measurement in events.values():
        angle = measurement.get('angle_deg')
        weight = measurement.get('linearity')

        if (measurement.get('usable', False) and angle is not None
                and weight is not None and num.isfinite(angle)
                and num.isfinite(weight) and 0. < weight <= 1.):
            pairs.append((angle, weight))

    if not pairs:
        return summary

    angles, weights = num.asarray(pairs, dtype=float).T
    total_weight = num.sum(weights)
    vector = num.sum(weights * num.exp(2j * num.pi * angles / period)) / total_weight

    if abs(vector) < 1.e-12:
        return summary

    half = period / 2.
    mean = (float(num.angle(vector)) * period / (2. * num.pi) + half) % period - half
    deviations = (angles - mean + half) % period - half
    order = num.argsort(deviations)
    cumulative = num.cumsum(weights[order])
    index = int(num.searchsorted(cumulative, total_weight / 2.))
    median_deviation = deviations[order[index]]

    if (index + 1 < angles.size
            and num.isclose(cumulative[index], total_weight / 2., rtol=1.e-12, atol=0.)):
        median_deviation = (median_deviation + deviations[order[index + 1]]) / 2.

    summary.update(
        linearity_weighted_mean_angle_deg=float(mean),
        linearity_weighted_median_angle_deg=float((mean + median_deviation + half) % period - half),
        linearity_weighted_stddev_deg=(
            float(num.sqrt(num.sum(weights * deviations**2) / total_weight))
            if angles.size > 1 else None))

    return summary


def orientation_summary(events, period):
    """Calculate angular summaries, with linearity weights for polarization."""
    
    angles = num.asarray([item['angle_deg'] for item in events.values()
                          if item['usable'] and item['angle_deg'] is not None
                          and num.isfinite(item['angle_deg'])], dtype=float)

    summary = dict(median_angle_deg=None, mean_angle_deg=None,
                   scatter_deg=None, n_events=int(angles.size))

    if period == 180.:
        summary.update(_linearity_weighted_summary(events, period))

    if not angles.size:

        return summary

    vector = num.mean(num.exp(2j * num.pi * angles / period))

    if abs(vector) < 1.e-12:
        return summary

    half = period / 2.
    mean = (float(num.angle(vector)) * period / (2. * num.pi) + half) % period - half
    deviations = (angles - mean + half) % period - half
    summary.update(
        median_angle_deg=float((mean + num.median(deviations) + half) % period - half),
        mean_angle_deg=float(mean),
        scatter_deg=float(num.sqrt(num.mean(deviations**2))) if angles.size > 1 else None)

    return summary


def save_orientation_results(directory, filename, event_time, measurements, period):
    """Save an accepted event and station summaries as Guts YAML."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, filename)
    saved = _orientation_mapping(guts.load(filename=path)) if os.path.exists(path) else {}

    for item in saved.values():
        item['events'].pop(event_time, None)

    for nsl, measurement in measurements.items():
        key = '.'.join(nsl)
        item = saved.setdefault(key, dict(summary={}, events={}))
        item['events'][event_time] = dict(measurement)

    for key in sorted(saved):
        events = dict(sorted(saved[key]['events'].items()))
        saved[key] = dict(summary=orientation_summary(events, period), events=events)

    # Replace atomically so an interrupted write does not truncate accepted results.
    descriptor, temporary = tempfile.mkstemp(prefix='.' + filename, dir=directory)
    os.close(descriptor)
    try:
        _rounded_yaml_output(_orientation_document(saved, period)).dump(filename=temporary)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

    logging.info('Saved TeleCheck orientation results: %s', path)


def plot_orientation_results(directory, relative_orient_alignment=False):
    """Plot station rose diagrams and event-time angles from saved YAML."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    methods = [(polarization_filename, 'p_polarization', 'P polarization', 180.)]

    if relative_orient_alignment:
        methods.append((relative_filename, 'relative_alignment',
                        'Relative network alignment', 360.))

    for filename, suffix, label, period in methods:
        path = os.path.join(directory, filename)

        if not os.path.exists(path):
            continue

        saved = _orientation_mapping(guts.load(filename=path))

        for station, item in saved.items():
            events = [(time, measurement['angle_deg'])
                      for time, measurement in sorted(item['events'].items())
                      if measurement['usable'] and measurement['angle_deg'] is not None
                      and num.isfinite(measurement['angle_deg'])]

            if not events:
                continue

            times, values = zip(*events)
            angles = num.asarray(values, dtype=float)
            summary = item['summary']
            mean, median, scatter = [summary[name] for name in
                                    ('mean_angle_deg', 'median_angle_deg', 'scatter_deg')]
            weighted_mean = summary.get('linearity_weighted_mean_angle_deg')
            weighted_median = summary.get('linearity_weighted_median_angle_deg')
            weighted_stddev = summary.get('linearity_weighted_stddev_deg')

            if suffix == 'p_polarization' and weighted_mean is None:
                weighted = _linearity_weighted_summary(item['events'], period)
                weighted_mean = weighted['linearity_weighted_mean_angle_deg']
                weighted_median = weighted['linearity_weighted_median_angle_deg']
                weighted_stddev = weighted['linearity_weighted_stddev_deg']

            caption = 'n=%s, median=%s°, scatter=%s°' % (
                summary['n_events'],
                'undefined' if median is None else '%.1f' % median,
                'undefined' if scatter is None else '%.1f' % scatter)

            fig = Figure(figsize=(7, 6))
            FigureCanvasAgg(fig)
            
            ax = fig.add_subplot(111, projection='polar')
            edges = num.arange(0., period + rose_bin_deg / 2., rose_bin_deg)
            counts, _ = num.histogram(angles % period, bins=edges)
            centres = (edges[:-1] + edges[1:]) / 2.
            # Show both ends of an axial polarization direction, without doubling n.
            offsets = (0., 180.) if period == 180. else (0.,)

            for offset in offsets:
                ax.bar(num.deg2rad(centres + offset), counts,
                       width=num.deg2rad(rose_bin_deg), alpha=.7, edgecolor='black')

            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            radius = ax.get_rmax()

            for value, color, style, legend in (
                    (mean, 'C1', '-', 'Mean'), (median, 'black', '--', 'Median')):

                if value is None:
                    continue

                for i, offset in enumerate(offsets):
                    theta = num.deg2rad(value + offset)
                    ax.plot([theta, theta], [0., radius], color=color,
                            linestyle=style, label=legend if i == 0 else None)

            if mean is not None and scatter is not None:

                for offset in offsets:
                    theta = num.deg2rad(num.linspace(
                        mean + offset - scatter, mean + offset + scatter, 100))
                    ax.fill_between(theta, 0., radius, color='C1', alpha=.15)

            if suffix == 'p_polarization':
                for value, style, legend in (
                        (weighted_mean, '-', 'Linearity-weighted mean'),
                        (weighted_median, '--', 'Linearity-weighted median')):

                    if value is None:
                        continue

                    for i, offset in enumerate(offsets):
                        theta = num.deg2rad(value + offset)
                        ax.plot([theta, theta], [0., radius], color='C2',
                                linestyle=style, linewidth=2., label=legend if i == 0 else None)

                if weighted_mean is not None and weighted_stddev is not None:

                    for i, offset in enumerate(offsets):
                        theta = num.deg2rad(num.linspace(
                            weighted_mean + offset - weighted_stddev,
                            weighted_mean + offset + weighted_stddev, 100))
                        ax.fill_between(theta, 0., radius, color='C2', alpha=.18,
                                        label='Weighted mean ± stddev' if i == 0 else None)

            ax.set_title('%s: %s\n%s' % (station, label, caption), pad=25)

            if (mean is not None or median is not None
                    or weighted_mean is not None or weighted_median is not None):
                ax.legend(loc='lower left', bbox_to_anchor=(0., -.2))

            stem = station.replace('.', '_') + '_' + suffix

            fig.savefig(os.path.join(directory, stem + '_rose.png'),
                        dpi=plot_dpi, bbox_inches='tight')
            fig.clear()

            fig = Figure(figsize=(10, 3.5))
            FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            dates = [datetime.datetime.fromtimestamp(
                util.str_to_time(time), datetime.timezone.utc) for time in times]

            if suffix == 'p_polarization':
                linearities = num.asarray([
                    item['events'][time].get('linearity')
                    if item['events'][time].get('linearity') is not None else num.nan
                    for time in times], dtype=float)
                
                known_linearity = num.isfinite(linearities)
                areas = (linearity_symbol_area_min
                         + (linearity_symbol_area_max - linearity_symbol_area_min)
                         * num.clip(linearities[known_linearity], 0., 1.)**linearity_symbol_exponent)
                ax.scatter(num.asarray(dates)[known_linearity], angles[known_linearity],
                           c='black', s=areas, linewidths=0.)

                for value in linearity_legend_values:
                    area = (linearity_symbol_area_min
                            + (linearity_symbol_area_max - linearity_symbol_area_min)
                            * value**linearity_symbol_exponent)
                    ax.scatter([], [], c='black', s=area, linewidths=0.,
                               label='Linearity %.1f' % value)

                if not known_linearity.all():
                    ax.scatter(num.asarray(dates)[~known_linearity], angles[~known_linearity],
                               facecolors='none', edgecolors='black',
                               s=linearity_symbol_area_min, label='Linearity unavailable')

            else:
                ax.scatter(dates, angles, c='black', s=15)

            if mean is not None:
                ax.axhline(mean, c='black', ls='--', label='Mean')

            if median is not None:
                ax.axhline(median, c='red', ls=':', label='Median')

            if mean is not None and scatter is not None:
                # Wrap the uncertainty band at the limits of the angular interval.
                half = period / 2.

                for offset in (-period, 0., period):
                    low = max(-half, mean - scatter + offset)
                    high = min(half, mean + scatter + offset)

                    if low < high:
                        ax.axhspan(low, high, color='gray', alpha=.2)

            if suffix == 'p_polarization':
                if weighted_mean is not None:
                    ax.axhline(weighted_mean, c='C2', ls='-', linewidth=2.,
                               label='Linearity-weighted mean')

                if weighted_median is not None:
                    ax.axhline(weighted_median, c='C2', ls='--',
                               label='Linearity-weighted median')

                if weighted_mean is not None and weighted_stddev is not None:
                    half = period / 2.
                    band_label = 'Weighted mean ± stddev'

                    for offset in (-period, 0., period):
                        low = max(-half, weighted_mean - weighted_stddev + offset)
                        high = min(half, weighted_mean + weighted_stddev + offset)

                        if low < high:
                            ax.axhspan(low, high, color='C2', alpha=.18, label=band_label)
                            band_label = None

            ax.set_ylim(-period / 2., period / 2.)
            ax.set_xlabel('Event date (UTC)')
            ax.set_ylabel('Angle [°]')
            ax.set_title('%s: %s\n%s' % (station, label, caption))

            if mean is not None or median is not None or suffix == 'p_polarization':
                ax.legend()

            fig.autofmt_xdate()
            fig.savefig(os.path.join(directory, stem + '_overtime.png'),
                        dpi=plot_dpi, bbox_inches='tight')
            fig.clear()


class Correction_Statistics(Object):
    code = String.T()
    median_factor = Float.T()
    mean_factor = Float.T()
    stdev_factor = Float.T()
    n_ev = Int.T()


def get_nslc(st):
    return '%s.%s' % (st.network, st.station)


def get_correction_statistcs(station_list, filename_list, output_directory=None):
    
    mapping_stations_index = {}
    for i_st, st in enumerate(station_list):
        mapping_stations_index[get_nslc(st)] = i_st 

    corrections = num.empty((len(mapping_stations_index), len(filename_list)))
    corrections.fill(num.nan)

    for i_fn, fn in enumerate(filename_list):
        st_cor = grond.load_station_corrections(filename=fn)

        for sc in st_cor:
            if not sc.codes[-1] == 'Z':
                continue
            
            ii = mapping_stations_index['%s.%s' % (sc.codes[0], sc.codes[1])]
            corrections[ii, i_fn] = sc.factor


    #ii_test = mapping_stations_index['CH.FUORN']
    #print(corrections[ii_test, :])

    mean_correction_st = num.nanmean(corrections, axis=1)
    #print(mean_correction_st[ii_test])
    stdev_correction_st = num.nanstd(corrections, axis=1)
    median_correction_st = num.nanmedian(corrections, axis=1)

    is_nan = num.isnan(corrections[:,:]).sum(axis=1)

    corrstats_list = []
    for m, i_m in mapping_stations_index.items():
        logging.info('m %s' % m)
        logging.info('median %s' % median_correction_st[i_m])
        logging.info('mean %s' % mean_correction_st[i_m])
        logging.info('std %s' % stdev_correction_st[i_m])
        logging.info('len %s' % (len(filename_list)-is_nan[i_m]))

        corrstats_list.append(Correction_Statistics(code=m,
                                                    median_factor=median_correction_st[i_m],
                                                    mean_factor=mean_correction_st[i_m],
                                                    stdev_factor=stdev_correction_st[i_m],
                                                    n_ev=len(filename_list)-is_nan[i_m]))
    

    dump_all(_rounded_yaml_output(corrstats_list),
             filename=os.path.join(output_directory or 'results/tele_check',
                                   'tele_check_stats.yaml'),
             stream=None)
