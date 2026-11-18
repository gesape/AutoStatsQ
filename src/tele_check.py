# https://github.com/gesape/AutoStatsQ - GPLv3
#
# AutoStasQ is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
#
# The AutoStatsQ Developers, 21st Century
# ---|P------/S----------~Lg----------

'''
File containing tele check functions.
'''

import os, io, logging
import numpy as num

from pyrocko.gui_util import PhaseMarker
from pyrocko.snuffling import Snuffling, Choice, Param
from pyrocko.gui import marker as pmarker
from pyrocko import cake, trace, util
from pyrocko.plot import mpl_color

from . import call_tele_check as tele


d2r = num.pi / 180.
r2d = 1.0 / d2r
relative_alignment_color_threshold_deg = 20.0
polarization_axis_length = 0.9


class NoArrival(Exception):
    pass


class NotFound(Exception):
    pass


def get_trace(trs, func):
    trs_match = [tr for tr in trs if func(tr)]
    if len(trs_match) == 1:
        return trs_match[0]

    raise NotFound()


def invert_relative_amplitudes(pair_corrs):
    nslcs = sorted(set(x[0] for x in pair_corrs))
    nslc_to_i = dict((nslc, i) for (i, nslc) in enumerate(sorted(nslcs)))

    ncols = len(nslcs)

    rows = []
    data = []
    weights = []
    for a_nslc, b_nslc, _, relamp, _, _ in pair_corrs:
        a_i = nslc_to_i[a_nslc]
        b_i = nslc_to_i[b_nslc]
        row = num.zeros(ncols)
        row[a_i] = 1.
        row[b_i] = -1.
        rows.append(row)
        data.append(num.log(abs(relamp)))
        weights.append(1.0)

    for nslc in nslcs:
        i = nslc_to_i[nslc]
        row = num.zeros(ncols)
        row[i] = 1.0
        rows.append(row)
        data.append(0.0)
        weights.append(0.01)

    data = num.array(data, dtype=float)
    mat = num.vstack(rows)
    weights = num.array(weights)

    mod, res = num.linalg.lstsq(
        mat * weights[:, num.newaxis], data*weights,rcond=None)[0:2]

    relamps_rel = num.exp(mod)

    relamps_rel /= num.median(relamps_rel)

    corrs = {}
    for nslc in nslcs:
        i = nslc_to_i[nslc]
        factor = relamps_rel[i]
        #print(nslc, factor)
        corrs[nslc] = factor

    return corrs


class TeleCheck(Snuffling):

    def setup(self):
        self.set_name('Tele-Check')
        self.add_parameter(
            Choice(
                'Phase',
                'phasename',
                'P', ['P', 'S', 'PP', 'PPP', 'SS']))

        self.add_parameter(
            Choice(
                'Channels for polarization analysis',
                'channels_polar',
                'R, T', ['R, T', 'N, E']))

        self.add_parameter(
            Choice(
                'Channels for relative amplitude analysis',
                'channels_relamp',
                'Z', ['Z', 'R', 'T', 'E', 'N', 'All']))

        self.add_parameter(
            Param(
                'Minimum cross-correlation value',
                'cc_min',
                0.8, 0.1, 0.95))

        self.add_trigger(
            'Save Grond Corrections',
            self.save_grond_corrections)

        self.fframes = {}
        self.set_live_update(False)
        self._nslc_to_relamp = None
        self.relative_orient_alignment = False
        self.output_directory = os.path.join('results', 'tele_check')
        self._pending_orientation = None
        self._polar_display_filter = None
        self._previous_quick_filter = None
        self._polar_display_channels = None

    def update_channel_display(self, restore=False):
        """Filter displayed channels without changing the waveform pile."""
        viewer = self.get_viewer()

        if self.channels_polar in ('R, T', 'N, E') and not restore:

            if (self._polar_display_filter is not None
                    and viewer.quick_filter is self._polar_display_filter
                    and self._polar_display_channels == self.channels_polar):

                return

            # Switching channel selections must not chain incompatible filters.
            previous = (self._previous_quick_filter
                        if self._polar_display_filter is not None
                        and viewer.quick_filter is self._polar_display_filter
                        else viewer.quick_filter)
            self._previous_quick_filter = previous

            channels = self.channels_polar
            self._polar_display_filter = lambda tr: (
                (tr.channel.endswith(('Z', 'R', 'T')) if channels == 'R, T'
                 else not tr.channel.endswith(('R', 'T')))
                and (previous is None or previous(tr)))
            self._polar_display_channels = channels
            viewer.set_quick_filter(self._polar_display_filter)
        else:

            if (self._polar_display_filter is not None
                    and viewer.quick_filter is self._polar_display_filter):
                viewer.set_quick_filter(self._previous_quick_filter)

            self._polar_display_filter = None
            self._previous_quick_filter = None
            self._polar_display_channels = None

        if not viewer.closing:
            viewer.update()

    def choose_on_snuffling_panel(self, ident, state):
        # """Update the channel display immediately, without running analysis."""
        """Refresh traces when polarization channels change."""
        super().choose_on_snuffling_panel(ident, state)

        if ident == 'channels_polar':
            self.update_channel_display()

    def pre_destroy(self):
        # """Restore the prior display filter when TeleCheck is removed."""
        """Restore the display filter and release resources."""

        if self._viewer is not None:
            self.update_channel_display(restore=True)

        super().pre_destroy()

    def _capture_polarization_view(self, pending, axes):
        """Keep a PNG snapshot of the accepted polarization panel."""
        pending['polarization_view'] = self._capture_panel_view(axes)

    def _capture_panel_view(self, axes):
        """Capture a result panel as PNG bytes."""
        figure = axes.figure
        figure.canvas.draw()
        bounds = axes.get_tightbbox(figure.canvas.get_renderer())
        bounds = bounds.transformed(figure.dpi_scale_trans.inverted()).expanded(1.05, 1.10)

        with io.BytesIO() as buffer:
            figure.savefig(buffer, format='png', dpi=tele.plot_dpi, bbox_inches=bounds)

            return buffer.getvalue()


    def save_grond_corrections(self):
        """Save accepted orientation results and gain corrections."""
        # Save only a complete Run associated with the currently active event.
        viewer = self.get_viewer()
        event = viewer.get_active_event()
        pending = self._pending_orientation

        if (event is None or pending is None
                or pending['event_time'] != util.time_to_str(event.time)):
            self.error('Run TeleCheck for the active event before saving.')

            return

        os.makedirs(self.output_directory, exist_ok=True)

        view_path = os.path.join(
            self.output_directory, '%s_polarization_view.png' % event.name)

        with open(view_path, 'wb') as stream:
            stream.write(pending['polarization_view'])

        logging.info('Saved TeleCheck polarization view: %s', view_path)

        tele.save_gain_results(
            self.output_directory, pending['event_time'], pending['gain_stations'],
            pending['gain_factors'], pending['gain_components'])
        gain_view_path = os.path.join(
            self.output_directory, '%s_gain_view.png' % event.name)

        with open(gain_view_path, 'wb') as stream:
            stream.write(pending['gain_view'])

        logging.info('Saved TeleCheck gain view: %s', gain_view_path)

        tele.save_orientation_results(
            self.output_directory, tele.polarization_filename,
            pending['event_time'], pending['polarization'], 180.)

        if self.relative_orient_alignment:
            tele.save_orientation_results(
                self.output_directory, tele.relative_filename,
                pending['event_time'], pending['relative'], 360.)

        if self._nslc_to_relamp:
            self.cleanup()
            viewer = self.get_viewer()
            event = viewer.get_active_event()
            fn = os.path.join(self.output_directory, '%s.cor' % event.name)

            import grond

            corrs = []

            for nslc in sorted(self._nslc_to_relamp.keys()):
                for c in 'ZRT':
                    corr = grond.StationCorrection(
                        codes=nslc[:3] + (c,),
                        delay=0.0,
                        factor=float(self._nslc_to_relamp[nslc]))
                    corrs.append(corr)

            grond.dump_station_corrections(corrs, fn)


    def call(self):
        # Invalidate a previous run before any operation that could fail.
        self._pending_orientation = None
        self._nslc_to_relamp = None
        self.cleanup()
        viewer = self.get_viewer()
        event = viewer.get_active_event()
        stations = self.get_stations()

        #for s in stations:
        #    print(s.nsl())

        nsl_to_station = dict(
            (s.nsl(), s) for s in stations)

        if event is None:
            self.error('No active event set.')

        markers = self.get_selected_markers()
        if len(markers) != 1:
            self.error('Exactly one marker must be selected.')

        marker = markers[0]

        try:
            nslc = marker.one_nslc()
        except pmarker.MarkerOneNSLCRequired:
            self.error('Marker must be picked on a single trace.')

        marker_station = nsl_to_station[nslc[:3]]

        mod = cake.load_model()

        def traveltime(station):
            dist = event.distance_to(station)
            arrivals = mod.arrivals(
                zstart=event.depth,
                zstop=0.,
                distances=[dist*cake.m2d],
                phases=[cake.PhaseDef(self.phasename)])

            if not arrivals:
                raise NoArrival()

            return arrivals[0].t

        try:
            tt_marker_station = traveltime(marker_station)
        except NoArrival:
            self.error('Selected phase does not arrive at station.')

        nsl_to_delay = {}
        for station in stations:
            nsl_to_delay[station.nsl()] = \
                traveltime(station) - tt_marker_station

        pile = self.get_pile()

        nsl_to_traces = {}
        nsl_to_tspan = {}
        fs = [f for f in [viewer.lowpass, viewer.highpass] if f is not None]
        if fs:
            tpad = 2.0 / min(fs)

        else:
            tpad = 0.0

        deltats = set()
        for nsl in nsl_to_delay.keys():
            delay = nsl_to_delay[nsl]
            tmin = marker.tmin + delay
            tmax = marker.tmax + delay
            nsl_to_tspan[nsl] = (tmin, tmax)

            trs = pile.all(
                tmin=tmin, tmax=tmax, tpad=tpad,
                trace_selector=lambda tr: tr.nslc_id[:3] == nsl,
                want_incomplete=False)

            for tr in trs:
                if viewer.lowpass is not None:
                    tr.lowpass(4, viewer.lowpass)

                if viewer.highpass is not None:
                    tr.highpass(4, viewer.highpass)

                tr.chop(tr.wmin, tr.wmax)
                deltats.add(tr.deltat)

            if trs:
                nsl_to_traces[nsl] = trs

        if len(deltats) != 1:
            self.error('All traces must have same sampling rate.')

        # add markers
        for nsl in nsl_to_traces.keys():
            tmin, tmax = nsl_to_tspan[nsl]

            for tr in nsl_to_traces[nsl]:
                mark = PhaseMarker(
                        [tr.nslc_id],
                        tmin=tmin,
                        tmax=tmax,
                        kind=1,
                        phasename=self.phasename)

                self.add_marker(mark)

        # cross correlations
        nsls = sorted(list(nsl_to_traces.keys()))
        pair_corrs = []
        for nsl_a in nsls:
            trs_a = nsl_to_traces[nsl_a]

            if self.channels_relamp == 'All':
                comps = sorted(set([tr.channel[-1] for tr in trs_a]))
            else:
                comps = [c.strip() for c in self.channels_relamp.split(',')]

            for nsl_b in nsls:
                trs_b = nsl_to_traces[nsl_b]

                for comp in comps:
                    try:
                        tr_a = get_trace(
                            trs_a, lambda tr: tr.channel.endswith(comp))
                        tr_b = get_trace(
                            trs_b, lambda tr: tr.channel.endswith(comp))

                    except NotFound:
                        continue

                    if tr_a is tr_b:
                        continue

                    tr_cor = trace.correlate(
                        tr_a, tr_b, mode='full', normalization='normal')

                    delaymax, ccmax = tr_cor.max()
                    delaymin, ccmin = tr_cor.min()

                    delay_syn = nsl_to_delay[nsl_b] - nsl_to_delay[nsl_a]

                    if abs(ccmin) < abs(ccmax):
                        delay = delaymax
                        ccabsmax = abs(ccmax)
                        ccsignedmax = ccmax
                    else:
                        delay = delaymin
                        ccabsmax = abs(ccmin)
                        ccsignedmax = ccmin

                    tr_b_shifted = tr_b.copy()
                    tr_b_shifted.shift(-delay)

                    tmin_com = max(tr_b_shifted.tmin, tr_a.tmin)
                    tmax_com = min(tr_b_shifted.tmax, tr_a.tmax)

                    tr_a_chopped = tr_a.chop(
                        tmin_com, tmax_com, inplace=False)
                    tr_b_chopped = tr_b_shifted.chop(
                        tmin_com, tmax_com, inplace=False)

                    ya = tr_a_chopped.ydata
                    yb = tr_b_chopped.ydata

                    relamp1 = num.sum(ya*yb) / num.sum(yb**2)
                    relamp2 = num.sum(ya*yb) / num.sum(ya**2)

                    if nsl_a[1] == 'LYKK':
                        print(
                            ccabsmax, relamp1, relamp2,
                            abs((relamp1 / (1.0/relamp2) - 1.0)))

                    if ccabsmax < self.cc_min:
                        continue

                    if abs((relamp1 / (1.0/relamp2) - 1.0)) > 0.2:
                        continue

                    relamp = (relamp1+1./relamp2) * 0.5

                    pair_corrs.append(
                        (tr_a.nslc_id, tr_b.nslc_id,
                         ccsignedmax, relamp, delay, delay_syn))

        # Orientation remains available when no gain pairs pass the quality checks.
        nslc_to_relamp = invert_relative_amplitudes(pair_corrs) if pair_corrs else {}
        self._nslc_to_relamp = nslc_to_relamp

        nsl_to_xy = {}
        for nsl in nsl_to_traces.keys():
            trs = nsl_to_traces[nsl]
            try:
                cc = [c.strip() for c in self.channels_polar.split(',')]

                tr_y, tr_x = [
                    get_trace(trs, lambda tr: tr.channel.endswith(c))
                    for c in cc]

                x = tr_x.get_ydata()
                y = tr_y.get_ydata()

                nsl_to_xy[nsl] = (x, y)

            except NotFound:
                pass

        nsls = sorted(list(nsl_to_xy.keys()))
        n = len(nsls)

        # Original combined polarization and relative-alignment implementation:
        # xs_l = [nsl_to_xy[nsl][0] for nsl in nsls]
        # ys_l = [nsl_to_xy[nsl][1] for nsl in nsls]
        # nsamp = min(min(x.size for x in xs_l), min(y.size for y in ys_l))
        # xs = num.vstack([x[:nsamp] for x in xs_l])
        # ys = num.vstack([y[:nsamp] for y in ys_l])
        # amps = num.sqrt(xs**2 + ys**2)
        # amp_maxs = num.max(amps, axis=1)
        # xs = xs / amp_maxs[:, num.newaxis]
        # ys = ys / amp_maxs[:, num.newaxis]
        # nphi = 361 ##
        # phis = num.linspace(-180., 180., nphi)
        # d = num.zeros((n, n, nphi))
        # print('Result Event %s' % str(util.tts(event.time)))
        # for ia in range(n):
        #     for iphi, phi in enumerate(phis):
        #         x = xs[ia, :]
        #         y = ys[ia, :]
        #         xrot = num.cos(phi*d2r) * x + num.sin(phi*d2r) * y
        #         yrot = - num.sin(phi*d2r) * x + num.cos(phi*d2r) * y
        #         d[ia, :, iphi] = num.sqrt(num.sum(
        #                 (xrot[num.newaxis, :] - xs)**2
        #                 + (yrot[num.newaxis, :] - ys)**2, axis=1))
        #     data = num.zeros((2,len(x)))
        #     data[0] = x
        #     data[1] = y
        #     cov = num.cov(data)
        #     evals, evecs = num.linalg.eigh(cov)
        #     pc = evecs[:, -1]
        #     azimuth = r2d*num.arctan2(pc[1], pc[0])
        #     azimuth = (-azimuth) % 180. - 90.
        #     print('POL AZI (assumes R,T): %-20s %8.0f' % ('.'.join(nsls[ia]), azimuth))
        # imins = num.argmin(d, axis=2)
        # dmins = num.min(d, axis=2)
        # dmin_median = num.median(dmins)
        # phimins = phis[imins]
        # nsl_to_rot = {}
        # for nsl in nsls:
        #     nsl_to_rot[nsl] = 0.
        # failed = set()
        # for i in range(n):
        #     mean_min_error = num.mean(dmins[i, :]/dmin_median)
        #     #print(mean_min_error, nsls[i])
        #     if mean_min_error > 3.0:
        #         failed.add(nsls[i])
        # while True:
        #     ia_worst = num.argmax(num.mean(num.abs(phimins), axis=1))
        #     phimod = (
        #         (phis[num.newaxis, :] + phimins[ia_worst, :, num.newaxis]
        #          + 180.) % 360.) - 180.
        #     phirot = phis[num.argmin(num.mean(num.abs(phimod), axis=0))]
        #     if abs(phirot) < 1.:
        #         print('PHI ROT <1: %-20s %8.0f' % ('.'.join(nsl), phirot))
        #         break
        #     nsl = nsls[ia_worst]
        #     mean_min_error = num.mean(dmins[ia_worst, :]/dmin_median)
        #     phimins[ia_worst, :] = (
        #         (phimins[ia_worst, :] + phirot) + 180.) % 360. - 180.
        #     phimins[:, ia_worst] = (
        #         (phimins[:, ia_worst] - phirot) + 180.) % 360. - 180.
        #     if nsl not in failed:
        #         print('PHI ROT: %-20s %8.0f' % ('.'.join(nsl), phirot))
        #         nsl_to_rot[nsl] += phirot

        polarization = {
            nsl: dict(angle_deg=None, linearity=None, usable=False,
                      status='missing_channels') for nsl in nsl_to_station}
        relative = {}

        if nsls:
            nsamp = min(min(x.size, y.size) for x, y in nsl_to_xy.values())
            xs = num.vstack([nsl_to_xy[nsl][0][:nsamp] for nsl in nsls])
            ys = num.vstack([nsl_to_xy[nsl][1][:nsamp] for nsl in nsls])

            for ia, nsl in enumerate(nsls):
                polarization[nsl] = tele.polarization_measurement(xs[ia], ys[ia])

            # No relative-alignment work or allocation when disabled.
            if self.relative_orient_alignment and nsamp > 0:
                relative = tele.relative_orientation(nsls, xs, ys)

        if self.relative_orient_alignment:
            for nsl in nsl_to_station:
                relative.setdefault(nsl, dict(
                    angle_deg=None, relative_waveform_misfit=None,
                    usable=False, status='missing_channels'))

        for nsl, result in polarization.items():
            result.update(phase=self.phasename, channels=self.channels_polar)
            # Only P on R,T supports interpretation as a radial P-axis deviation.
            if self.phasename != 'P' or self.channels_polar != 'R, T':
                result.update(usable=False, status='requires_P_on_RT')

        # for measurements in (polarization, relative):
        #     for nsl, result in measurements.items():
        #         result.update(
        #             window_start=util.time_to_str(nsl_to_tspan[nsl][0]),
        #             window_end=util.time_to_str(nsl_to_tspan[nsl][1]),
        #             highpass_hz=viewer.highpass, lowpass_hz=viewer.lowpass)
        for method, measurements in (
                ('P polarization', polarization), ('Relative network alignment', relative)):
            for nsl, result in measurements.items():
                logging.debug(
                    'TeleCheck %s: station=%s event=%s phase=%s channels=%s '
                    'window_start=%s window_end=%s status=%s usable=%s '
                    'highpass_hz=%s lowpass_hz=%s',
                    method, '.'.join(nsl), util.time_to_str(event.time),
                    self.phasename, self.channels_polar,
                    util.time_to_str(nsl_to_tspan[nsl][0]),
                    util.time_to_str(nsl_to_tspan[nsl][1]),
                    result['status'], result['usable'], viewer.highpass, viewer.lowpass)

        pending = dict(event_time=util.time_to_str(event.time),
                       polarization=polarization, relative=relative)
        gain_components = (
            sorted({tr.channel[-1] for traces in nsl_to_traces.values() for tr in traces})
            if self.channels_relamp == 'All'
            else [component.strip() for component in self.channels_relamp.split(',')])
        pending.update(gain_stations=list(nsl_to_station),
                       gain_factors=dict(nslc_to_relamp), gain_components=gain_components)

        fframe = self.figure_frame()
        fig = fframe.gcf()
        fig.clf()


        ncols = 1
        while ncols**2 < n:
            ncols += 1

        nrows = ncols

        axes = fig.add_subplot(1, 2, 1, aspect=1.0)
        axes.axison = False
        axes.set_xlim(-0.05 - ncols, ncols + 0.05)
        axes.set_ylim(-0.05 - nrows, nrows + 0.05)

        axes.set_title('P polarization: %s (black line)%s' % (
            event.name, '\nRelative network alignment enabled (gray/orange line)'
            if self.relative_orient_alignment else ''))

        for insl, nsl in enumerate(nsls):
            irow = insl // ncols
            icol = insl % ncols
            trs = nsl_to_traces[nsl]
            try:
                x, y = nsl_to_xy[nsl]
                cc = [c.strip() for c in self.channels_polar.split(',')]

                tr_y, tr_x = [
                    get_trace(trs, lambda tr: tr.channel.endswith(c))
                    for c in cc]

                xpos = icol*2-ncols+1
                ypos = -irow*2+nrows-1
                x = tr_x.get_ydata()
                y = tr_y.get_ydata()

                a = num.sqrt(x**2 + y**2)
                amax = num.max(a)

                if not num.isfinite(amax) or amax <= 0.:
                    continue

                axes.plot(x / amax + xpos, y / amax + ypos,
                          color='black', alpha=0.7)
                result = polarization[nsl]
                angle = result['angle_deg']

                annotation = '.'.join(part for part in nsl if part)

                if angle is not None and num.isfinite(angle):
                    # The angle is measured from the first component (R or N).
                    direction_x = num.sin(angle * d2r)
                    direction_y = num.cos(angle * d2r)
                    ends = num.asarray([-polarization_axis_length, polarization_axis_length])
                    axes.plot(xpos + ends * direction_x, ypos + ends * direction_y,
                              color='black', linestyle=':', linewidth=1.5, zorder=4)

                    if self.relative_orient_alignment and relative[nsl]['usable']:
                        phi = relative[nsl]['angle_deg']
                        aligned_x = num.cos(phi*d2r) * direction_x - num.sin(phi*d2r) * direction_y
                        aligned_y = num.sin(phi*d2r) * direction_x + num.cos(phi*d2r) * direction_y
                        alignment_color = (
                            'gray' if abs(phi) < relative_alignment_color_threshold_deg
                            else mpl_color('orange2'))
                        axes.plot(xpos + ends * aligned_x, ypos + ends * aligned_y,
                                  color=alignment_color, linestyle='--', linewidth=1.5, zorder=5)

                axes.annotate(
                    annotation, xy=(xpos, ypos - 1.), xycoords='data',
                    xytext=(0, 0), textcoords='offset points',
                    verticalalignment='center', horizontalalignment='center',
                    fontsize=8, rotation=0.,
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.75))

            except NotFound:
                pass

        polarization_axes = axes
        axes = fig.add_subplot(1, 2, 2)
        axes.set_title('Relative gain: %s (%s)' % (event.name, self.channels_relamp))
        nslcs = sorted(nslc_to_relamp.keys())
        pdata = []
        for inslc, nslc in enumerate(nslcs):
            nsl = nslc[:3]
            cha = nslc[3]
            tr = get_trace(
                nsl_to_traces[nsl],
                lambda tr: tr.channel == cha).copy()

            tr.shift(
                -(event.time + tt_marker_station + nsl_to_delay[nsl]))

            relamp = nslc_to_relamp[nslc]
            tr.ydata /= relamp

            color = 'black'
            if abs(num.log10(relamp)) > num.log10(1.1):
                color = mpl_color('chocolate2')

            if abs(num.log10(relamp)) > num.log10(2.0):
                color = mpl_color('orange2')

            if abs(num.log10(relamp)) > num.log10(10.0):
                color = mpl_color('scarletred2')

            pdata.append((tr, relamp, color, inslc, nslc))

        if not pdata:
            axes.text(.5, .5, 'No usable gain pairs', ha='center',
                      transform=axes.transAxes)
            fframe.draw()
            self._capture_polarization_view(pending, polarization_axes)
            pending['gain_view'] = self._capture_panel_view(axes)
            self._pending_orientation = pending

            return

        ranges = trace.minmax([aa[0] for aa in pdata], lambda tr: None)
        ymin, ymax = ranges[None]
        yabsmax = max(abs(ymin), abs(ymax))

        for (tr, relamp, color, inslc, nslc) in pdata:
            axes.plot(
                tr.get_xdata(),
                inslc + tr.get_ydata()/yabsmax,
                color=color)

            axes.annotate(
                '.'.join(_ for _ in nslc if _),
                xy=(0, inslc),
                xycoords=('axes fraction', 'data'),
                xytext=(-5, 0),
                textcoords='offset points',
                verticalalignment='center',
                horizontalalignment='right',
                rotation=0.,
                color=color)

            axes.annotate(
                'x %g' % (1.0/relamp),
                xy=(1., inslc),
                xycoords=('axes fraction', 'data'),
                xytext=(+5, 0),
                textcoords='offset points',
                verticalalignment='center',
                horizontalalignment='left',
                rotation=0.,
                color=color)

        axes.get_yaxis().set_visible(False)
        for which in ['top', 'right', 'left']:
            axes.spines[which].set_visible(False)

        axes.set_xlabel('Time [s]')
        fframe.draw()
        self._capture_polarization_view(pending, polarization_axes)
        pending['gain_view'] = self._capture_panel_view(axes)
        self._pending_orientation = pending


def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''

    return [TeleCheck()]
