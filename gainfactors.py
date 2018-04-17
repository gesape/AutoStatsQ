#!/use/bin/env python

import os
import logging
import gc
import numpy as num
from collections import defaultdict
from pyrocko import cake, util
from pyrocko.orthodrome import distance_accurate50m
from pyrocko.guts import Object, Dict, String, Float

import gainplots as gp

# based on https://github.com/HerrMuellerluedenscheid/autogain

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
km = 1000.


class Gains(Object):
    trace_gains = Dict.T(String.T(), Float.T())
    ref_stats = String.T(optional=True)


class PhasePie():
    '''Calculates and caches phase arrivals'''

    def __init__(self, mod='prem-no-ocean.f'):
        '''
        :param mod: Name of the model to be used.
        '''
        self.model = cake.load_model(mod)
        self.arrivals = defaultdict(dict)
        self.which = None

    def t(self, phase_selection, z_dist):
        '''
        :param phase_selection: phase names speparated by vertical bars
        :param z_dist: tuple with (depth, distance)
        '''
        if 'first' in phase_selection:
            self.which = 'first'

        if 'last' in phase_selection:
            self.which = 'last'

        if self.which:
            phase_selection = self.strip(phase_selection)

        z, dist = z_dist

        if (phase_selection, dist, z) in self.arrivals.keys():
            return self.arrivals[(phase_selection, dist, z)]

        phases = [cake.PhaseDef(pid) for pid in phase_selection.split('|')]

        arrivals = self.model.arrivals(distances=[dist*cake.m2d],
                                       phases=phases,
                                       zstart=z)

        if arrivals == []:
            logger.info('none of defined phases at d=%s, z=%s' % (dist, z))
            return

        else:
            want = self.phase_selector(arrivals)
            self.arrivals[(phase_selection, dist, z)] = want
            return want

    def phase_selector(self, _list):
        if self.which == 'first':
            return min(_list, key=lambda x: x.t).t
        if self.which == 'last':
            return max(_list, key=lambda x: x.t).t

    def strip(self, ps):
        ps = ps.replace(self.which, '')
        ps = ps.rstrip(')')
        ps = ps.lstrip('(')
        return ps


class StaticLengthWindow():
    def __init__(self, static_length, phase_position):
        '''
        phase_position: 0-> start ... 0.5 -> center ... 1.0 -> end'''
        self.phase_position = phase_position
        self.static_length = static_length

    def t(self):
        return self.static_length*self.phase_position,\
               self.static_length*1.0-self.phase_position


def guess_nsl_template(code):
    if len(code) == 1 or isinstance(code, str):
        return '*.%s.*.*' % (code)
    elif len(code) == 2:
        return '*.%s.%s.*' % (code)
    elif len(code) == 3:
        return '%s.%s.%s.*' % (code)


class Section():
    ''' Related to one event.
    All traces scale relative to average median abs
    max, to one reference station or to 1.0'''
    def __init__(self, event, stations):
        self.stations = stations
        self.event = event
        self.reference_scale = None
        self.max_tr = {}
        self.relative_scalings = {}
        self.finished = False

    def finish(self, method, fband, taper, ev_counter):
        if len(method) == 2:
            reference_nsl = method[1][1]
            # print(reference_nsl)
            reference_nslc = list(filter(
                lambda x: util.match_nslc(guess_nsl_template(reference_nsl), x),
                                          self.max_tr.keys()))
            self.____reference_nslc = reference_nslc

            if not len(reference_nslc) == 1:
                logger.info('no reference trace available. ' +
                            'remains unfinished: %s' % self.event)
                self.finished = False
            else:
                self.reference_scale = self.max_tr[reference_nslc[0]]
                self.set_relative_scalings()
                self.finished = True

        elif method == 'scale_one' or method == 'median_all_avail':
            self.reference_scale = 1.
            self.set_relative_scalings()
            self.finished = True

    def set_relative_scalings(self):
        for nslc_id, maxs in self.max_tr.items():
            self.relative_scalings[nslc_id] = self.reference_scale/maxs

    def get_gained_traces(self):
        gained = []
        for tr in self.traces:
            tr = tr.copy()
            tr.ydata *= self.relative_scalings[tr.nslc_id]
            tr.set_location('G')
            gained.append(tr)
        return gained

    def get_ungained_traces(self):
        return self.traces

    def iter_scalings(self):
        for nslc_id, scaling in self.relative_scalings.items():
            yield (nslc_id, scaling)


class AutoGain():
    def __init__(self, data_pile, stations, events,
                 component='Z', gain_rel_to='scale_one',
                 phase_selection='first(p|P)'):
        '''
       :param phase_selection: follows the logic of
                               fomosto's Store phase definitions
       :param gain_rel_to: gain relative to options:
                           * 'scale_one' for reference amplitude = 1.
                           *  tuple ('reference_nsl', refernce_id) - gain
                              relative to one specific station,
                              reference_id e.g. 'BFO'
                           * 'median_all_avail': first assesses all abs.
                              amplitudes, in the end chooses a reference
                              station of those that recorded
                              most events based on median gain of first
                              event
       :param events: List of pyrocko events
       :param stations: List of pyrocko stations
       :param phase_selection: Phases for arrival time calculation, pyrocko
                               cake naming. Needs to be adjusted for 
                               teleseismic events.
       :param data_pile: Pyrocko data pile 

        '''
        self.method = gain_rel_to
        self.component = component
        self.data_pile = data_pile
        self.stations = stations
        self.events = events
        self.phaser = None
        self.phase_selection = phase_selection
        self.all_nslc_ids = set()
        self.minmax = {}
        self.scaling_factors = {}
        self.sections = []
        self.results = None
        self._mean = None

    def process(self, fband, taper):
        no_events = len(self.events)
        for i_ev, event in enumerate(self.events):
            tr_nslc_ids = []
            print('Processing event %s of %s' % (i_ev, no_events))
            section = Section(event, self.stations)
            skipped = 0
            unskipped = 0
            for i_s, s in enumerate(self.stations):
                dist = distance_accurate50m(event, s)
                arrival = self.phaser.t(self.phase_selection,
                                        (event.depth, dist))
                if arrival is None:
                    skipped += 1
                    continue
                else:
                    unskipped += 1

                selector = lambda tr: util.match_nslc('%s.*%s' %
                                                      (s.nsl_string(),
                                                       self.component),
                                                      tr.nslc_id)

                window_min, window_max = self.window.t()
                tr_generator = self.data_pile.chopper(tmin=event.time+arrival -
                                                      window_min,
                                                      tmax=event.time+arrival +
                                                      window_max,
                                                      trace_selector=selector,
                                                      load_data=True)

                for tr in tr_generator:
                    if not len(tr) > 1 and tr:
                        tr = tr[0]
                        if len(tr.ydata) > 0 and num.max(num.abs(tr.get_ydata())) != 0:
                            dtype = type(tr.ydata[0])
                            tr.ydata -= dtype(tr.get_ydata().mean())
                            tr.highpass(fband['order'], fband['corner_hp'])
                            tr.taper(taper, chop=False)
                            tr.lowpass(fband['order'], fband['corner_lp'])
                            if num.max(num.abs(tr.get_ydata())) != 0:
                                section.max_tr[tr.nslc_id] = num.max(num.abs(tr.get_ydata()))
                                tr_nslc_ids.append(tr.nslc_id)

            logger.debug('skipped %s/%s' % (skipped, unskipped))

            section.finish(self.method, fband, taper, i_ev)

            self.all_nslc_ids.update(tr_nslc_ids)
            gc.collect()

            self.sections.append(section)

            if self.method == 'median_all_avail' and i_ev == no_events-1:
                self.handle_median_stats_option()

    def congregate(self):
        indx = dict(zip(self.all_nslc_ids, num.arange(len(self.all_nslc_ids))))
        self.results = num.empty((len(self.sections), len(self.all_nslc_ids)))
        self.results[:] = num.nan
        for i_sec, section in enumerate(self.sections):
            for nslc_id, scaling in section.iter_scalings():
                self.results[i_sec, indx[nslc_id]] = scaling

    def set_phaser(self, phaser):
        self.phaser = phaser

    def set_window(self, window):
        self.window = window

    def handle_median_stats_option(self):
        indx = dict(zip(self.all_nslc_ids, num.arange(len(self.all_nslc_ids))))
        results_all = num.empty((len(self.sections), len(self.all_nslc_ids)))
        results_all[:] = num.nan

        for i_ev, section in enumerate(self.sections):  # loop over events
            for nslc_id, scaling in section.iter_scalings():
                try:
                    results_all[i_ev, indx[nslc_id]] = scaling
                except:
                    continue

        stats_list = list(self.all_nslc_ids)

        lst_st_avail = []
        no_stats = results_all.shape[1]
        no_ev = results_all.shape[0]

        for i_st in range(no_stats):
            no_avail = 0
            for i_ev in range(no_ev):
                if not num.isnan(results_all[i_ev, i_st]):
                    no_avail += 1
            lst_st_avail.append(no_avail)

        ind_most_events = num.argwhere(lst_st_avail == num.amax(lst_st_avail))
        # indices of all stations that recorded max. number of events

        candidates_ref_values = [results_all[0, i_st]
                                 for i_st in ind_most_events]
        # gain values according to those indices

        i_median_stat_rel = num.argsort(candidates_ref_values)\
         [len(candidates_ref_values)//2]
        # indice of median gain station, with respect to list ind_most_events

        i_median_stat = ind_most_events[int(i_median_stat_rel)]

        median_stat = stats_list[i_median_stat[0]]

        results_rel = num.empty((results_all.shape))
        results_rel[:] = num.nan

        for i_row, row in enumerate(results_all):
            for i_element, element in enumerate(row):
                results_rel[i_row, i_element] = element / row[i_median_stat[0]]

        self.method = ('reference_nsl_med', median_stat)
        self.results = results_rel

    @property
    def mean(self):
        if self.results is None:
            self.congregate()
        if self._mean is None and self.method != 'syn':
            self._mean = dict(zip(map(lambda x: '.'.join(x),
                                      self.all_nslc_ids),
                                  num.nanmean(self.results, axis=0)))
            # might rise a warning if for one station only nan are in results,
            # but works fine.
        if self._mean is None and self.method == 'syn':
            nslcs = [st.split(' ')[1]+'.'+st.split(' ')[2]+'..'+'HHZ'
                     for st in self.stations]
            self._mean = dict(zip(nslcs,
                                  num.nanmean(self.results, axis=0)))
        return self._mean

    def save_mean(self, fn, directory):
        g = Gains()
        g.trace_gains = self.mean
        if len(self.method) == 2:
            g.ref_stats = '%s %s' % (self.method[1][0], self.method[1][1])
        g.regularize()
        g.validate()
        g.dump(filename=directory+fn)

    def save_single_events(self, fn, directory, plot=False):
        ''' Save table with all gains (for all events + stations)
        '''
        if not self.method[0] == 'reference_nsl_med' and not self.method == 'syn':
            indx = dict(zip(self.all_nslc_ids,
                            num.arange(len(self.all_nslc_ids))))
            results_all = num.empty((len(self.sections),
                                     len(self.all_nslc_ids)))
            results_all[:] = num.nan
            for i_ev, section in enumerate(self.sections):  # loop over events
                for nslc_id, scaling in section.iter_scalings():
                    results_all[i_ev, indx[nslc_id]] = scaling

        else:
            results_all = self.results

        stats_list = self.all_nslc_ids
        if self.method == 'syn':
            stats_list = self.stations

        with open(directory+fn, 'w') as outfile:
            outfile.write('Gain relative to station %s. \n'
                          % (str(self.method[1])))
            outfile.write('Station:')
            for st in stats_list:
                if not self.method == 'syn':
                    outfile.write(', ' + str(st[0]) + ' ' + str(st[1]))
                else:
                    outfile.write(', ' + st)

            outfile.write('\n')
            for i_line, line in enumerate(results_all):
                if not self.method == 'syn':
                    outfile.write(util.time_to_str(
                                  self.sections[i_line].event.time))
                else:
                    outfile.write(self.events[i_line])
                for item in line:
                    outfile.write(', ' + str(item))
                outfile.write('\n')

        if plot:
            gp.plot_allgains(self, results_all, stats_list, directory, fn)
