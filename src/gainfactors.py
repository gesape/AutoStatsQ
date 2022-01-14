#!/use/bin/env python

import os
import logging
import gc
import numpy as num
# from collections import defaultdict
from pyrocko import cake, util, trace
from pyrocko.gui import marker as pm
from pyrocko.orthodrome import distance_accurate50m
from pyrocko.guts import Object, Dict, String, Float, Int

from .gainplots import plot_allgains

# based on https://github.com/HerrMuellerluedenscheid/autogain

km = 1000.


class Gains(Object):
    trace_gains_mean = Dict.T(String.T(), Float.T())
    trace_gains_median = Dict.T(String.T(), Float.T())
    trace_gains_stdev = Dict.T(String.T(), Float.T())
    n_ev_used = Dict.T(String.T(), Int.T())
    ref_stats = String.T(optional=True)


def guess_nsl_template(code):
    if len(code) == 1 or isinstance(code, str):
        return '*.%s.*.*' % (code)
    elif len(code) == 2:
        return '*.%s.%s.*' % (code)
    elif len(code) == 3:
        return '%s.%s.%s.*' % (code)


class Section():
    """ Related to one event.
    All traces scale relative to average median abs
    max, to one reference station or to 1.0"""
    def __init__(self, event, stations):
        # Set Logger name
        self.logs = logging.getLogger('Section')

        self.stations = stations
        self.event = event
        self.reference_scale = None
        self.max_tr = {}
        self.relative_scalings = {}
        self.finished = False
        self.max_tr_syn = {}

    def finish(self, method, fband, taper, ev_counter):
        self.logs.debug('METHOD %s' % method)
        self.logs.debug(len(method))

        if len(method) == 2:
            reference_nsl = method[1][1]
            reference_nslc = list(filter(
                lambda x: util.match_nslc(guess_nsl_template(reference_nsl), x),
                                          self.max_tr.keys()))

            self.____reference_nslc = reference_nslc

            if not len(reference_nslc) == 1:
                self.logs.info(' No reference trace available. ' +
                               'remains unfinished: %s' % self.event.name)
                self.finished = False
            else:
                self.reference_scale = self.max_tr[reference_nslc[0]]
                self.set_relative_scalings()
                self.finished = True

        elif method == 'scale_one' or method == 'median_all_avail':
            self.reference_scale = 1.
            self.set_relative_scalings()
            self.finished = True

        elif method == 'syn_comp':
            #print(self.max_tr_syn)
            #print(self.max_tr)
            for nslc_id, maxA in self.max_tr.items():
                try:
                    #print(nslc_id[0:2])
                    self.relative_scalings[nslc_id] = maxA / self.max_tr_syn[nslc_id[0:2]]
                except:
                    self.logs.warning(' Data or synthetic data for comparison missing: %s' 
                                      % str(nslc_id[0:2]))
            #print(self.relative_scalings)
            self.finished = True


    def set_relative_scalings(self):
        for nslc_id, maxs in self.max_tr.items():
            self.relative_scalings[nslc_id] = maxs/self.reference_scale

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
    def __init__(self, data_pile, stations, events, arrT, snr_thresh,
                 component='Z', gain_rel_to='scale_one',
                 syn_data_pile=None):
        """
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

       :param data_pile: Pyrocko data pile 

        """
        # Set Logger name
        self.logs = logging.getLogger('AutoGain')

        self.method = gain_rel_to
        self.snr_thresh = snr_thresh
        self.component = component
        self.data_pile = data_pile
        self.stations = stations
        self.events = events
        self.all_nslc_ids = set()
        self.minmax = {}
        self.scaling_factors = {}
        self.sections = []
        self.results = None
        self._mean = None
        self._median = None  
        self._stdev = None  
        self._n_ev = None    
        self.arrT = arrT
        self.syn_data_pile = syn_data_pile

    def process(self, fband, taper, twd, debug):
        no_events = len(self.events)

        for i_ev, event in enumerate(self.events):
            tr_nslc_ids = []
            self.logs.info(' Processing event %s of %s' % (i_ev, no_events))
            section = Section(event, self.stations)
            skipped = 0
            unskipped = 0
            for i_s, s in enumerate(self.stations):
                dist = distance_accurate50m(event, s)
                arrival = self.arrT[i_ev, i_s]
                if num.isnan(arrival):
                    skipped += 1
                    self.logs.warning('skipped %s.%s %s' % (s.network, s.station, event.time))
                    continue
                else:
                    unskipped += 1
                selector = lambda tr: (s.network, s.station, self.component)\
                                   == (tr.network, tr.station, tr.channel) 

                tr_generator = self.data_pile.chopper(tmin=arrival -
                                                      twd[0],
                                                      tmax=arrival +
                                                      twd[1],
                                                      trace_selector=selector,
                                                      load_data=True)
                if self.method == 'syn_comp':
                    tr_syn_generator = self.syn_data_pile.chopper(tmin=arrival -
                                                      twd[0],
                                                      tmax=arrival +
                                                      twd[1],
                                                      trace_selector=selector,
                                                      load_data=True)

                for tr in tr_generator:
                    if not len(tr) > 1 and tr:
                        tr = tr[0]
                        if len(tr.ydata) > 0 and num.max(num.abs(tr.get_ydata())) != 0:
                            dtype = type(tr.ydata[0])
                            tr.ydata -= dtype(tr.get_ydata().mean())
                            # make SNR threshold here!
                            st_s = num.argmax(num.abs(tr.ydata))-10
                            snr = num.mean([y*y for y in tr.ydata[st_s:st_s+60]])/\
                                  num.mean([y*y for y in tr.ydata[0:60]])
                            if snr < self.snr_thresh:
                                continue
                            # mean(A*A_signal)/mean(A*A_noise)
                            tr.highpass(fband['order'], fband['corner_hp'])
                            tr.taper(taper, chop=False)
                            tr.lowpass(fband['order'], fband['corner_lp'])
                            
                            if debug is True:
                                self.logs.debug('SNR %s' % snr)
                                self.logs.debug('arrival time %s' % util.time_to_str(arrival))
                                trace.snuffle(tr, markers=[pm.Marker(nslc_ids=[tr.nslc_id], tmin=arrival, tmax=arrival+3)])

                            if num.max(num.abs(tr.get_ydata())) != 0:
                                section.max_tr[tr.nslc_id] = num.max(num.abs(tr.get_ydata()))
                                tr_nslc_ids.append(tr.nslc_id)                            


                    else:
                        for t in tr:
                            tt = t#[0]
                            if len(tt.ydata) > 0 and num.max(num.abs(tt.get_ydata())) != 0:
                                dtype = type(tt.ydata[0])
                                # print(tr.ydata, type(tr.ydata))
                                tt.ydata -= dtype(tt.get_ydata().mean())
                                st_s = num.argmax(num.abs(tt.ydata))-10
                                snr = num.mean([y*y for y in tt.ydata[st_s:st_s+60]])/\
                                      num.mean([y*y for y in tt.ydata[0:60]])
                                      
                                # print('SNR', snr) 
                                if snr < self.snr_thresh:
                                    continue                               
                                tt.highpass(fband['order'], fband['corner_hp'])
                                tt.taper(taper, chop=False)
                                tt.lowpass(fband['order'], fband['corner_lp'])

                                if debug is True:
                                    self.logs.debug('SNR %s' % snr)
                                    self.logs.debug('arrival time %s' % util.time_to_str(arrival))
                                    trace.snuffle(tt, markers=[pm.Marker(nslc_ids=[tt.nslc_id], tmin=arrival, tmax=arrival+3)])

                                if num.max(num.abs(tt.get_ydata())) != 0:
                                    section.max_tr[tt.nslc_id] = num.max(num.abs(tt.get_ydata()))
                                    tr_nslc_ids.append(tt.nslc_id)

                if self.method == 'syn_comp':
                    for tr in tr_syn_generator:
                        if not len(tr) > 1 and tr:
                            tr = tr[0]
                            if len(tr.ydata) > 0 and num.max(num.abs(tr.get_ydata())) != 0:
                                dtype = type(tr.ydata[0])
                                tr.ydata -= dtype(tr.get_ydata().mean())
                                st_s = num.argmax(num.abs(tr.ydata))-10
                                snr = num.mean([y*y for y in tr.ydata[st_s:st_s+60]])/\
                                      num.mean([y*y for y in tr.ydata[0:60]])
                                tr.highpass(fband['order'], fband['corner_hp'])
                                tr.taper(taper, chop=False)
                                tr.lowpass(fband['order'], fband['corner_lp'])
                                
                                if debug is True:
                                    self.logs.debug('SNR %s' % snr)
                                    self.logs.debug('arrival time %s' % util.time_to_str(arrival))
                                    trace.snuffle(tr, markers=[pm.Marker(nslc_ids=[tr.nslc_id], tmin=arrival, tmax=arrival+3)])

                                if num.max(num.abs(tr.get_ydata())) != 0:
                                    section.max_tr_syn[tr.nslc_id[0:2]] = num.max(num.abs(tr.get_ydata()))
                                    # tr_nslc_ids_syn.append(tr.nslc_id)


                        else:
                            for t in tr:
                                tt = t #[0]
                                if len(tt.ydata) > 0 and num.max(num.abs(tt.get_ydata())) != 0:
                                    dtype = type(tt.ydata[0])
                                    # print(tr.ydata, type(tr.ydata))
                                    tt.ydata -= dtype(tt.get_ydata().mean())
                                    st_s = num.argmax(num.abs(tt.ydata))-10
                                    snr = num.mean([y*y for y in tt.ydata[st_s:st_s+60]])/\
                                          num.mean([y*y for y in tt.ydata[0:60]])
                                          
                                    # print('SNR', snr) 
                                    if snr < self.snr_thresh:
                                        continue                               
                                    tt.highpass(fband['order'], fband['corner_hp'])
                                    tt.taper(taper, chop=False)
                                    tt.lowpass(fband['order'], fband['corner_lp'])

                                    if debug is True:
                                        self.logs.debug('SNR %s' % snr)
                                        self.logs.debug('arrival time %s' % util.time_to_str(arrival))
                                        trace.snuffle(tt, markers=[pm.Marker(nslc_ids=[tt.nslc_id], tmin=arrival, tmax=arrival+3)])

                                    if num.max(num.abs(tt.get_ydata())) != 0:
                                        section.max_tr_syn[tt.nslc_id[0:2]] = num.max(num.abs(tt.get_ydata()))
                                        #tr_nslc_ids_syn.append(tt.nslc_id)                    



                    #else:
                    #    print('no trace', s.network, s.station, tr, util.time_to_str(event.time))
                #break
                # print(i_s)
            self.logs.debug('skipped %s/%s' % (skipped, unskipped))

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

    @property
    def median(self):
        if self.results is None:
            self.congregate()
        if self._median is None and self.method != 'syn':
            self._median = dict(zip(map(lambda x: '.'.join(x),
                                      self.all_nslc_ids),
                                  num.nanmedian(self.results, axis=0)))
            # might rise a warning if for one station only nan are in results,
            # but works fine.
        if self._median is None and self.method == 'syn':
            nslcs = [st.split(' ')[1]+'.'+st.split(' ')[2]+'..'+'HHZ'
                     for st in self.stations]
            self._median = dict(zip(nslcs,
                                  num.nanmedian(self.results, axis=0)))
        return self._median

    @property
    def stdev(self):
        if self.results is None:
            self.congregate()
        if self._stdev is None and self.method != 'syn':
            self._stdev = dict(zip(map(lambda x: '.'.join(x),
                                      self.all_nslc_ids),
                                  num.nanstd(self.results, axis=0)))
            # might rise a warning if for one station only nan are in results,
            # but works fine.

        return self._stdev

    @property
    def n_ev(self):
        if self.results is None:
            self.congregate()
        if self._n_ev is None and self.method != 'syn':            
            nev_used = []
            for n_st in range(self.results.shape[1]):
                cnt = 0
                for n_ev in range(self.results.shape[0]):
                    if not num.isnan(self.results[n_ev, n_st]): 
                        cnt+=1
                nev_used.append(cnt)

            self._n_ev = dict(zip(map(lambda x: '.'.join(x),
                                      self.all_nslc_ids),
                                      nev_used))    
        return self._n_ev

    def save_mean(self, fn, directory):
        g = Gains()
        g.trace_gains_mean = self.mean
        if len(self.method) == 2:
            g.ref_stats = '%s %s' % (self.method[1][0], self.method[1][1])
        g.regularize()
        g.validate()
        g.dump(filename=os.path.join(directory, fn))

    def save_median(self, fn, directory):
        g = Gains()
        g.trace_gains_median = self.median
        if len(self.method) == 2:
            g.ref_stats = '%s %s' % (self.method[1][0], self.method[1][1])
        g.regularize()
        g.validate()
        g.dump(filename=os.path.join(directory, fn))

    def save_median_and_mean_and_stdev(self, fn, directory):
        g = Gains()
        g.trace_gains_median = self.median
        g.trace_gains_mean = self.mean
        g.trace_gains_stdev = self.stdev
        g.n_ev_used = self.n_ev 
        if len(self.method) == 2:
            g.ref_stats = '%s %s' % (self.method[1][0], self.method[1][1])
        g.regularize()
        g.validate()
        g.dump(filename=os.path.join(directory, fn))

    def save_single_events(self, fn, directory, plot=False):
        """ Save table with all gains (for all events + stations)
        """
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

        with open(os.path.join(directory, fn), 'w') as outfile:
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
            plot_allgains(self, results_all, stats_list, directory, fn)
