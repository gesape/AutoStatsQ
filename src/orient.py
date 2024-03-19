import numpy as num
import math, sys
import logging
import os
import datetime
from pyrocko import trace, pile
from pyrocko import util
from pyrocko.orthodrome import distance_accurate50m_numpy, azibazi
import matplotlib as mpl 
mpl.use('agg')
from matplotlib import pyplot as plt
from pyrocko.guts import Object, Dict, String, Float, List, Int, Tuple, load
from pyrocko.plot.automap import Map
import matplotlib.dates as mdates
from pyrocko.gui import marker as pm

import matplotlib as make_pile
from mpl_toolkits.axes_grid1 import make_axes_locatable


class dict_stats_rota(Object):
    CorrectAngl_perStat_median = Dict.T(String.T(), Float.T())
    CorrectAngl_perStat_mean = Dict.T(String.T(), Float.T())
    CorrectAngl_perStat_stdd = Dict.T(String.T(), Float.T())
    n_events = Dict.T(String.T(), Int.T())


class rota_ev_by_stat(Object):
    station = Tuple.T(3, String.T())
    ev_rota = Dict.T(String.T(), Int.T())


class dict_stats_all_rota(Object):
    dict_stats_all = List.T(rota_ev_by_stat.T())


class Event_sw(Object):
    station = Tuple.T(3, String.T())
    name = String.T()
    time = String.T()
    max_cc_angle = Int.T()
    max_cc_coef = Float.T()


class Polarity_switch(Object):
    switched_by_stat = List.T(Event_sw.T())


def plot_corr_time(nsl, filename, dir_ro):
    """
    Plot angle of max. cr-corr vs event time for each station.
    Results below cr-corr threshold are ignored.

    """
    angles_fromfile = load(filename=os.path.join(dir_ro, filename))
    # y_lim = (-180., 180.)

    for item in angles_fromfile.dict_stats_all:
        st = item.station        
        ev_list = []
        angle_list = []
        for i_ev, (ev, angle) in enumerate(item.ev_rota.items()): 
            ev_list.append(float(util.str_to_time(ev)))
            angle_list.append(float(angle))
        
        if ev_list and angle_list:
            ev_list_d = [datetime.datetime.strptime(util.time_to_str(ev), '%Y-%m-%d %H:%M:%S.%f') for ev in ev_list]
            absmax = max(angle_list, key=abs)
            if abs(absmax) < 60:
                absmax=60
            _median = num.median(angle_list)
            _mean = num.mean(angle_list)
            min_ev = num.min(ev_list)
            max_ev = num.max(ev_list)
            xvals = num.linspace(min_ev-1000, max_ev+1000, 100)
            xvals_d = [datetime.datetime.strptime(util.time_to_str(x), '%Y-%m-%d %H:%M:%S.%f') for x in xvals]

            fig, ax = plt.subplots(nrows=1, figsize=(10, 3))

            ax.plot(ev_list_d, angle_list, 'bo')

            ax.plot(xvals_d, [_median for i in range(len(xvals))], 'g--', label='median')
            ax.plot(xvals_d, [_mean for i in range(len(xvals))], 'r--', label='mean')
            
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.set_ylim((-abs(absmax)-1,abs(absmax)+1))
            ax.set_title('%s.%s' % (st[0],st[1]))
            ax.set_xlabel('Event date')
            ax.set_ylabel('Correction angle [°]')
            # xticks = ax.get_xticks()
            # print(xticks)
            # ax.set_xticklabels([util.time_to_str(x)[0:16] for x in xticks], rotation=70)
            plt.tick_params(labelsize=12)
            plt.tight_layout(rect=[0,0,0.8,1])
            # plt.show()
            
            fig.savefig(os.path.join(dir_ro, '%s_%s_overtime.png' % (st[0], st[1])))
            plt.close(fig)


def plot_corr_baz(nsl, filename_all, filename_stats, dir_ro, events, stations):
    """
    Plot angle of max. cr-corr vs event baz for each station.
    """

    stats_fromfile = load(filename=os.path.join(dir_ro, filename_stats))
    #print(stats_fromfile)
    angles_fromfile = load(filename=os.path.join(dir_ro, filename_all))

    times = [ev.time for ev in events]
    cmap = plt.get_cmap('viridis')
    ticks = [min(times), num.mean(times), max(times)]
    ticklabels = [util.tts(min(times))[:10], util.tts(num.mean(times))[:10], util.tts(max(times))[:10]]

    for item in angles_fromfile.dict_stats_all:
        st = item.station
        st_pyr = [s for s in stations if s.station == st[1] and s.network == st[0]][0]

        angle_list = []
        baz_list = []
        t_list = []

        for i_ev, (ev, angle) in enumerate(item.ev_rota.items()):
            #print(ev)
            angle_list.append(float(angle))
            try:
                ev_pyr = [e for e in events if abs((e.time - float(util.str_to_time(ev)))) < 100][0]
            except IndexError:
                continue
            bazi = azibazi( ev_pyr.lat, ev_pyr.lon, st_pyr.lat, st_pyr.lon)[1]
            t_list.append(ev_pyr.time)
            if bazi < 0:
                bazi = 360 + bazi
            baz_list.append(bazi)

        if angle_list and baz_list:
            absmax = max(angle_list, key=abs)
            if abs(absmax) < 60:
                absmax = 60
            #_median = num.median(angle_list)
            #_mean = num.mean(angle_list)
            #xvals = num.linspace(-180, 180, 100)
            try:
                _median = stats_fromfile.CorrectAngl_perStat_median['%s %s %s' % (st[0],st[1], st[2])]
                _mean =  stats_fromfile.CorrectAngl_perStat_mean['%s %s %s' % (st[0], st[1], st[2])]
                _stdev = stats_fromfile.CorrectAngl_perStat_stdd['%s %s %s' % (st[0], st[1], st[2])]
            
            except:
                continue
            
            xvals = num.linspace(-0, 360, 100)

            fig, ax = plt.subplots(figsize=(10, 3))
            ax.plot(xvals, [_mean for i in range(len(xvals))], c='black', ls='--', label='mean')
            ax.plot(xvals, [_median for i in range(len(xvals))], c='red', ls=':', label='median')
            ax.legend(loc='center right', bbox_to_anchor=(1, 0.5))
            ax.fill_between(xvals, [_mean-_stdev for i in range(len(xvals))], [_mean+_stdev for i in range(len(xvals))], facecolor='gray', alpha=0.2)
            im = ax.scatter(baz_list, angle_list, c=t_list, vmin=min(times), vmax=max(times), s=8)
            ax.set_ylim((-abs(absmax)-20,abs(absmax)+20))
            ax.set_xlim((0,+360))
            ax.set_title('%s.%s' % (st[0],st[1]))
            ax.set_xlabel('Station-Event BAZ [°]')
            ax.set_ylabel('Sensor correction angle [°]')
            plt.tick_params(labelsize=10)
            plt.yticks(num.arange(-60,70,20.))
            plt.tight_layout(rect=[0,0,0.8,1])
            #fig.colorbar(times, ax=ax)
            cbar = fig.colorbar(im, ax=ax, ticks=ticks)
            cbar.ax.set_yticklabels(ticklabels, fontsize=8)
            fig.savefig(os.path.join(dir_ro, '%s_%s_baz.png' % (st[0], st[1])))
            plt.close(fig)
        #sys.exit()         


def write_output(list_median_a, list_mean_a, list_stdd_a, list_switched,
                 n_ev, used_stats, dir_ro, ccmin):

    if list_switched:
        # write to yaml
        l_sw = [Event_sw(station=(s[0], s[1], s[2]), name=s[3],
                         time=s[4], max_cc_angle=s[5], max_cc_coef=s[6])
                for l in list_switched for s in l if s and s[6] > ccmin]

        sw = Polarity_switch(switched_by_stat=l_sw)
        sw.regularize()
        sw.validate()
        sw.dump(filename=os.path.join(dir_ro, 'Switched_Polarities.yaml'))

    if list_median_a and list_mean_a and list_stdd_a and n_ev:
        ns_list = ['%s %s %s' % (st[0], st[1], st[2]) for st in used_stats]
        perStat_median = dict(zip(map(lambda x: x, ns_list), list_median_a))
        perStat_mean = dict(zip(map(lambda x: x, ns_list), list_mean_a))
        perStat_std = dict(zip(map(lambda x: x, ns_list), list_stdd_a))
        perStat_nev = dict(zip(map(lambda x: x, ns_list), n_ev))

        dicts_rota = dict_stats_rota(CorrectAngl_perStat_median=perStat_median,
                                     CorrectAngl_perStat_mean=perStat_mean,
                                     CorrectAngl_perStat_stdd=perStat_std,
                                     n_events=perStat_nev)

        dicts_rota.regularize()
        dicts_rota.validate()
        dicts_rota.dump(filename=os.path.join(dir_ro, 'CorrectionAngles.yaml'))


def write_all_output_csv(list_all_angles, used_stats, dir_ro):
    list_rrr = []

    for st, ev_dict in zip(used_stats, list_all_angles):
        rrr = rota_ev_by_stat(station=(st[0], st[1], st[2]),
                              ev_rota=ev_dict)
        list_rrr.append(rrr)

    dict_save_rot_st_ev = dict_stats_all_rota(dict_stats_all=list_rrr)

    dict_save_rot_st_ev.dump(filename=os.path.join(dir_ro, 'AllCorrectionAngles.yaml'))


def get_m_angle_switched(cc_i_ev_vs_rota, catalog, st, ccmin):
    """
    1) did polarity swith occur? for single events, say if
    max c-c is closer to 180 deg than to 0 deg! (only trusted if
    coef > 0.80)
    2) provide median of correction angle associated to max-cc,
    only use those with max cc-coef > 0.85 or sth like that
    """
    angles = []
    values = []
    switched = []

    for i_ev, ev in enumerate(catalog):
        maxcc_value = num.max(cc_i_ev_vs_rota[i_ev, :])

        if not num.isnan(maxcc_value):
            maxcc_angle = -180 + num.argmax(cc_i_ev_vs_rota[i_ev, :])
            angles.append(maxcc_angle)
            values.append(maxcc_value)

            if abs(0. - maxcc_angle) > 90:
                switched.append((st.network, st.station, st.location,
                                 ev.name,
                                 util.time_to_str(ev.time),
                                 maxcc_angle, maxcc_value))

    list_v_above_ccmin = [a for (a, v) in zip(angles, values) if v > ccmin]
    median_a = num.median(list_v_above_ccmin)
    # use a vector mean instead!
    # mean_a = num.mean(list_v_above_ccmin)
    if len(list_v_above_ccmin) > 0:

        sum_v = num.asarray((0, 0))

        for a in list_v_above_ccmin:
            x = num.cos(num.deg2rad(a))
            y = num.sin(num.deg2rad(a))

            sum_v = sum_v + num.asarray((x,y)) #/ num.linalg.norm((x,y))) 
        
        # normalize
        sum_v = sum_v / len(list_v_above_ccmin)

        if sum_v[0] > 0:
            mean_a = num.rad2deg(num.arctan(sum_v[1]/sum_v[0]))
        elif sum_v[0] < 0 and sum_v[1] >= 0:
            mean_a = num.rad2deg(num.arctan(sum_v[1]/sum_v[0])) + 180.
        elif sum_v[0] < 0 and sum_v[1] < 0:
            mean_a = num.rad2deg(num.arctan(sum_v[1]/sum_v[0])) - 180.
        elif sum_v[0] == 0 and sum_v[1] > 0:
            mean_a = 90.
        elif sum_v[0] == 0 and sum_v[1] < 0:
            mean_a = 270.

        # print('mean_a', mean_a)       

        # vector standard deviation
        # std_a = num.std([a for (a, v) in zip(angles, values) if v > ccmin])
        n_ev = len(list_v_above_ccmin)
        # std_a = num.sqrt((sum([(x_i - mean_a) * (x_i - mean_a) for x_i in list_v_above_ccmin])) / n_ev)
        # print(list_v_above_ccmin)
        if len(list_v_above_ccmin) > 1:
            sum_d_xi_xm = 0
            for x_i in list_v_above_ccmin:
                # print(x_i)

                x_i_x = num.cos(num.deg2rad(x_i))
                x_i_y = num.sin(num.deg2rad(x_i))
                mean_x = num.cos(num.deg2rad(mean_a))
                mean_y = num.sin(num.deg2rad(mean_a))
                len_x_i = num.sqrt(x_i_x*x_i_x + x_i_y*x_i_y)
                len_mean = num.sqrt(mean_x*mean_x + mean_y*mean_y)

                # print(num.asarray((num.cos(num.deg2rad(x_i)),
                #                       num.sin(num.deg2rad(x_i)))))
                # print(num.asarray((num.cos(num.deg2rad(mean_a)),
                #                       num.sin(num.deg2rad(mean_a)))))
                # d_xi_m = num.abs(num.asarray((num.cos(num.deg2rad(x_i)),
                #                       num.sin(num.deg2rad(x_i)))) -\
                #          num.asarray((num.cos(num.deg2rad(mean_a)),
                #                       num.sin(num.deg2rad(mean_a)))))
                # print(d_xi_m)
                # if d_xi_m[0] > 0:
                #     phi_d = num.rad2deg(num.arctan(d_xi_m[1]/d_xi_m[0]))
                # elif d_xi_m[0] < 0 and d_xi_m[1] >= 0:
                #     phi_d = num.rad2deg(num.arctan(d_xi_m[1]/d_xi_m[0])) + 180.
                # elif d_xi_m[0] < 0 and d_xi_m[1] < 0:
                #     phi_d = num.rad2deg(num.arctan(d_xi_m[1]/d_xi_m[0])) - 180.
                # elif d_xi_m[0] == 0 and d_xi_m[1] > 0:
                #     phi_d = 90.
                # elif d_xi_m[0] == 0 and d_xi_m[1] < 0:
                #     phi_d = 270.
                # elif d_xi_m[0] == 0 and d_xi_m[1] == 0:
                #     phi_d = 0.
                
                phi_d = num.rad2deg(num.arccos(x_i_x/len_x_i * mean_x/len_mean
                                               + x_i_y/len_x_i * mean_y/len_mean))
                # Fallunterscheidungen?

                # print(phi_d)
                sum_d_xi_xm += phi_d*phi_d
                # print(sum_d_xi_xm)
            std_a = num.sqrt(sum_d_xi_xm / n_ev)
            # print('new std', std_a)
            # print('old std',num.std([a for (a, v) in zip(angles, values) if v > ccmin]))
        else:
            std_a = num.nan    

    else:
        median_a = num.nan
        mean_a = num.nan
        std_a = num.nan
        n_ev = 0

    return median_a, mean_a, std_a, switched, n_ev


def get_m_angle_all(cc_i_ev_vs_rota, catalog, st, ccmin):

    dict_ev_angle = {}

    for i_ev, ev in enumerate(catalog):
        maxcc_value = num.max(cc_i_ev_vs_rota[i_ev, :])

        if not num.isnan(maxcc_value):
            maxcc_angle = -180 + num.argmax(cc_i_ev_vs_rota[i_ev, :])

            if maxcc_value > ccmin:
                dict_ev_angle[util.time_to_str(ev.time)] = int(maxcc_angle)

    return dict_ev_angle


def get_tr_by_cha(pile, start_twd, end_twd, loc, cha):
    tr = pile.all(
        tmin=start_twd,
        tmax=end_twd,
        trace_selector=lambda tr: tr.nslc_id[3] == cha and tr.nslc_id[2] == loc,
        want_incomplete=True)

    return tr


def max_or_min(c):
    """
    for testing... not used because valid cc is used which returns only
    one value (no time shifts!)
    Get time and value of maximum of the absolute of data.
    """
    tmi, mi = c.min()
    tma, ma = c.max()
    if abs(mi) > abs(ma):
        return tmi, mi
    else:
        return tma, ma


def plot_ccdistr_each_event(cc_i_ev_vs_rota, catalog, rot_angles, st, loc, dir_ro):
    """
    Plots max. cc-coef vs. rotation angle for each event in subplots.
    rather for debugging than for including into normal testing workflow!

    :param cc_i_ev_vs_rota: num.array [i_ev, rot_a] containing max. cc-coef for
                            each rotation angle (rot_a) and event (i_ev)
    :param catalog: list of pyrocko events used for analysis
    :param rot_angles: range of rotation angles
    :param st: current station (pyrocko Station)
    :param dir_ro: output directory
    """
    n_ev = len(catalog)
    nrows = math.ceil(n_ev / 5)
    if not n_ev < 5:
        ncols = 5
    else:
        ncols = n_ev

    y_lim = (-1., 1.)
    x_lim = (-180, 180)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*2, nrows*2))

    for (i_row, row), ev in zip(enumerate(cc_i_ev_vs_rota), catalog):
        i_x = int(i_row/ncols)
        i_y = int(i_row % ncols)
        ev_time_str = util.time_to_str(ev.time)[0:10]

        if nrows != 1:
            ax[i_x, i_y].set_title(ev_time_str, fontsize=10)
            if i_x == nrows-1:
                ax[i_x, i_y].set_xlabel('Correction angle [deg]', fontsize=8)
            if i_y == 0:
                ax[i_x, i_y].set_ylabel('C-c coef.', fontsize=8)
            ax[i_x, i_y].plot(rot_angles, row, 'k')
            ax[i_x, i_y].set_xlim(x_lim)
            ax[i_x, i_y].set_ylim(y_lim)
            ax[i_x, i_y].set_xticks([-180, 0, 180])
            ax[i_x, i_y].tick_params(labelsize=6)

        elif nrows == 1 and ncols != 1:
            ax[i_y].set_title(ev_time_str, fontsize=10)
            if i_x == nrows-1:
                ax[i_y].set_xlabel('Correction angle [deg]', fontsize=8)
            if i_y == 0:
                ax[i_y].set_ylabel('C-c coef.', fontsize=8)
            ax[i_y].plot(rot_angles, row, 'k')
            ax[i_y].set_xlim(x_lim)
            ax[i_y].set_ylim(y_lim)
            ax[i_y].set_xticks([-180, 0, 180])
            ax[i_y].tick_params(labelsize=6)

        elif nrows == 1 and ncols == 1:
            ax.set_title(ev_time_str, fontsize=10)
            ax.set_xlabel('Correction angle [deg]', fontsize=8)
            ax.set_ylabel('C-c coef.', fontsize=8)
            ax.plot(rot_angles, row, 'k')
            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)
            ax.set_xticks([-180, 0, 180])
            ax.tick_params(labelsize=6)


    if nrows*ncols > n_ev:
        dif = nrows*ncols - n_ev
        for i in range(dif):
            if nrows != 1:
                ax[i_x, i_y+i+1].set_xlabel('Correction angle [deg]', fontsize=8)
                ax[i_x, i_y+i+1].set_xticks([])
                ax[i_x, i_y+i+1].set_yticks([])
                ax[i_x, i_y+i+1].axis('off')
            else:
                ax[i_x].set_xlabel('Correction angle [deg]', fontsize=8)
                ax[i_x].set_xticks([])
                ax[i_x].set_yticks([])
                ax[i_x].axis('off')
    plt.tight_layout()
    # plt.show()
    fig.savefig(os.path.join(dir_ro, '%s_%s_%s_distr.png' % (st.network, st.station, loc)))
    plt.close(fig)


def prep_orient(datapath, st, i_st, nst, loc, catalog, dir_ro, v_rayleigh,
                bp, dt_start, dt_stop, ccmin=0.80,
                plot_heatmap=False,  plot_distr=False,
                debug=False):
    """
    Perform orientation analysis using Rayleigh waves, main function.

    time wdw: 20s before 4.0km/s arrival and 600 s afterwards
    (Stachnik et al. 2012)
    - compute radial component for back values of 0 to 360 deg
    - for each c-c of hilbert(R) with Z comp.
    - call plotting functions and/or write results to file

    :param datapath: path to rrd data
    :param st: current station (pyrocko station object)
    :param catalog: list of pyrocko events used for analysis
    :param dir_ro: output directory
    :param plot_heatmap: bool, optional
    :param plot_distr: bool, optional
    """
    logs = logging.getLogger('prep_orient')
    #logs.setLevel('DEBUG')
    st_data_pile = pile.make_pile(datapath, regex='%s_%s_' % (st.network, st.station),
                                  show_progress=False)
    n_ev = len(catalog)

    if st_data_pile.tmin is not None and st_data_pile.tmax is not None:

        # calculate dist between all events and current station
        r_arr_by_ev = num.empty(n_ev)
        ev_lats = num.asarray([ev.lat for ev in catalog])
        ev_lons = num.asarray([ev.lon for ev in catalog])
        dists = distance_accurate50m_numpy(a_lats=ev_lats, a_lons=ev_lons,
                                           b_lats=st.lat, b_lons=st.lon,
                                           implementation='c')
        r_arr_by_ev = (dists/1000.) / v_rayleigh
        cc_i_ev_vs_rota = num.empty((n_ev, 360))
        rot_angles = range(-180, 180, 1)
        for i_ev, ev in enumerate(catalog):
            arrT = ev.time + r_arr_by_ev[i_ev]

            start_twd1 = ev.time
            end_twd1 = arrT + 1800

            trZ = get_tr_by_cha(st_data_pile, start_twd1, end_twd1, loc, 'Z')
            trR = get_tr_by_cha(st_data_pile, start_twd1, end_twd1, loc, 'R')
            trT = get_tr_by_cha(st_data_pile, start_twd1, end_twd1, loc, 'T')

            start_twd2 = ev.time + r_arr_by_ev[i_ev] - dt_start
            end_twd2 = arrT + dt_stop

            if len(trZ) == 1 and len(trR) == 1 and len(trT) == 1:
                trZ = trZ[0]
                trR = trR[0]
                trT = trT[0]
                # debugging - window selection:
                if debug is True:
                    trace.snuffle([trZ,trR,trT], markers=
                        [pm.Marker(nslc_ids=[trZ.nslc_id, trR.nslc_id, trT.nslc_id],
                                   tmin=start_twd2, tmax=end_twd2),
                         pm.Marker(nslc_ids=[trZ.nslc_id, trR.nslc_id, trT.nslc_id],
                                   tmin=arrT, tmax=arrT+3)])

            else:
                cc_i_ev_vs_rota[i_ev, :] = num.nan
                continue

            try:
                trZ.bandpass(bp[0], bp[1], bp[2])
                trZ.chop(tmin=start_twd2, tmax=end_twd2)
            except trace.NoData:
                logs.warning('no data %s %s %s' % (trZ, trR, trT))
                continue

            for i_r, r in enumerate(rot_angles):
                print('Station: %5d/%s, Event: %5d/%s, rotation angle [deg]: %5d' 
                       % (i_st, nst, i_ev, n_ev, r), end='\r')
                try:
                    rot_2, rot_3 = trace.rotate(traces=[trR, trT], azimuth=r,
                                            in_channels=['R', 'T'],
                                            out_channels=['2', '3'])
                except ValueError:
                    logs.warning('Rotation failed, %s.%s, %s' % (st.network, st.station, r))
                    cc_i_ev_vs_rota[i_ev, i_r] = num.nan
                    continue

                rot_2_y = rot_2.ydata
                rot_2_hilb = num.imag(trace.hilbert(rot_2_y, len(rot_2_y)))
                rot_2_hilb_tr = trace.Trace(deltat=rot_2.deltat,
                                            ydata=rot_2_hilb,
                                            tmin=rot_2.tmin)
                # problem: rot_2 and rot_2_hilb look exactly the same!
                # --> no phase shift. why? should be num.imag!!!
                # trace.snuffle([rot_2, rot_2_hilb_tr])
                rot_2_hilb_tr.bandpass(bp[0], bp[1], bp[2])
                rot_2_hilb_tr.chop(tmin=start_twd2, tmax=end_twd2)

                # if st.station == 'RORO' and r == 0:
                #     trace.snuffle([rot_2_hilb_tr, trZ])
                # normalize traces
                trZ.ydata /= abs(max(trZ.ydata))
                rot_2_hilb_tr.ydata /= abs(max(rot_2_hilb_tr.ydata))

                c = trace.correlate(trZ, rot_2_hilb_tr,
                                    mode='valid',
                                    normalization='normal')

                t, coef = c.max()
                t2, coef2 = max_or_min(c)
                '''
                if st.station == 'MATE' and r == 0:
                    print(i_ev, ev.name, ev.depth)
                    print(r, t, coef, t2, coef2)
                    trace.snuffle([trZ, trR, rot_2_hilb_tr])
                '''
                cc_i_ev_vs_rota[i_ev, i_r] = coef
        '''
        if st.station == 'MATE':
            for i_ev in range(n_ev):
                print(num.argmax(cc_i_ev_vs_rota[i_ev,:]),
                      num.max(cc_i_ev_vs_rota[i_ev,:]))
        '''

        if plot_heatmap is True:
            logs.debug('Plotting heatmap for station %s.%s' % (st.network, st.station))
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 2))

            cax = ax.imshow(cc_i_ev_vs_rota, interpolation='nearest',
                            vmin=-1.0, vmax=1.0,
                            aspect='auto', extent=[-180, 180, n_ev, 0],
                            cmap='binary')
            ax.set_ylabel('i_ev')
            ax.set_xlabel('Correction angle (deg)')
            ax.set_title('%s %s' % (st.network, st.station))
            cbar = fig.colorbar(cax, ticks=[0, 0.5, 1.0],
                                orientation='horizontal',
                                fraction=0.05, pad=0.5)
            cbar.ax.set_xticklabels(['0', '0.5', '1.0'])
            plt.tight_layout()
            # plt.show(fig)
            fig.savefig(os.path.join(dir_ro, '%s_%s_%s_rot_cc_heatmap.png' % (st.network, st.station, loc)))
            plt.close()

        if plot_distr is True:
            logs.debug('Plotting distributions for station %s.%s' % (st.network, st.station))
            plot_ccdistr_each_event(cc_i_ev_vs_rota, catalog,
                                    rot_angles, st, loc, dir_ro)

        median_a, mean_a, std_a, switched, n_ev =\
            get_m_angle_switched(cc_i_ev_vs_rota, catalog, st, ccmin)

        dict_ev_angle = get_m_angle_all(cc_i_ev_vs_rota, catalog, st, ccmin)

        return median_a, mean_a, std_a, switched, dict_ev_angle, n_ev


def plot_corr_angles(ns, st_lats, st_lons, orientfile, dir_orient,
                     pl_options, pl_topo, mapsize, outformat, ls=False):
    '''
    Plot correction angles of all stations on a map. nans are igored.
    Values for plotting are read from file which was automatically prepared in
    the orient section.

    :param ls: label position [lon (symbol), lat (symbol + label), lon (label)]
    '''

    logs = logging.getLogger('plot_corr_angles')

    gmtconf = dict(
                   MAP_TICK_PEN_PRIMARY='1.25p',
                   MAP_TICK_PEN_SECONDARY='1.25p',
                   MAP_TICK_LENGTH_PRIMARY='0.2c',
                   MAP_TICK_LENGTH_SECONDARY='0.6c',
                   FONT_ANNOT_PRIMARY='20p,1,black',
                   FONT_LABEL='20p,1,black',
                   PS_CHAR_ENCODING='ISOLatin1+',
                   MAP_FRAME_TYPE='fancy',
                   FORMAT_GEO_MAP='D',
                   PS_PAGE_ORIENTATION='portrait',
                   MAP_GRID_PEN_PRIMARY='thinnest,0/50/0',
                   MAP_ANNOT_OBLIQUE='6')

    angles_fromfile = load(filename=os.path.join(dir_orient, orientfile))

    angle_no_nan = []
    lat_no_nan = []
    lon_no_nan = []
    angle_no_nan_u = []
    lat_no_nan_u = []
    lon_no_nan_u = []
    stats_no_nan = []
    stats_no_nan_u = []


    for i_ns, ns_now in enumerate(ns):
        for l in ['00', '', '01', '10', '60']:
            try:
                ns_now_ = '%s %s %s' % (ns_now[0], ns_now[1], l)
                a = angles_fromfile.CorrectAngl_perStat_median[ns_now_]
                nev = angles_fromfile.n_events[ns_now_]
                if a > -181.0 and a < 180.0:  # not nan
                    if nev >= 5:
                        angle_no_nan.append(0.0-a)
                        stats_no_nan.append(ns_now[1])
                        lat_no_nan.append(st_lats[i_ns])
                        lon_no_nan.append(st_lons[i_ns])
                    else:
                        angle_no_nan_u.append(0.0-a)
                        stats_no_nan_u.append(ns_now[1])                        
                        lat_no_nan_u.append(st_lats[i_ns])
                        lon_no_nan_u.append(st_lons[i_ns])
            except KeyError:
                continue

    # Generate the basic map
    m = Map(
        lat=pl_options[0],
        lon=pl_options[1],
        radius=pl_options[2],
        width=mapsize[0],
        height=mapsize[1],
        show_grid=False,
        show_topo=pl_topo,
        #color_dry=(143, 188, 143), #(238, 236, 230),
        topo_cpt_wet='white_sea_land',#'light_sea_uniform',
        topo_cpt_dry='light_land_uniform',
        illuminate=True,
        illuminate_factor_ocean=0.15,
        show_rivers=True,
        show_plates=False,
        gmt_config=gmtconf)

    # Draw some larger cities covered by the map area
    # m.draw_cities()

    # Draw max. amplitudes at station locations as colored circles
    cptfile = 'tempfile2.cpt'
    abs_angs = list(num.abs(angle_no_nan))
    try:
        m.gmt.makecpt(
                    C=pl_options[3],
                    T='%f/%f' % (0.1, 180.),
                    out_filename=cptfile)#, suppress_defaults=True)
    except:
        try:
            m.gmt.makecpt(
                    C='split',
                    T='%f/%f' % (0.1, 180.),
                    out_filename=cptfile)#, suppress_defaults=True)
            logs.warning('Could not find gmt cptfile, using split instead.')
        except:
            logs.error('Could not make gmt cpt file for map.')

    # m.gmt.makecpt(
    #             C='/home/gesap/Documents/CETperceptual_GMT/CET-D8.cpt',
    #             T='%f/%f' % (0.1, 180.),
    #             out_filename=cptfile)#, suppress_defaults=True)


    # same length for every vector:
    length = [1.5 for a in range(len(lat_no_nan))]
    length_u = [0.7 for a in range(len(lat_no_nan_u))]

    # angle_zero = [1.5 for a in range(len(lat_no_nan))]
    # angle_zero_u = [1.5 for a in range(len(lat_no_nan_u))]

    # plot obtained rotation vectors:
    #m.gmt.psxy(in_columns=(lon_no_nan, lat_no_nan, angle_zero, length),
    #           S='V0.5c+jc', W='0.07c,black',
    #           *m.jxyr)
    #m.gmt.psxy(in_columns=(lon_no_nan_u, lat_no_nan_u, angle_zero_u, length_u),
    #       S='V0.5c+jc', W='0.07c,black',
    #       *m.jxyr)

    m.gmt.psxy(in_columns=(lon_no_nan_u, lat_no_nan_u, angle_no_nan_u, length_u),
           S='V0.4c+jc+eA', W='0.05c,black',
           *m.jxyr)
    
    m.gmt.psxy(in_columns=(lon_no_nan, lat_no_nan, abs_angs, angle_no_nan, length),
               C=cptfile, S='V0.7c+jc+eA', W='0.1c+cl',
               *m.jxyr)


    # add handmade label
    if ls:
        m.gmt.psxy(in_columns=([ls[0]], [ls[1]], [ls[3]], [0.9]),
                   S='V0.5c+eA', W='0.07c,red', *m.jxyr)
        labels = ['Sensor orientation']
        lat_lab = [ls[1]]
        lon_lab = [ls[2]]
        for i in range(len(labels)):
            m.add_label(lat_lab[i], lon_lab[i], labels[i])
    
    # add a colorbar
    B_opt_psscale = 'xaf'
    m.gmt.psscale(
                 B=B_opt_psscale+'+l abs. misorientation [deg]',
                 D='x9c/6c+w12c/0.5c+jTC+h', # 'x9c/17c+w12c/0.5c+jTC+h'
                 C=cptfile)
    
    # add station labels

    has_label = []

    for i in range(len(stats_no_nan)):
        if stats_no_nan[i] not in has_label:
            m.add_label(lat_no_nan[i], lon_no_nan[i], stats_no_nan[i])    
            has_label.append(stats_no_nan[i])

    has_label = []

    for i in range(len(stats_no_nan_u)):
        if stats_no_nan_u[i] not in has_label:
            m.add_label(lat_no_nan_u[i], lon_no_nan_u[i], stats_no_nan_u[i])    
            has_label.append(stats_no_nan_u[i])

    m.save(os.path.join(dir_orient, 'map_orient.%s' % outformat))
