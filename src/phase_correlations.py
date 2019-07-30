# read textfile with event_name, phase, tmin, tmax (t with repect to cake using ak135)
# compute cc between all stations for each event
# if Z comp low, check.

from pyrocko import cake, util, trace, orthodrome, pile
import numpy as num
import matplotlib.pyplot as plt
from math import ceil


def get_ev_phases(events, infile):

    use_events = []

    phase_dict = {}
    with open(infile, 'r') as f:
        for line in f.readlines():
            name, phase, tmin, tmax = line.split(', ')
            phase_dict[name] = (phase, tmin, tmax)

    use_events = [ev for ev in events if ev.name in phase_dict.keys()]

    return use_events, phase_dict


def prepare_trace(comp, pile, st, ev_time):

    tr = pile.all(tmin=ev_time, tmax=ev_time+3600,
                            trace_selector=lambda tr: tr.station == st.station and 
                                                      tr.channel == comp and
                                                      tr.location != '10')
    if not len(tr) == 1 or tr == []:
        print(tr)
        tr_ = []
    else:
        tr_ = tr[0].copy()
        tr_.bandpass(3, 0.01, 0.12)

    return tr_


def chop_and_correlate(tr1, tr2, tmin1, tmax1, tmin2, tmax2):
    try:
        tr1_ = tr1.chop(tmin=tmin1, tmax=tmax1, inplace=False)
        tr2_ = tr2.chop(tmin=tmin2, tmax=tmax2, inplace=False)
        cc = trace.correlate(tr1_, tr2_, mode='full', normalization='normal')
        #print(cc.absmax())
        #trace.snuffle([tr1_, tr2_])

    except:
        cc = None
    return cc


def rotate_and_cc(trs_rot, trs_fix, tmin_rot, tmax_rot, tmin_fix, tmax_fix):#, 
    #tr1_R, tr2_R,
    #tmin1, tmax1, tmin2, tmax2):

    tr_fix_R = [tr for tr in trs_fix if tr.channel == 'R'][0].chop(
                                                   tmin=tmin_fix, tmax=tmax_fix,
                                                   inplace=False)
    #tr_fix_T = [tr for tr in trs_fix if tr.channel == 'T'][0].chop(
    #                                               tmin=tmin_fix, tmax=tmax_fix)

    coefs = num.empty(len(range(0, 360)))
    coefs.fill(num.nan)


    #tr1_ = tr1_R.chop(tmin=tmin1, tmax=tmax1, inplace=False)
    #tr2_ = tr2_R.chop(tmin=tmin2, tmax=tmax2, inplace=False)

    for a in range(0, 360):#[0:1]:
        #print('a', a)

        trs = trace.rotate(trs_rot, a, in_channels=['R', 'T'], out_channels=['Rx', 'Tx'])
        tr_rot_Rx = [tr for tr in trs if tr.channel == 'Rx'][0].chop(
                                                   tmin=tmin_rot, tmax=tmax_rot,
                                                   inplace=False)
        #tr_rot_Tx = [tr for tr in trs if tr.channel == 'Tx'][0].chop(
        #                                           tmin=tmin_fix, tmax=tmax_fix)

        cc = trace.correlate(tr_rot_Rx, tr_fix_R, mode='full', normalization='normal')
        t, coef = cc.max()
        #if a == 0:
        #    print(coef)
        #    trace.snuffle([tr_rot_Rx, tr_fix_R, trs_rot[0], tr1_, tr2_])

        coefs[a] = coef
    #print(tr_rot_Rx.station, tr_fix_R.station)
    #print(num.nanmax(coefs), num.nanargmax(coefs))

    #plt.plot(coefs)
    #plt.show()
    return max(coefs), num.argmax(coefs)-1




def get_ccs(stations, events, phase_dict, work_dir, rotate=False):
    print('Starting correlation computations.')

    model = cake.load_model('ak135-f-continental.m')

    ccs = num.empty((3, len(events), len(stations), len(stations)))
    ccs.fill(num.nan)

    t_ex = 20

    if rotate:
        angles = num.empty(((len(events), len(stations), len(stations))))
        angles.fill(num.nan)
        ccs_rot = num.empty((len(events), len(stations), len(stations)))
        ccs_rot.fill(num.nan)
        #plt.imshow(angles[1,2,:,:])
        #plt.show()


    for i_ev, ev in enumerate(events):

        dists = num.asarray([float(orthodrome.distance_accurate50m_numpy(
                 ev.lat, ev.lon, st.lat, st.lon))
                 for st in stations])
        
        phases = [phase_dict[ev.name][0]]
        
        rays = model.arrivals(distances=list(dists*cake.m2d),
                                  phases=phases,
                                  zstart=ev.depth)
        arr_times = [ray.t for ray in rays]

        p = pile.make_pile('%s/rrd/%s' % (work_dir, 
                           util.time_to_str(ev.time).replace(' ', '_')))

        for i_st, st in enumerate(stations):
            tr1_Z = prepare_trace('Z', p, st, ev.time)
            tr1_R = prepare_trace('R', p, st, ev.time)
            tr1_T = prepare_trace('T', p, st, ev.time)

            tmin1 = ev.time+arr_times[i_st]+float(phase_dict[ev.name][1])-t_ex
            tmax1 = ev.time+arr_times[i_st]+float(phase_dict[ev.name][2])+t_ex

            for i_st2, st2 in enumerate(stations):
                tr2_Z = prepare_trace('Z', p, st2, ev.time)
                tr2_R = prepare_trace('R', p, st2, ev.time)
                tr2_T = prepare_trace('T', p, st2, ev.time)

                tmin2 = ev.time+arr_times[i_st2]+float(phase_dict[ev.name][1])-t_ex
                tmax2 = ev.time+arr_times[i_st2]+float(phase_dict[ev.name][2])+t_ex
                
                cc = chop_and_correlate(tr1_Z, tr2_Z, tmin1, tmax1, tmin2, tmax2)
                if cc is None:
                    continue
                t, coef = cc.absmax()
                ccs[0, i_ev, i_st, i_st2] = coef

                cc = chop_and_correlate(tr1_R, tr2_R, tmin1, tmax1, tmin2, tmax2)
                if cc is None:
                    continue                
                t, coef = cc.absmax()
                ccs[1, i_ev, i_st, i_st2] = coef

                cc = chop_and_correlate(tr1_T, tr2_T, tmin1, tmax1, tmin2, tmax2)
                if cc is None:
                    continue
                t, coef = cc.absmax()
                ccs[2, i_ev, i_st, i_st2] = coef

                if rotate:
                    #if '%s.%s' % (tr1_Z.network, tr1_Z.station) == rotate:
                    trs_rot = [tr2_R, tr2_T]
                    trs_fix = [tr1_R, tr1_T]
                    tmin_rot = tmin2
                    tmax_rot = tmax2
                    tmin_fix = tmin1
                    tmax_fix = tmax1
                    #elif '%s.%s' % (tr2_Z.network, tr2_Z.station) == rotate:
                    #    trs_rot = [tr2_R, tr2_T]
                    #    trs_fix = [tr1_R, tr1_T]
                    #    tmin_rot = tmin2
                    #    tmax_rot = tmax2
                    #    tmin_fix = tmin1
                    #    tmax_fix = tmax1
                    #else:

                    #else:
                    #    continue
                    coef, a_max = rotate_and_cc(trs_rot, trs_fix,
                                                tmin_rot, tmax_rot,
                                                tmin_fix, tmax_fix)#,
                                                #tr1_R, tr2_R, tmin1, tmax1, tmin2, tmax2)                      
                    angles[i_ev, i_st, i_st2] = a_max
                    ccs_rot[i_ev, i_st, i_st2] = coef


    if not rotate:
        return ccs
    else:
        return ccs, ccs_rot, angles


def plot_ccs(ccs, stations, events):
    # one plot for each component, one subplot for each event
    #print(ccs.shape)

    for i_comp, comp in enumerate(['Z', 'R', 'T']):
        fig, ax = plt.subplots(ncols=3, nrows=ceil(ccs.shape[1]/3), figsize=(10,10))
        for i_ev, ev in enumerate(events):
            if i_ev < ceil(ccs.shape[1]/3):
                i_c = 0
                i_r = i_ev
            elif i_ev < 2*ceil(ccs.shape[1]/3):
                i_c = 1
                i_r = i_ev - ceil(ccs.shape[1]/3)
            elif i_ev < 3*ceil(ccs.shape[1]/3):
                i_c = 2
                i_r = i_ev - 2*ceil(ccs.shape[1]/3)
            #print(i_ev, i_r, i_c)
            a = ax[i_r, i_c].imshow(ccs[i_comp, i_ev, :, :], interpolation='nearest',
                       vmin=0, vmax=+1, origin='lower')
            ax[i_r, i_c].set_title('%s' % util.time_to_str(ev.time)[0:10], fontsize=10)

            ax[i_r, i_c].set_yticks(num.arange(len(stations)))
            ax[i_r, i_c].set_yticklabels(['%s.%s' % (st.network, st.station) 
                                          for st in stations], fontsize=10)
            ax[i_r, i_c].set_xticks(num.arange(len(stations)))
            ax[i_r, i_c].set_xticklabels(['%s.%s' % (st.network, st.station) 
                                          for st in stations], fontsize=10,
                                          rotation=90)
        
        if i_ev == ccs.shape[1]-1:
            cbar = plt.colorbar(a, ax=ax[i_r, i_c])
            ticks = cbar.ax.get_yticklabels()
            cbar.ax.set_yticklabels(ticks, fontsize=10)

        plt.tight_layout()
        fig.savefig('cc_phases_%s.png' % comp)

        #plt.show()
        plt.close()


def plot_ccs_rot(ccs_rot, angles, stations, events):

    comp = 'R'

    fig, ax = plt.subplots(ncols=3, nrows=ceil(ccs_rot.shape[1]/3), figsize=(10,10))

    #x_start = 0
    #x_end = len(stations)-1
    #y_start = x_start
    #y_end = x_end
    #size = len(events)

    #jump_x = 1#(x_end - x_start) / (2.0 * size)
    #jump_y = 1#(y_end - y_start) / (2.0 * size)
    #x_positions = num.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
    #y_positions = num.linspace(start=y_start, stop=y_end, num=size, endpoint=False)
    #print(x_positions, y_positions)

    for i_ev, ev in enumerate(events):

        if i_ev < ceil(ccs_rot.shape[0]/3):
            i_c = 0
            i_r = i_ev
        elif i_ev < 2*ceil(ccs_rot.shape[0]/3):
            i_c = 1
            i_r = i_ev - ceil(ccs_rot.shape[0]/3)
        elif i_ev < 3*ceil(ccs_rot.shape[0]/3):
            i_c = 2
            i_r = i_ev - 2*ceil(ccs_rot.shape[0]/3)

        a = ax[i_r, i_c].imshow(ccs_rot[i_ev, :, :], interpolation='nearest',
                   vmin=0, vmax=+1, origin='lower')
        ax[i_r, i_c].set_title('%s' % util.time_to_str(ev.time)[0:10], fontsize=10)

        ax[i_r, i_c].set_yticks(num.arange(len(stations)))
        ax[i_r, i_c].set_yticklabels(['%s.%s' % (st.network, st.station) 
                                      for st in stations], fontsize=10)
        ax[i_r, i_c].set_xticks(num.arange(len(stations)))
        ax[i_r, i_c].set_xticklabels(['%s.%s' % (st.network, st.station) 
                                      for st in stations], fontsize=10,
                                      rotation=90)

        # add text (best rotation angles):
        #for y_index, y in enumerate(y_positions):
        #    for x_index, x in enumerate(x_positions):
        for i_st, st in enumerate(stations):
            for i_st2, st2 in enumerate(stations):
                label = angles[i_ev, i_st, i_st2]
                text_x = i_st #+ 0.5
                text_y = i_st #+ 0.5
                if not num.isnan(label):
                    ax[i_r, i_c].text(i_st, i_st2, label, color='black', ha='center', va='center')

        
        if i_ev == ccs_rot.shape[0]-1:
            cbar = plt.colorbar(a, ax=ax[i_r, i_c])
            ticks = cbar.ax.get_yticklabels()
            cbar.ax.set_yticklabels(ticks, fontsize=10)

    plt.tight_layout()
    fig.savefig('cc_phases_rot_%s.png' % comp)

    plt.show()