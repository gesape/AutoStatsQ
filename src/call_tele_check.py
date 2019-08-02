import grond
import numpy as num
from pyrocko.guts import Object, List, Int, Float, String, dump_all

#--- !grond.StationCorrection                                                    
#codes: [CH, FUORN, '', Z]                                                       
#delay: 0.0                                                                      
#factor: 0.6542812498507474


class Correction_Statistics(Object):
    code = String.T()
    median_factor = Float.T()
    mean_factor = Float.T()
    stdev_factor = Float.T()
    n_ev = Int.T()   

def get_nslc(st):
    return '%s.%s' % (st.network, st.station)


def get_correction_statistcs(station_list, filename_list):
    
    mapping_stations_index = {}
    for i_st, st in enumerate(station_list):
        mapping_stations_index[get_nslc(st)] = i_st 

    corrections = num.empty((len(mapping_stations_index), len(filename_list)))
    corrections.fill(num.nan)

    for i_fn, fn in enumerate(filename_list):
        st_cor = grond.load_station_corrections(filename=fn)

        for sc in st_cor:
            #print(sc.codes[-1])
            if not sc.codes[-1] == 'Z':
                continue
            
            ii = mapping_stations_index['%s.%s' % (sc.codes[0], sc.codes[1])]
            corrections[ii, i_fn] = sc.factor


    ii_test = mapping_stations_index['CH.FUORN']
    print(corrections[ii_test, :])
    mean_correction_st = num.nanmean(corrections, axis=1)
    print(mean_correction_st[ii_test])
    stdev_correction_st = num.nanstd(corrections, axis=1)
    median_correction_st = num.nanmedian(corrections, axis=1)

    is_nan = num.isnan(corrections[:,:]).sum(axis=1)

    corrstats_list = []
    for m, i_m in mapping_stations_index.items():

        corrstats_list.append(Correction_Statistics(code=m,
                                                    median_factor=median_correction_st[i_m],
                                                    mean_factor=mean_correction_st[i_m],
                                                    stdev_factor=stdev_correction_st[i_m],
                                                    n_ev=len(filename_list)-is_nan[i_m]))

    #for c in corrstats_list:
    #    print(c)
    #print(corrstats_list)

    dump_all(corrstats_list, filename='test.yaml', stream=None)
    # dump...