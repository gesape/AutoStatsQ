import os.path as op
import numpy as num
import hashlib

from io import BytesIO

from pyrocko.guts import Object, String
from pyrocko import cake, spit, util
from pyrocko.gf import Earthmodel1D


class Earthmodel(Object):
    id = String.T()


class CakeEarthmodel(Earthmodel):
    earthmodel_1d = Earthmodel1D.T()


def ttt_path(ct, ehash):
    return op.join(ct.dir_ttt, ehash + ".spit")


def ttt_hash(ct, earthmodel, phases, x_bounds, x_tolerance, t_tolerance):
    f = BytesIO()
    earthmodel.earthmodel_1d.profile("vp").dump(f)
    earthmodel.earthmodel_1d.profile("vs").dump(f)
    earthmodel.earthmodel_1d.profile("rho").dump(f)
    x_bounds.dump(f)
    x_tolerance.dump(f)
    num.array(t_tolerance).dump(f)
    s = f.getvalue()
    h = hashlib.sha1(s).hexdigest()
    p = ",".join(phase.definition() for phase in phases)
    f.close()
    return h + p


def get_ttt(ct, coords, depth_opt):
    x_bounds = num.array(
        [
            [ct.r_depth_min, ct.r_depth_max],
            # [ct.s_depth_min, ct.s_depth_max],
            [depth_opt[0], depth_opt[1]],
            [ct.dist_min, ct.dist_max],
        ],
        dtype=num.float,
    )

    # spatial tolerance for interpolating/ accuracy
    x_tolerance = num.array((ct.r_acc, ct.s_acc, ct.dist_acc))

    # temporal tolerance for interpolating/ accurancy
    t_tolerance = ct.t_acc  # s?

    earthmodel = CakeEarthmodel(
        id=ct.earthmodel_id,
        earthmodel_1d=cake.load_model(cake.builtin_model_filename(ct.earthmodel_id)),
    )

    interpolated_tts = {}

    for phase_def in ct.tabulated_phases:

        ttt_hash_ = ttt_hash(
            ct, earthmodel, phase_def.phases, x_bounds, x_tolerance, t_tolerance
        )

        fpath = ttt_path(ct, ttt_hash_)

        # make ttt for current phase, bounds, tolerance settings
        # or read if already existing!
        if not op.exists(fpath):

            def evaluate(args):
                receiver_depth, source_depth, x = args
                t = []
                rays = earthmodel.earthmodel_1d.arrivals(
                    phases=phase_def.phases,
                    distances=[x * cake.m2d],
                    zstart=source_depth,
                    zstop=receiver_depth,
                )

                for ray in rays:
                    t.append(ray.t)

                if t:
                    return min(t)
                else:
                    return None

            sptree = spit.SPTree(
                f=evaluate, ftol=t_tolerance, xbounds=x_bounds, xtols=x_tolerance
            )

            util.ensuredirs(fpath)
            sptree.dump(filename=fpath)
        else:
            sptree = spit.SPTree(filename=fpath)

        interpolated_tts["stored:" + str(phase_def.id)] = sptree
        intp_times = sptree.interpolate_many(coords)

        return intp_times
