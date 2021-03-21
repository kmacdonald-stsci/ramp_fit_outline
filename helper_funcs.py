#!/usr/bin/env python


from pprint import pprint
import sys

import numpy as np

from jwst.ramp_fitting.ramp_fit import ramp_fit
from jwst.ramp_fitting.ramp_fit import calc_num_seg
from jwst.datamodels import dqflags
from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]


DELIM = "-" * 70


def base_meta(model, ngroups, nrows, ncols, deltatime):
    model.meta.instrument.name = "MIRI"
    model.meta.instrument.detector = "MIRIMAGE"
    model.meta.instrument.filter = "F480M"

    model.meta.observation.date = "2015-10-13"

    model.meta.subarray.name = "FULL"
    model.meta.subarray.xstart = 1
    model.meta.subarray.ystart = 1
    model.meta.subarray.xsize = ncols
    model.meta.subarray.ysize = nrows

    model.meta.exposure.type = "MIR_IMAGE"
    model.meta.exposure.group_time = deltatime
    model.meta.exposure.frame_time = deltatime
    model.meta.exposure.ngroups = ngroups
    model.meta.exposure.group_time = deltatime
    model.meta.exposure.nframes = 1
    model.meta.exposure.groupgap = 0
    model.meta.exposure.drop_frames1 = 0


def zero_model(nints, ngroups, nrows, ncols, nframes, grouptime, deltatime):

    # Set up ramp model
    shape_4d = (nints, ngroups, nrows, ncols)
    shape_3d = (ngroups, nrows, ncols)
    shape_2d = (nrows, ncols)

    data = np.zeros(shape=shape_4d, dtype=np.float32)
    err = np.ones(shape=shape_4d, dtype=np.float32)
    pixdq = np.zeros(shape=shape_2d, dtype=np.uint32)
    gdq = np.zeros(shape=shape_4d, dtype=np.uint8)
    int_times = np.zeros((nints,))

    model = RampModel(
        data=data, err=err, pixeldq=pixdq, groupdq=gdq, int_times=int_times
    )

    # base_meta(model)  # Sets information in model.meta
    # Sets information in model.meta
    base_meta(model, ngroups, nrows, ncols, deltatime)

    return model


def gain_model(gain_level, nrows, ncols):
    # Set up gain model
    # gain_arr = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain_level
    gain_arr = np.full((nrows, ncols), gain_level, dtype=np.float64)
    gain = GainModel(data=gain_arr)

    gain.meta.instrument.name = "MIRI"
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = ncols
    gain.meta.subarray.ysize = nrows

    return gain


def rn_model(readnoise, nrows, ncols):
    # Set up read noise model
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float32)
    rnModel = ReadnoiseModel(data=read_noise)

    rnModel.meta.instrument.name = "MIRI"
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = ncols
    rnModel.meta.subarray.ysize = nrows

    return rnModel


# Need test for multi-ints near zero with positive and negative slopes
def setup_inputs(
        readnoise=10, gain_level=1,
        nints=1, ngroups=10, nrows=103, ncols=102, 
        nframes=1, grouptime=1.0, deltatime=1,):
    """
    Set up info.
    """
    model = zero_model(nints, ngroups, nrows, ncols, nframes, grouptime, deltatime)
    gain = gain_model(gain_level, nrows, ncols)
    rnModel = rn_model(readnoise, nrows, ncols)

    return model, rnModel, gain


if __name__=="__main__":
    rn, g, ni, ng, nr, nc, nf, gr, dt = 10, 1, 1, 3, 2, 2, 1, 1.0, 1
    model, rnModel, gain = setup_inputs(
        rn, g, ni, ng, nr, nc, nf, gr, dt)
