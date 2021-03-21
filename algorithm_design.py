#!/usr/bin/env python

import multiprocessing
import sys

import numpy as np

from helper_funcs import setup_inputs

from jwst.datamodels import dqflags
from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel

DELIM = "-" * 70
SDELIM = "*" * 70

"""
Receive 4-D RampModel with processing information:
"""


class IntermediateComputations:
    def __init__(self, data):
        self.orig_dim = data.shape
        # Add variables here

    def set_stuff(self, stuff):
        pass


class ModelArray:
    def __init__(self, orig_model, data, err, groupdq, pixeldq):
        self.orig_model = orig_model
        self.data = data
        self.err = err
        self.groupdq = groupdq
        self.pixeldq = pixeldq


def ramp_fit(model, proc_options):
    """
    Top level ramp fitting designed for 4-D processing  with dimensions:
        (integrations, groups, rows, columns)
    """
    icomp = IntermediateComputations(model.data)
    miri_answer = miri_correction(model)
    model_4d = ModelArray(model, model.data, model.err, model.groupdq, model.pixeldq)

    return ramp_fit_4d(model_4d, proc_options, icomp)


def miri_correction(model):
    return None


def ramp_fit_4d(model, proc_options, icomp):
    """
    Rampfitting a 4-D cube.
    """
    print("    *** ramp_fit_4d")
    num_slices = compute_multiprocessing_slices(proc_options)
    if num_slices == 1:
        return across_integration_ramp_fitting(model, proc_options, icomp)
    else:
        return ramp_fit_4d_multi(model, proc_options, icomp, num_slices)


def across_integration_ramp_fitting(model, proc_options, icomp):
    print("    *** across_integration_ramp_fitting")
    primary_int = []
    rateint_int = []
    optres_int = []
    print(f"The number of integrations: {model.data.shape[0]}")
    for integration in range(model.data.shape[0]):
        print(f"--> Integration = {integration}")
        model_3d, icomp_3d = get_integration(model, icomp, integration)
        primary, rateint, optres = ramp_fit_3d(model_3d, proc_options, icomp_3d)
        primary_int.append(primary)
        rateint_int.append(rateint)
        optres_int.append(optres)

    return combine_integrations(primary_int, rateint_int, optres_int)


def combine_integrations(primary_int, rateint_int, optres_int):
    return None, None, None
# *******************************************************************
def ramp_fit_3d(model, proc_options, icomp):
    """
    Rampfitting a 3-D cube.
    """
    print("    *** ramp_fit_3d")
    print(DELIM)
    print(model.data.shape)
    print(DELIM)
    print(model.data)
    print(DELIM)

    return None, None, None
    # return primary, rateint, opt_res


# *******************************************************************


def get_integration(model, icomp, n_int):
    """
    model: ModelArray
        Model information and image arrays.

    icomp: IntermediateComputations
        Intermediate computations used to fit ramp.

    n_int: int
        The current integration.
    """
    # Reduce by one dimension
    data = model.data[n_int, :, :, :]
    err = model.err[n_int, :, :, :]
    groupdq = model.groupdq[n_int, :, :, :]
    pixeldq = model.pixeldq

    model_3d = ModelArray(model.orig_model, data, err, groupdq, pixeldq)

    icomp_3d = IntermediateComputations(model_3d.data)

    return model_3d, icomp_3d


def compute_multiprocessing_slices(proc_options):
    """
    Compute multiprocessing slices per processor.
    """
    max_cores = proc_options["max_cores"]
    if max_cores == "none":
        number_slices = 1
    else:
        num_cores = multiprocessing.cpu_count()
        # log.debug(f'Found {num_cores} possible cores to use for ramp fitting')
        if max_cores == "quarter":
            number_slices = num_cores // 4 or 1
        elif max_cores == "half":
            number_slices = num_cores // 2 or 1
        elif max_cores == "all":
            number_slices = num_cores
        else:
            number_slices = 1

    return number_slices


def ramp_fit_4d_multi(model, proc_options, icomp, num_slices):
    print("    *** ramp_fit_4d_multi")


def get_row_slices(nrows, num_slices):
    """"""
    rows_per_slice = nrows // num_slices
    remainder = nrows % num_slices

    rslices = [rows_per_slice] * num_slices
    row_slices = np.array(rslices, dtype=np.int32)

    if remainder > 0:
        row_slices[:remainder] += 1

    return row_slices


# --------------------------------------------------------------------
def test_get_row_slices():
    nrows = 105
    nslices = 10

    rslices = get_row_slices(nrows, nslices)
    trows = sum(rslices)

    assert trows == nrows
    assert rslices[0] - rslices[-1] == 1


def test_compute_multiprocessing_slices():
    ncores = multiprocessing.cpu_count()

    proc_options = {"max_cores": "none"}
    num_slices = compute_multiprocessing_slices(proc_options)
    assert num_slices == 1

    proc_options = {"max_cores": "all"}
    num_slices = compute_multiprocessing_slices(proc_options)
    assert num_slices == ncores

    proc_options = {"max_cores": "half"}
    num_slices = compute_multiprocessing_slices(proc_options)
    assert num_slices == ncores // 2

    proc_options = {"max_cores": "quarter"}
    num_slices = compute_multiprocessing_slices(proc_options)
    assert num_slices == ncores // 4


def test_ramp_fit():
    ni, ng, nr, nc = 2, 1, 2, 2
    model, rnModel, gain = setup_inputs(nints=ni, ngroups=ng, nrows=nr, ncols=nc)
    print(DELIM)
    for val in range(ni):
        print(model.data)
        model.data[val, :, :, :] + (val+1)
        print(SDELIM)
        print(model.data)
        print(DELIM)

    return 0
    proc_options = {"max_cores": "none", "read_noise": rnModel, "gain": gain}
    print(DELIM)
    print(f"Original shape = {model.data.shape}")
    print(DELIM)
    ans = ramp_fit(model, proc_options)


# --------------------------------------------------------------------

if __name__ == "__main__":
    # test_get_row_slices()
    # test_compute_multiprocessing_slices()
    test_ramp_fit()
