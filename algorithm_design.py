#!/usr/bin/env python

import multiprocessing
import sys
import warnings

import numpy as np

from helper_funcs import setup_inputs
from helper_funcs import setup_inputs_simplified

from algorithm_design_classes import IntermediateComputations
from algorithm_design_classes import ModelArray
from algorithm_design_classes import RampSegment

from jwst.datamodels import dqflags
from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]

DELIM = "-" * 70
SDELIM = "*" * 70

"""
Receive 4-D RampModel with processing information:
"""

def ramp_fit(model, proc_options):
    """
    Top level ramp fitting designed for 4-D processing  with dimensions:
        (integrations, groups, rows, columns)

    Parameter
    ---------
    model: RampModel
        Contains the information needed to fit ramps.

    proc_options: dict
        Contains processing information

    Return
    ------
    TODO
    primary: RampModel (or ImageModel <--- will have to check this)
        The primary ramp fitting product.

    rateint: (Need to check what type this is)
        Flesh out details
    
    optres:
        Flesh out details
    """
    icomp = IntermediateComputations(model.data)
    miri_answer = miri_correction(model)

    # Create a model agnostic class that allows for dimension reduction in order to
    # separate across integration computations from within integration computations.
    model_4d = ModelArray(model, model.data, model.err, model.groupdq, model.pixeldq)

    return ramp_fit_4d(model_4d, proc_options, icomp)


# --------------------------------------------------------------------

def across_integration_ramp_fitting(model, proc_options, icomp):
    """
    Fits each ramp in each integration, then computes fit over all
    integrations.
    """
    slope_est = compute_slope_est(model, icomp)

    primary_int = []
    rateint_int = []
    optres_int = []
    for integration in range(model.data.shape[0]):
        model_3d, icomp_3d = get_integration(model, icomp, integration)
        primary, rateint, optres = ramp_fit_3d(model_3d, proc_options, icomp_3d)
        primary_int.append(primary)
        rateint_int.append(rateint)
        optres_int.append(optres)

    return combine_integrations(primary_int, rateint_int, optres_int)


def combine_integrations(primary_int, rateint_int, optres_int):
    return None, None, None


def compute_first_diffs(data, groupdq):
    """
    Compute the median first differences between groups not affected by
    saturation or cosmic arrays.

    Parameter
    ---------
    data: ndarray
        The 3-D array containing the pixel data for groups.

    groupdq: ndarray
        The 3-D array containing the data quality for each pixel in each
        group.

    Return
    ------
    first_diffs: ndarray
        The median of the first differences to estimate the slope of the ramp.
    """
    # Compute the first differences of all groups
    first_diffs = np.diff(data, axis=0)

    # If the dataset has only 1 group/integ, assume the 'previous group'
    #   is all zeros, so just use data as the difference
    if first_diffs.shape[0] == 0:
        first_diffs= data.copy()
    else:
        # Similarly, for datasets having >1 group/integ and having
        #   single-group segments, just use the data as the difference
        wh_nan = np.where(np.isnan(first_diffs[0, :, :]))

        if len(wh_nan[0]) > 0:
            first_diffs[0, :, :][wh_nan] = data[0, :, :][wh_nan]

        del wh_nan

        # Mask all the first differences that are affected by a CR,
        #   starting at group 1.  The purpose of starting at index 1 is
        #   to shift all the indices down by 1, so they line up with the
        #   indices in first_diffs.
        i_group, i_yy, i_xx, = np.where(np.bitwise_and(groupdq[1:, :, :], JUMP_DET))
        first_diffs[i_group - 1, i_yy, i_xx] = np.NaN

        del i_group, i_yy, i_xx

        # Check for pixels in which there is good data in 0th group, but
        #   all first_diffs for this ramp are NaN because there are too
        #   few good groups past the 0th. Due to the shortage of good
        #   data, the first_diffs will be set here equal to the data in
        #   the 0th group.
        wh_min = np.where(np.logical_and(
            np.isnan(first_diffs).all(axis=0), np.isfinite(data[0, :, :])))
        if len(wh_min[0] > 0):
            first_diffs[0, :, :][wh_min] = data[0, :, :][wh_min]

        del wh_min

    return first_diffs


def compute_multiprocessing_slices(proc_options):
    """
    Compute multiprocessing slices per processor.

    Parameter
    ---------
    proc_options: dict
        Dictionary must contain the key word 'max_cores' which is a string.  The
        string value will determine the portion of available cores to use for
        processing, "all", "half", or "quarter".  Other values default to single
        core processing.  By convention "none" denotes single core processing.

    Return
    ------
    number_slices: int
        The number of cores to use for multiprocessing
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


def compute_pixel_ramp(model, proc_options, icomp, row, col):
    """

    Parameter
    ---------
    model
    proc_options
    icomp
    row
    col

    Return
    ------
    """
    # HERE
    ramp = ModelArray(
        model, model.data[:, row, col], model.err[:, row, col], 
        model.groupdq[:, row, col], model.pixeldq[row, col])

    # compute segments
    ans = compute_segments_fits(ramp, proc_options, icomp, row, col)

    # compute overall ramps
    ans = pixel_fit(ramp, ans, proc_options, icomp, row, col)

    sys.exit(1)


def compute_segments_fits(ramp, proc_options, icomp, row, col):
    segments = get_segments(ramp)
    sys.exit(1)


def compute_slope_est(model, icomp):
    """
    Estimate the slope using the median differences over each integration.

    Parameter
    ---------
    model: ArrayModel
        Contains the full 4-D image data.

    Return
    ------
    slope_est: ndarray
        The slope estimates using the median differences between groups
        not affected by across all integrations.
    """
    number_of_integrations, ngroups, nrows, ncols = model.data.shape[0]
    slope_est = np.zeros((nrows, ncols))
    for integration in range(number_of_integrations):
        data = model.data[integration, :, :, :]
        groupdq = model.groupdq[integration, :, :, :]

        # Reset all saturated groups in the input data array to NaN
        where_sat = np.where(np.bitwise_and(groupdq, SATURATED))
        data[where_sat] = np.NaN
        del where_sat

        first_diffs = compute_first_diffs(data, groupdq)

        # All first differences affected by saturation and CRs have been set
        #  to NaN, so compute the median of all non-NaN first differences.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "All-NaN.*", RuntimeWarning)
            nan_med = np.nanmedian(first_diffs, axis=0)
        nan_med[np.isnan(nan_med)] = 0.  # if all first_diffs_sect are nans
        slope_est[:, :] += nan_med
        del nan_med

    # Save estimated slopes
    slope_est /= number_of_integrations
    slope_est /= model.orig_model.meta.exposure.group_time


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


def get_row_slices(nrows, num_slices):
    """
    Create a list with length equal to the number of slices.  Each element in the list is
    the number of rows to operate on in the original data cube.


    Parameter
    ---------
    nrows: int
        The number of rows in the image dimensions.

    num_slices: int
        The number of cores to use for multiprocessing

    Return
    ------
    row_slices: list
        The list of slices.  The length of the list is equal to the number of slices.  The
        elements of the list are the number of rows that slice is to process.
    """
    rows_per_slice = nrows // num_slices
    remainder = nrows % num_slices

    rslices = [rows_per_slice] * num_slices
    row_slices = np.array(rslices, dtype=np.int32)

    # Spread the remainder across the slices, so the difference between the number of
    # rows each core processes is no greater than 1.  For example, if an image has 16
    # and with 3 cores, the slice list will be [6, 5, 5]
    if remainder > 0:
        row_slices[:remainder] += 1

    return row_slices


def get_segments(ramp):
    """
    Need to figure this out.  Are groups labelled as JUMP_DET used or separated?
    How to break up a set of groups into segments?
    """
    pass


def miri_correction(model):
    return None


def pixel_fit(ramp, ans, proc_options, icomp, row, col):
    """
    Fit segments, then fit use segments to fit overall.
    """
    pass


def ramp_fit_3d(model, proc_options, icomp):
    """
    Rampfitting a 3-D cube.

    model: ModelArray
        Contains the 3-D cube arrays and the metadata needed for processing.

    proc_options: dict
        Processing options

    icomp: IntermediateComputations
        The intermediate computations.
    """
    # TODO
    '''
    ngroups, nrows, ncols = model.data.shape
    for row in range(nrows):
        for col in range(ncols):
            ans = compute_pixel_ramp(model, proc_options, icomp, row, col)
    '''
    ans = compute_pixel_ramp(model, proc_options, icomp, 0, 0)

    return None, None, None
    # return primary, rateint, opt_res


def ramp_fit_4d(model, proc_options, icomp):
    """
    Rampfitting a 4-D cube.  The data in the model is the 4-D data after the MIRI correction
    (if any).  The 4-D computations first perform within integration computations, then
    combine those computations to compute across integrations getting the final exposure
    level computations.  These computations can be done on a single process or in parallel.

    Parameter
    ---------
    model: RampModel
        Contains the information needed to fit ramps.

    proc_options: dict
        Contains processing information

    Return
    ------
    TODO
    primary: RampModel (or ImageModel <--- will have to check this)
        The primary ramp fitting product.

    rateint: (Need to check what type this is)
        Flesh out details
    
    optres:
        Flesh out details
    """
    # print("    *** ramp_fit_4d")
    num_slices = compute_multiprocessing_slices(proc_options)
    if num_slices == 1:
        return across_integration_ramp_fitting(model, proc_options, icomp)
    else:
        return ramp_fit_4d_multi(model, proc_options, icomp, num_slices)


def ramp_fit_4d_multi(model, proc_options, icomp, num_slices):
    """
    TODO: Needs implementation.  Break down the data cube by rows in the image dimension.
    """
    print("    *** ramp_fit_4d_multi")


# NEXT FUNCTION


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
    proc_options = {"max_cores":  "none"}

    dims = (2, 5, 1, 2)
    model, rnModel, gain = setup_inputs_simplified(dims=dims)

    print(DELIM)
    print(f"Dimensions = {dims}");
    ramp = np.array([k+1 for k in range(dims[1])])

    # Integration 0
    # model.data[0, :, 0, 0] = ramp * 5 + 5
    model.data[0, :, 0, 0] = np.array([10. ,15., 65., 70., 75.])
    model.groupdq[0, 2, 0, 0] = JUMP_DET
    model.data[0, :, 0, 1] = ramp * 8 + 3

    # Integration 1
    model.data[1, :, 0, 0] = ramp * 15 + 1
    model.data[1, :, 0, 1] = ramp * 18 + 9
    print(DELIM)

    ans = ramp_fit(model, proc_options)


# --------------------------------------------------------------------

if __name__ == "__main__":
    test_ramp_fit()
