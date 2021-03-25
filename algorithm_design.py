#!/usr/bin/env python

"""
Based on:
https://jwst-pipeline.readthedocs.io/en/latest/jwst/ramp_fitting/description.html
"""

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
from jwst.datamodels import GainModel
from jwst.datamodels import RampModel
from jwst.datamodels import ReadnoiseModel

DO_NOT_USE = dqflags.group["DO_NOT_USE"]
JUMP_DET = dqflags.group["JUMP_DET"]
SATURATED = dqflags.group["SATURATED"]

# INTERMEDIATE_DTYPE = np.64float
INTERMEDIATE_DTYPE = np.float32

DELIM = "-" * 70
SDELIM = "*" * 70


# --------------------------------------------------------------------
#                        Primary Function
# --------------------------------------------------------------------

"""
The primary output product is an ImageModel:
    A data model for 2D images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    area : numpy float32 array
         Pixel area map array

    pathloss_point : numpy float32 array
         Pathloss correction for point source

    pathloss_uniform : numpy float32 array
         Pathloss correction for uniform source


 ---> schema_url = "http://stsci.edu/schemas/jwst_datamodel/image.schema"

---------------------------------------------------------------------------
The rateint output product is a CubeModel
    A data model for 3D image cubes.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zero frame array

    area : numpy float32 array
         Pixel area map array

    int_times : numpy table
         table of times for each integration

    wavelength : numpy float32 array
         Wavelength array

    var_poisson : numpy float32 array
         Integration-specific variances of slope due to Poisson noise

    var_rnoise : numpy float32 array
         Integration-specific variances of slope due to read noise


---> schema_url = "http://stsci.edu/schemas/jwst_datamodel/cube.schema"

---------------------------------------------------------------------------
The optional results product is a RampFitOutputModel:
    A data model for the optional output of the ramp fitting step.

    In the parameter definitions below, `n_int` is the number of
    integrations, `max_seg` is the maximum number of segments that
    were fit, `nreads` is the number of reads in an integration, and
    `ny` and `nx` are the height and width of the image.

    Parameters
    __________

    slope : numpy float32 array (n_int, max_seg, ny, nx)
        Segment-specific slope

    sigslope : numpy float32 array (n_int, max_seg, ny, nx)
        Sigma for segment-specific slope

    var_poisson : numpy float32 array (n_int, max_seg, ny, nx)
        Variance due to poisson noise for segment-specific slope

    var_rnoise : numpy float32 array (n_int, max_seg, ny, nx)
        Variance due to read noise for segment-specific slope

    yint : numpy float32 array (n_int, max_seg, ny, nx)
        Segment-specific y-intercept

    sigyint : numpy float32 array (n_int, max_seg, ny, nx)
        Sigma for segment-specific y-intercept

    pedestal : numpy float32 array (n_int, max_seg, ny, nx)
        Pedestal array

    weights : numpy float32 array (n_int, max_seg, ny, nx)
        Weights for segment-specific fits

    crmag : numpy float32 array (n_int, max_seg, ny, nx)
        Approximate CR magnitudes


---> schema_url = "http://stsci.edu/schemas/jwst_datamodel/rampfitoutput.schema"
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
    primary:ImageModel
        The primary ramp fitting product.

    rateint: TODO (Need to check what type this is)
        Flesh out details

    optres: TODO
        Flesh out details
    """

    icomp = IntermediateComputations(model.data)
    miri_answer = miri_correction(model)
    # TODO check to make sure there is data left to continue

    # Create a model agnostic class that allows for dimension reduction in order to
    # separate across integration computations from within integration computations.
    model_4d = ModelArray(model, model.data, model.err, model.groupdq, model.pixeldq)

    primary, rateint, optres = ramp_fit_4d(model_4d, proc_options, icomp)

    return primary, rateint, optres


# --------------------------------------------------------------------
#                          Helper Functions
# --------------------------------------------------------------------
def across_integration_ramp_fitting(model, proc_options, icomp):
    """
    Fits each ramp in each integration, then computes fit over all
    integrations.

    Parameter
    ---------
    model: RampModel
        Contains the information needed to fit ramps.

    proc_options: dict
        Contains processing information

    Return
    ------
    primary:ImageModel
        The primary ramp fitting product.

    rateint: TODO (Need to check what type this is)
        Flesh out details

    optres: TODO
        Flesh out details
    """
    # This is needed to compute the variance due to Poisson noise and is computed
    # across all integrations.
    slope_est = compute_slope_est(model, icomp)

    primary_int = []
    rateint_int = []
    optres_int = []
    for integration in range(model.data.shape[0]):

        # Get 3-D cube for the integeration
        model_3d, icomp_3d = get_integration(model, icomp, integration)
        icomp_3d.slope_est = icomp.slope_est

        # Compute fit and variances on each cube, then save each integration
        # output to a list.
        primary, rateint, optres = ramp_fit_3d(model_3d, proc_options, icomp_3d)
        primary_int.append(primary)
        rateint_int.append(rateint)
        optres_int.append(optres)

    # Combine each integration to compute the overall exposure fit and variances.
    primary, rateint, optres = combine_integrations(
        primary_int, rateint_int, optres_int
    )

    return primary, rateint, optres


def combine_integrations(primary_int, rateint_int, optres_int):
    """
    primary_int: list
        Contains the primary product for each integration

    rateint_int: list
        Contains the rateint product for each integration

    optres_int: list
        Contains the optional results for each integration
    """
    """
    if len(primary_int) == 1:
        return primary_int[0], rateint_int[0], optres_int[0]
    """

    # TODO Combine integration level info into exposure level info
    # Create exposure level output models
    # Compute exposure level outputs from integration level computations
    # Populate exposure level outputs
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
        first_diffs = data.copy()
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
        (
            i_group,
            i_yy,
            i_xx,
        ) = np.where(np.bitwise_and(groupdq[1:, :, :], JUMP_DET))
        first_diffs[i_group - 1, i_yy, i_xx] = np.NaN

        del i_group, i_yy, i_xx

        # Check for pixels in which there is good data in 0th group, but
        #   all first_diffs for this ramp are NaN because there are too
        #   few good groups past the 0th. Due to the shortage of good
        #   data, the first_diffs will be set here equal to the data in
        #   the 0th group.
        wh_min = np.where(
            np.logical_and(
                np.isnan(first_diffs).all(axis=0), np.isfinite(data[0, :, :])
            )
        )
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
    model: ModelArray
        Contains 3-D cubes needed for ramp fitting.

    proc_options: dict
        Contains the processing options.

    icomp: IntermediateComuptations
        Contains intermediate computations.

    row: int
        Row of current pixel

    col: int
        Column of current pixel

    Return
    ------
    # TODO - change the return values to be
    #
    segments_fit: list

    slope_i: float

    var_ri: float

    var_pi: float

    var_ci: float
    """
    ramp = ModelArray(
        model.orig_model,
        model.data[:, row, col].astype(
            INTERMEDIATE_DTYPE
        ),  # astypes may be overkill here
        model.err[:, row, col].astype(INTERMEDIATE_DTYPE),
        model.groupdq[:, row, col],
        model.pixeldq[row, col],
    )

    # compute segments
    segment_fits = compute_segments_fits(ramp, proc_options, icomp, row, col)

    # Use segments to compute integration level ramp fitting
    var_ri, var_pi, var_ci, num_s, den_s = 0.0, 0.0, 0.0, 0.0, 0.0
    for seg_fit in segment_fits:
        slope_s, y_int, var_rs, var_ps, var_cs = seg_fit
        var_ri += 1.0 / var_rs
        var_pi += 1.0 / var_ps
        var_ci += 1.0 / (var_rs + var_ps)
        num_s += slope_s / var_cs
        den_s += 1.0 / var_cs

    # Complete integration level computations for the pixel
    var_ri = 1.0 / var_ri
    var_pi = 1.0 / var_pi
    var_ci = 1.0 / var_ci
    slope_i = num_s / den_s

    return segments_fit, slope_i, var_ri, var_pi, var_ci


def compute_segments_fits(ramp, proc_options, icomp, row, col):
    """
    Compute the segments in the ramp, then compute the fit of each ramp.

    Parameter
    ---------
    ramp: ModelArray
        Contains the ramp to be fitted.

    proc_options: dict
        Contains the processing option.

    icomp IntermediateComputation
        Contains the intermediate computations.

    row: int
        The row of the current ramp.

    col: int
        The column of the current ramp.

    Return
    ------
    segments_fits: list
        A list of tuples.  Each tuple corresponds to the fit or each
        segment in the ramp.
    """
    # HERE
    segments = get_segments(ramp)
    segments_fits = []
    """
    TODO: Add the following loop
    for seg in segments:
    """
    seg = segments[0]
    seg_fit = compute_segment_fit(seg, ramp, proc_options, icomp, row, col)

    segments_fits.append(seg_fit)

    return segments_fits


def compute_segment_fit(seg, ramp, proc_options, icomp, row, col):
    """
    Compute the fit for a segment in the ramp.

    Parameter
    ---------
    seg: tuple
        Contains the first and last group number in the ramp, as well as the length.

    proc_options: dict
        Contains the processing option.

    icomp IntermediateComputation
        Contains the intermediate computations.

    row: int
        The row of the current ramp.

    col: int
        The column of the current ramp.

    Return
    ------
    segment_fit: tuple
        The computed segment fit.  The tuple contains the slope estimate, y-intercept,
        variance of the slope due to read noise, variance of the slope due to the
        Poisson noise, and the combined variance due to noise.
    """
    # Check the fitting alorithm to use, either the OLS or GLS algorithm.
    # Use the key work "algorithm" in proc_options.
    slope_s, var_rs, var_ps, var_cs = None, None, None, None
    if seg[2] == 1:
        slope_s, y_int = compute_segment_slope_len_1(seg, ramp, proc_options)
        var_rs = compute_segment_rnoise_len_1(seg, ramp, proc_options, row, col)
        var_ps = compute_segment_pnoise_len_1(seg, ramp, proc_options, icomp, row, col)
        var_cs = var_rs + var_ps

    elif seg[2] == 2:
        slope_s, y_int = compute_segment_slope_len_2(seg, ramp, proc_options)
        var_rs = compute_segment_rnoise_len_2(seg, ramp, proc_options, row, col)
        var_ps = compute_segment_pnoise_len_2(seg, ramp, proc_options, icomp, row, col)
        var_cs = var_rs + var_ps

    elif seg[2] > 2:
        slope_s, y_int = compute_segment_slope_len_3(seg, ramp, proc_options)
        var_rs = compute_segment_rnoise_len_3(seg, ramp, proc_options, row, col)
        var_ps = compute_segment_pnoise_len_3(seg, ramp, proc_options, icomp, row, col)
        var_cs = var_rs + var_ps

    segment_fit = (slope_s, y_int, var_rs, var_ps, var_cs)
    return segment_fit


def compute_segment_slope_len_1(seg, ramp, proc_options):
    """
    Compute slope with only one group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    Return
    ------
    slope
    y_int
    """
    return None


def compute_segment_rnoise_len_1(seg, ramp, proc_options, row, col):
    """
    Compute readnoise variance with only one group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    row: int
        The row of the ramp.

    col: int
        The column of the ramp.

    Return
    ------
    var_rs
    """
    return None


def compute_segment_pnoise_len_1(seg, ramp, proc_options, icomp, row, col):
    """
    Compute Poisson variance with only one group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    row: int
        The row of the ramp.

    col: int
        The column of the ramp.

    Return
    ------
    var_ps
    """
    return None


def compute_segment_slope_len_2(seg, ramp, proc_options):
    """
    Compute slope with only two group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    Return
    ------
    slope
    y_int
    """
    return None


def compute_segment_rnoise_len_2(seg, ramp, proc_options, row, col):
    """
    Compute read noise variance with only two group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    row: int
        The row of the ramp.

    col: int
        The column of the ramp.

    Return
    ------
    var_rs
    """
    return None


def compute_segment_pnoise_len_2(seg, ramp, proc_options, icomp, row, col):
    """
    Compute Poisson variance with only two group.

    Parameter
    ---------
    seg: tuple
        Contains the segment information - (first group, last group, segment length).

    ramp: ArrayModel
        Contains the ramp in which the segment lies.

    proc_options: dict
        Processing options.

    row: int
        The row of the ramp.

    col: int
        The column of the ramp.

    Return
    ------
    var_ps
    """
    return None


def compute_segment_pnoise_len_3(seg, ramp, proc_options, icomp, row, col):
    """
    Compute Poisson variance for segments with 3 or more groups.

    Parameter
    ---------
    seg
    ramp
    proc_options
    icomp
    row
    col

    Return
    ------
    var_ps
    """
    ng = seg[2]
    gain = proc_options["gain"].data[row, col]
    group_time = ramp.orig_model.meta.exposure.group_time

    var_ps = icomp.slope_est[row, col] / (group_time * gain * (ng - 1))

    return var_ps


def compute_segment_rnoise_len_3(seg, ramp, proc_options, row, col):
    """
    Compute readnoise variance for segments with 3 or more groups.

    Parameter
    ---------
    seg
    ramp
    proc_options
    row
    col

    Return
    ------
    var_rs
    """
    ng = seg[2]  # Number of groups in segment
    R = proc_options["readnoise"].data[row, col]  # Read noise for pixel
    group_time = ramp.orig_model.meta.exposure.group_time

    var_rs = 12 * R ** 2 / ((ng ** 3 - ng) * group_time ** 2)

    return var_rs


def compute_segment_slope_len_3(seg, ramp, proc_options):
    """
    Fit slope to segment of length 3 or greater.

    Parameter
    ---------
    seg
    ramp
    proc_options

    Return
    ------
    slope
    y_int
    """
    ng = seg[2]
    group_time = ramp.orig_model.meta.exposure.group_time
    X = np.array([k + 1 for k in range(ng)])
    X *= group_time
    Y = ramp.data[seg[0] : seg[1] + 1]

    sum_x = X.sum()
    sum_xx = (X * X).sum()

    sum_y = Y.sum()
    sum_yy = (Y * Y).sum()

    sum_xy = (X * Y).sum()

    # Compute the slope and y-intercept
    cov_xy = sum_xy - (sum_x * sum_y) / ng
    var_x = sum_xx - (sum_x ** 2) / ng

    slope = cov_xy / var_x
    y_int = sum_y / ng - slope * sum_x / ng

    # TODO
    # Do weighted least squares.

    return slope, y_int


def compute_slope_est(model, icomp):
    """
    Estimate the slope using the median group differences over each integration.

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
    # Get 4-D dimensions
    number_of_integrations, ngroups, nrows, ncols = model.data.shape

    # Compute slope estimates over all groups and integrations.  This is a more
    # robust measure to use when computing the Poisson variance in the slope.
    slope_est = np.zeros((nrows, ncols), dtype=INTERMEDIATE_DTYPE)
    icomp.slope_est = slope_est

    # Loop over each integration to do all computations on the 3-D cubes independently.
    for integration in range(number_of_integrations):

        # Get 3-D information from the
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
        nan_med[np.isnan(nan_med)] = 0.0  # if all first_diffs_sect are nans
        slope_est[:, :] += nan_med

        del nan_med
        del first_diffs

    # Save estimated slopes
    slope_est /= number_of_integrations
    slope_est /= model.orig_model.meta.exposure.group_time

    return slope_est


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
    data = model.data[n_int, :, :, :].astype(INTERMEDIATE_DTYPE)
    err = model.err[n_int, :, :, :].astype(INTERMEDIATE_DTYPE)
    groupdq = model.groupdq[n_int, :, :, :]
    pixeldq = model.pixeldq

    model_3d = ModelArray(model.orig_model, data, err, groupdq, pixeldq)
    icomp_3d = IntermediateComputations(model_3d.data)

    return model_3d, icomp_3d


def get_row_slices(nrows, num_slices):
    """
    Create a list with length equal to the number of slices.  Each element in the
    list is the number of rows to operate on in the original data cube.


    Parameter
    ---------
    nrows: int
        The number of rows in the image dimensions.

    num_slices: int
        The number of cores to use for multiprocessing

    Return
    ------
    row_slices: list
        The list of slices.  The length of the list is equal to the number of slices.
        The elements of the list are the number of rows that slice is to process.
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
    max_seg_length = 0
    segments = []
    # Examine ramp.groupdq to determine where the segment boundaries are
    # Get next segment
    # If segment length is 1 and max_seg_length > 1 discard segment
    # If segment length is greater than max_seg_length update max_seg_length
    # If not a discarded length 1 segment, append segment to segments

    # TODO Need to add other cases
    #      Right now assume the whole ramp is a valid segment.
    return [(0, len(ramp.data) - 1, len(ramp.data))]


def miri_correction(model):
    """
    Corrects data accounting for MIRI artifacts.  There is no return value.
    This function modifies the shape of the model.

    Parameter
    ---------
    model: RampModel
    """
    return None


def ramp_fit_3d(model, proc_options, icomp):
    """
    Rampfitting a 3-D cube.

    Parameter
    ---------
    model: ModelArray
        Contains the 3-D cube arrays and the metadata needed for processing.

    proc_options: dict
        Processing options

    icomp: IntermediateComputations
        The intermediate computations.

    Return
    ------
    primary: ImageModel

    rateint:

    opt_res:
    """
    ngroups, nrows, ncols = model.data.shape
    # TODO Create output models
    for row in range(nrows):
        for col in range(ncols):
            print(f"\n{DELIM}")
            print(f"    -----> ({row}, {col})")
            ans = compute_pixel_ramp(model, proc_options, icomp, row, col)
            # TODO unpack 'ans' and populate output models

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


def test_ramp_fit_simple_linear_ramp():

    dims = (1, 5, 1, 1)
    model, rnModel, gain = setup_inputs_simplified(dims=dims)

    rlist = [k + 1 for k in range(dims[1])]
    ramp = np.array(rlist)

    model.data[0, :, 0, 0] = ramp * 5

    proc_options = {
        "max_cores": "none",
        "opt_save": True,
        "gain": gain,
        "readnoise": rnModel,
    }
    print(DELIM)
    primary, rateint, optres = ramp_fit(model, proc_options)
    print(DELIM)


def test_ramp_fit_small():
    dims = (2, 5, 2, 2)
    model, rnModel, gain = setup_inputs_simplified(dims=dims)

    rlist = [k + 1 for k in range(dims[1])]
    ramp = np.array(rlist)
    model.data[0, :, 0, 0] = ramp * 5
    model.data[0, :, 1, 0] = ramp * 2
    model.data[0, :, 0, 1] = ramp * 7
    model.data[0, :, 1, 1] = ramp * 1

    model.data[1, :, 0, 0] = ramp * 5 + 2
    model.data[1, :, 1, 0] = ramp * 2 + 2
    model.data[1, :, 0, 1] = ramp * 7 + 2
    model.data[1, :, 1, 1] = ramp * 1 + 2

    proc_options = {
        "max_cores": "none",
        "opt_save": True,
        "gain": gain,
        "readnoise": rnModel,
    }
    print(DELIM)
    primary, rateint, optres = ramp_fit(model, proc_options)
    print(DELIM)


# --------------------------------------------------------------------

if __name__ == "__main__":
    # test_ramp_fit_simple_linear_ramp()
    test_ramp_fit_small()
