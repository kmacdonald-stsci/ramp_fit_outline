#!/usr/bin/env python

import numpy as np

import ramp_fitting_classes

from ramp_fitting_classes import RampInput4d
from ramp_fitting_classes import RampOutput4d
from ramp_fitting_classes import RampInput3d
from ramp_fitting_classes import RampOutput3d


DELIM = "-" * 70


def ramp_fitting_4d(input_data):
    """
    Computes ramp fitting for a 4-D data cube.  It computes the ramp fit
    for each pixel in each integration then computes the total ramp fit
    across all integrations.

    Parameter
    ---------
    input_data: RampInput4d
        Contains all the arrays and parameters needed to compute ramp fitting.  This
        contains a 4-D image cube with dimensions (integrations, groups, rows, columns).

    Return
    ------
    ramp_out_4d: RampOutput4d
        Contains all the information for the computed ramp fit of the 4-D cube.
    """
    integrations = input_data.data_4d.shape[0]
    integration_data = []
    for integration in integrations:
        input_cube = RampInput3d(None, input_data, integration)
        return_cube_data = ramp_fitting_3d(input_cube)
        integration_data.append(return_data)

    ramp_out_4d = compute_across_integrations(input_data, integration_data)
    return ramp_out_4d


def ramp_fitting_3d(cube_data):
    """
    Compute ramp fit for each pixel over the groups.  This is done taking into account
    saturation and cosmic rays.  A pixel segment over a group is defined as a set of
    contiguous groups between cosmic rays.  Ramp fitting for a pixel first computes a
    ramp fit for each segment, then using the data for each segment, computes a ramp
    fit for the entire group.

    Parameter
    ---------
    cube_data: RampInput3d
        Contains all the arrays and parameters needed to compute ramp fitting.  This
        contains a 3-D image cube with dimensions (groups, rows, columns).  For ramp
        fitting done for a 4-D image cube with dimensions (groups, rows, columns), this
        variable is the data pertaining to a particular integration.

    Return
    ------
    ramp_out_3d: RampOutput3d
        Contains all the computations of ramp fitting of all pixels acros all groups.
    """
    # Maybe take advantage of
    if cube_data.multiprocessing != "None":
        raise ValueError("Ramp fitting multiprocessing not implemented, yet.")

        number_slices = compute_number_slices(max_cores)
        rows_per_slice = cube_data.data_3d.shape[1] / (number_slices - 1)

        # Break the cube up by rows and multiprocess row sets.
        slices = [
            (k * rows_per_slice, (k + 1) * row_per_slice)
            for k in range(number_slices - 1)
        ]
        last_slice = slices[-1]
        last_row = last_slice[0]
        last_rows = cube_data.data_3d.shape[1]
        slices.append((last_row, last_row + last_rows))

        return None  # Use below when implmemented
        # ramp_out_3d = ramp_fitting_pixels_multi(cube_data, slices)
        # return ramp_out_3d

    # Process all rows on a single process.
    ramp_out_3d = ramp_fitting_pixels(cube_data, 0, cube_data.data_3d.shape[1])
    return ramp_out_3d


def ramp_fitting_pixels(cube_data, row_start, row_end):
    """
    Computes ramp fitting on pixels in a set of rows.  The rows processed are all rows
    such that row_start <= row < row_end.

    Parameter
    ---------
    cube_data: RampInput3d
        Contains all the arrays and parameters needed to compute ramp fitting.  This
        contains a 3-D image cube with dimensions (groups, rows, columns).  For ramp
        fitting done for a 4-D image cube with dimensions (groups, rows, columns), this
        variable is the data pertaining to a particular integration.
    row_start: int
        The first row to be processed.
    row_end: int
        The last row exclusive to be processed.  The

    """
    output = RampOutput3d()
    # Create output data structure/class
    for row in range(row_start, row_end, 1):
        # Process each pixel in a row
        for col in cube_data.data_3d.shape[2]:
            # Proces each pixel
            ret = ramp_fitting_pixels(output_3d, cube_data, row, col)
    return output


def ramp_fitting_pixel(output_3d, cube_data, row, col):
    pixel = cube_data.data_3d[:, row, col]
    # Find segments
    # Compute least squares for each segment
    # Compute overall least squares for total group
    # Compute variances
    return True


def compute_across_integrations(input_data, integration_output):
    """
    Using the computed least squares for each pixel in each integration, compute the
    total least squares across the integrations.

    Parameter
    ---------
    input_data: RampInput4d
        Contains all the arrays and parameters needed to compute ramp fitting.  This
        contains a 4-D image cube with dimensions (integrations, groups, rows, columns).

    integration_output: list of RampOutput4d

    Return
    ------
    output4d: RampOutput4d
    """
    output4d = RampOutput4d()
    # For each pixel, compute the average ramp fit across each integration
    """
    for row input_data.data_4d.shape[2]:
        for col input_data.data_4d.shape[3]:
            pixel_ramp = sum(integration_output[:, row, col]) / 3.0
    """
    return output4d


def compute_number_slices(max_cores):
    # Determine number of slices to use for multi-processor computations
    if max_cores == "none":
        number_slices = 1
    else:
        num_cores = multiprocessing.cpu_count()
        log.debug(f"Found {num_cores} possible cores to use for ramp fitting")
        if max_cores == "quarter":
            number_slices = num_cores // 4 or 1
        elif max_cores == "half":
            number_slices = num_cores // 2 or 1
        elif max_cores == "all":
            number_slices = num_cores
        else:
            number_slices = 1
    return number_slices
