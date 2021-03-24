#!/usr/bin/env python

import numpy as np


class IntermediateComputations:
    """
    The class contains all intermediate comptutation values to allow for more compact function
    prototypes, as well as an easy to use container from which to recall data across function
    calls.
    """
    def __init__(self, data):
        """
        Captures initial values needed for logging, then sets future values to None.
        """
        self.orig_dim = data.shape

        self.slope_est = None


    def set_stuff(self, stuff):
        """
        Stub for setting intermediate variables.
        """
        pass


class ModelArray:
    """
    Ramp fitting calls for looping over integrations, computing groups within an
    integration first, then compbining these computations.  This effectively
    reduces the dimensions of relevant arrays from 4-D to 3-D.  The RampModel
    schema expects arrays to be 4-D, so this intermediate model is used to allow
    for more flexibility.  The original RampModel data is carried around because
    it contains metadata relevant to computation, such as exposure data, instrument
    data, etc.
    """
    def __init__(self, orig_model, data, err, groupdq, pixeldq):
        """
        orig_model: RampModel
            This is the original ramp model containing the full data set.

        data: ndarray
            Separated data array to avoid model schema dimension restrictions.  This
            allows for dimension reduction that is not possible using schemas.

        err: ndarray
            Separated err array to avoid model schema dimension restrictions.  This
            allows for dimension reduction that is not possible using schemas.

        groupdq: ndarray
            Separated GROUPDQ array to avoid model schema dimension restrictions.  This
            allows for dimension reduction that is not possible using schemas.

        pixeldq: ndarray
            Separated PIXELDQ array to avoid model schema dimension restrictions.  This
            allows for dimension reduction that is not possible using schemas.

        """
        self.orig_model = orig_model
        self.data = data
        self.err = err
        self.groupdq = groupdq
        self.pixeldq = pixeldq
        self.max_seg_length = 0


class RampSegment:
    def __init__(self, start, end):
        self.start = start  # The first group in the segment
        self.end = end  # The end group in the segment
        self.length = end - start + 1  # Add one because it's inclusive
        self.slope = 0.  # The slope of the segment
        self.var_p = 0.  # The variance of the slope due to Poisson noise
        self.var_r = 0.  # The variance of the slope due to read noise
        self.var_c = 0.  # The combined variance of the slope due to noise
