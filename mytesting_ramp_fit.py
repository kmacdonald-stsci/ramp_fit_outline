#!/usr/bin/env python

import numpy as np

import ramp_fitting
import ramp_fitting_classes

from ramp_fitting_classes import RampInput4d
from ramp_fitting_classes import RampOutput4d
from ramp_fitting_classes import RampInput3d
from ramp_fitting_classes import RampOutput3d


# --------------------------------------------------------------------


def generate_4d_image():
    """
    Generate a random 4-D image cube with integer values 0<= x < 65,535.
    """
    rid = RampInput4d()
    integrations, groups, rows, cols = 2, 5, 5, 5
    data = np.random.rand(integrations, groups, rows, cols)

    int_range = 65_536  # Make it look like real data
    data = data * int_range
    data = np.uint16(data)

    return data


def generate_3d_image():
    """
    Generate a random 3-D image cube with integer values 0<= x < 65,535.
    """
    rid = RampInput4d()
    integrations, groups, rows, cols = 1, 5, 5, 5
    data_4d = np.random.rand(integrations, groups, rows, cols)
    data = data_4d[0, :, :, :]

    int_range = 65_536  # Make it look like real data
    data = data * int_range
    data = np.uint16(data)
    return data


def basic_testing():
    rid = RampInput4d()
    rid.set_save_option(False)
    print(repr(rid))


# --------------------------------------------------------------------

if __name__ == "__main__":
    print(generate_3d_image())
