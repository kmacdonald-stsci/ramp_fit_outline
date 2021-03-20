#!/usr/bin/env python

import numpy as np

import ramp_fitting_classes

from ramp_fitting_classes import RampInput4d
from ramp_fitting_classes import RampOutput4d
from ramp_fitting_classes import RampInput3d
from ramp_fitting_classes import RampOutput3d


DELIM = "-" * 70

"""
"""

def segment_slope(pixel_seg, group_time):
    """
    Compute the group slope of a pixel. 
    """
    if len(pixel_seg) == 1:
        slope = pixel_seg[0] / group_time
        return slope
    elif len(pixel_seg) == 2:
        slope = (pixel_seg[1] - pixel_seg[0]) / group_time
        return slope

    X = group_time_seg(group_time, len(pixel_seg))
    Y = pixel_seg
    sx, sy, xx, yy, xy = 0, 0, 0, 0, 0
    n = len(X)
    for k in range(n):
        sx += X[k]          # Sum of X's
        sy += Y[k]          # Sum of Y's
        xx += X[k] * X[k]   # Sum of X squares
        yy += Y[k] * Y[k]   # Sum of Y squares
        xy += X[k] * Y[k]   # Sum of XY's

    slope = (n * xy - sx * sy) / (n * xx - sx * sx)
    y_int = (sy - slope * sx) / n

    return slope, y_int


def group_time_seg(group_time, length):
    x = np.arange(length) + 1
    return x * group_time


def three_segs():
    length = 15
    x = np.arange(length) + 1

    m1, m2, m3 = 5, 10, 20
    s1 = length // 3
    s2 = s1 + s1
    x[:s1] *= m1
    x[s1:s2] *= m
    x[s2:] *= m3
    return x

def two_segs():
    length = 10
    x = np.arange(length) + 1

    m1, m2 = 5, 10
    mid = length // 2
    x[:mid] *= m1
    x[mid:] *= m2
    return x


def simple_slope():
    m = 5
    x = (np.arange(10) + 1) * m
    return x


def main():
    x = simple_slope()
    slope, y_int = segment_slope(x, 1.)
    print(DELIM)
    print(x)
    print(DELIM)
    print(f"slope = {slope}, y_int = {y_int}")
    print(DELIM)

if __name__=="__main__":
    main()
