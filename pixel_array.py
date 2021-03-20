#!/usr/bin/env python

import numpy as np

import ramp_fitting_classes

from ramp_fitting_classes import RampInput4d
from ramp_fitting_classes import RampOutput4d
from ramp_fitting_classes import RampInput3d
from ramp_fitting_classes import RampOutput3d


DELIM = "-" * 70

"""
https://jwst-pipeline.readthedocs.io/en/latest/jwst/ramp_fitting/description.html
                SPECIAL CASES
1. If for each n_int ngroups == 1:
       if pixel is not saturated:
           rate_count = pixel / group_time

3. If for each n_int ngroups == 2:
       if pixel is not saturated in each group:
           rate_count = (pixel[1] - pixel[0]) / group_time

3. If for each n_int ngroups > 1:
       1. Segment has ngroups == 1
              if nsegs == 1:
                  rate_count = pixel / group_time
              if nsegs > 1 and ngroups > 1 in other segment:
                  segement with length 1 is ignored

       2. If pixel[0] is valid and pixel[k] is saturated for all k > 0:
              first_diff = pixel[0]

       3. If pixel[1] is saturated:
              rate_count = pixel[0] / group_time

4. MIRI correction step removes pixel[0] and pixel[-1], so the actual
   ngroups = orig_ngroups - 2.
       if ngroups < 2:
           log warning



                PRODUCTS
    Primary product:
1. 2-D: SCI - pixel slopes computed using weighted segment slopes and integrations.
2. 2-D: GROUPDQ - bitwise "or" over integrations and groups.
3. 2-D: The 3-D VAR_POISSON averaged over all integrations to 2-D.
4. 2-D: The 3-D VAR_RNOISE averaged over all integrations to 2-D.

    rateints product (per integration products):
1. 3-D: SCI - pixel slopes computed using weighted segment slopes.
2. 2-D: GROUPDQ - bitwise "or" over groups and integrations.
3. 3-D: The 4-D VAR_POISSON averaged over all fit segments in an integration to 3-D.
4. 3-D: The 4-D VAR_RNOISE averaged over all fit segments in an integration to 3-D.
5. 3-D: ERR - 

    Optional output (from save_opt flag); data from each segement for each pixel:
    * should "segment" above read "group"?
1. 4-D: SLOPE - slopes
2. 4-D: SIGSLOPE - uncertainty in the slopes
3. 4-D: YINT - y-intecepts
4. 4-D: SIGYINT - uncertainty in the y-intecepts
5. 4-D: WEIGHTS - fitting weights
6. 4-D: VAR_POISSON - variance of the slope due to poisson noise
7. 4-D: VAR_NOISE - variance of the slope due to read noise
8. 3-D: PEDESTAL - the signal at zero exposure time for each pixel
            - Computed using fitted data, setting exposure time to 0. 
              (*note: How is this different from y-intercepts?)
            - Any saturated pixel, this value is 0. 
9. 4-D: CRMAG - the magnitude of each group flagged as having a cosmic ray hit


                SLOPE AND VARIANCE CALCULATIONS
Notation:
    R - read noise
    P - Poisson noise
    C - combined noise
    S - segment
    I - integration
    O - overall


S (SNR) = (data * gain) / sqrt(rnoise**2 + (data * gain)**2)
w[k] (weighting) = w[k] = (k - k[midpoint])**P (*note: I don't understand this)
    - P = exponent of the weight determined by S

S ranges            P (weight exponent)
---------------------------------------
[0, 5)              0
[5, 10)             0.4
[10, 20)            1
[20, 50)            3
[50, 100)           6
[100, Inf)          10


**********************************************************************
                        SEGMENT CALCULATIONS
**********************************************************************
    ---- Read noise variance of the slope of a segment:

var[R, S] = 12 * R**2 / ((ngroups[S]**3 - ngroups[S]) * tgroups**2)

    - R = the noise in the difference between 2 frames (*note: frame or group?)
    - S = segment
    - ngroups[S] = the number of groups in the segment
    - tgroup = group time in seconds


    ---- Poisson noise variance of the slope of a segment:

var[P, S] = est_slope / (trgroup * gain * (ngroups[S] - 1))

    - P = the noise in the difference between 2 frames (*note: frame or group?)
    - S = segment
    - gain is the gain for the pixel (from reference file)
    - est_slope = median(first differences unaffected by saturation and cosmic 
                         rays in integrations)
    - ngroups[S] = the number of groups in the segment
    - tgroup = group time in seconds


    ---- Combined variance of the slope of a segment:

var[C, S] = var[R, S] + var[P, S]


**********************************************************************
                        INTEGRATION CALCULATIONS
**********************************************************************
    ---- Read noise variance of a slope over an integration:

var[R, I] = 1 / (sum(over S){1 / var[R, S]})


    ---- Poisson noise variance of a slope over an integration:

var[P, I] = 1 / (sum(over S){1 / var[P, S]})


    ---- Combined variance of a slope over an integration:

var[C, I] = 1 / (sum(over S){1 / (var[R, S] + var[P, S])})


    ---- Slope of an integration:

N = sum(over S){slope[S] / var[C, S]}
D = sum(over S){1 / var[C, S]}
slope[I] = N / D


**********************************************************************
                        EXPOSURE CALCULATIONS
**********************************************************************
    ---- Read noise variance over all integrations:

var[R, O] = 1 / sum(over I){1 / var[R, I]}


    ---- Poisson noise variance over all integrations:

var[P, O] = 1 / sum(over I){1 / var[P, I]}


    ---- Combined variance of the slope over all integrations:

var[C, O] = var[R, O] + var[P, O]


    ---- Slope over all integrations:
This is ambiguously defined in the documentation, since no [I, S] index is
defined for slopes and variances.  I think the I index is implied above, though.

N = sum(over S, over I){slope[I, S] / var[C, I, S]}
D = sum(over S, over I){1 / var[C, I, S]}
slope[I] = N / D


**********************************************************************
                        ERROR PROPOGATION
**********************************************************************
        Integration level:
ERR - the square root of the integration-level combined variances
VAR_POISSON
VAR_RNOISE

        Exposure level:
ERR - the square root of the exposure-level combined variances
VAR_POISSON
VAR_RNOISE
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
