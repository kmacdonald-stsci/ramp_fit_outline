#!/usr/bin/env python


import numpy as np


DELIM = '-' * 70


def random_matrix_u16(integrations=2, groups=5, rows=4, columns=4):
    """

    Parameter
    ---------
    integrations: int
        The number of integrations

    groups: int
        The number of groups

    rows: int
        The number of rows

    columns: int
        The number of columns

    Return
    ------
    cube: ndarray
        A 4-D cube noisy ramp.
    """
    mat = np.random.rand(integrations, groups, rows, columns)
    u16 = 2**16
    mat = np.uint16(mat * u16)
    return mat


def random_image_slopes(nbits=8, rows=4, columns=4):
    """

    Parameter
    ---------
    nbits: int
        The maximum size, in bits, of the slopes.

    rows: int
        The number of rows

    columns: int
        The number of columns

    Return
    ------
    mat: ndarray
        A 2-D array of pixel slopes.
    """
    mat = np.random.rand(rows, columns) # random array of [0, 1)
    limit = 2**nbits
    mat = np.uint16(mat * limit)
    return mat


def random_noisy_ramp(slope_array=None, nreads=5, noise_bits=4):
    """
    Given ramp slopes for an image, create noisy 3-D ramp data.

    Parameter
    ---------
    slope_array: ndarray
        Slopes for each pixel


    nreads: int
        Number of reads

    noise_bits: int
        Desired level of noise (in bits, i.e., 4 bits gives random noise 
        from 0 to 15.

    Return
    ------
    cube: ndarray
        A 3-D cube noisy ramp.
    """
    if slope_array is None:
        return None
    rows, cols = slope_array.shape
    cube = []
    for k in range(1, nreads+1):
        rnd = random_image_slopes(nbits=noise_bits, rows=rows, columns=cols)
        cur = slope_array * k + rnd
        cube.append(cur)
    return np.uint16(cube)


def random_noisy_4d(slope_array=None, integrations=2, nreads=5, noise_bits=4):
    """
    Given ramp slopes for an image, create noisy 4-D ramp data.

    Parameter
    ---------
    slope_array: ndarray
        Slopes for each pixel

    integrations: int
        Number of integrations

    nreads: int
        Number of reads

    noise_bits: int
        Desired level of noise (in bits, i.e., 4 bits gives random noise 
        from 0 to 15.

    Return
    ------
    data_4d: ndarray
        A 4-D array of multiple integrations of noisy ramps.
    """
    if slope_array is None:
        return None

    data_4d = []
    for k in range(integrations):
        ramp = random_noisy_ramp(slope_array=mat, nreads=nreads, noise_bits=noise_bits)
        data_4d.append(ramp)
    return np.uint16(data_4d)


def simple_3d(mat, groups):
    """
    From a simple two dimensional array, create a simple linear cube ramp, where
    each group in the third dimension are simply a linear multiple of 'mat'.

    Parameter
    ---------
    mat: list
        A two dimensional array from which to build a cube.

    ngroups: int
        The number of groups to create from 'mat'.

    Return
    ------
    cube: ndarry
        A linear three dimensional cube.
    """
    cube_3d = []
    np_mat = np.array(mat, dtype=np.uint16)

    for g in range(groups):
        next_group = np_mat * (g+1)
        cube_3d.append(next_group)

    cube = np.array(cube_3d, dtype=np.uint16)
    return cube

def simple_4d():
    """
    Using a simple 2x2 matrix with integer entries, convert to a 4-D array
    mimicking 4-D RampModel data array.  The dimensions are
    (integrations, groups, rows, columns).

    From the simple 2x2 matrix, the first integration is a set of groups such
    that each successive group is simply the next linear multiple of 'mat'.

    From the simple 2x2 matrix, the first integration is a set of groups such
    that each successive group is simply the next linear multiple of 'mat', with
    a constant added to change the y-intecept.

    For each pixel, each integration will have the same slope, but with different
    y-intercept.

    Return
    ------
    cube_4d: ndarray
        A 4-D array with each integration a simple linear ramp of the first group
        in each integration.
    """
    n_int, ngroups, rows, cols = 2, 10, 2, 2
    basic_mat = [[5, 10], [15, 9]]

    # Create two simple ramps with the same slope, but different y intercept.
    int_1 = simple_3d(basic_mat, ngroups)
    int_2 = simple_3d(basic_mat, ngroups)
    int_2 = int_2 + 6

    cube_4d = np.array([int_1, int_2], dtype=np.uint16)

    return cube_4d

if __name__=="__main__":
    mat = simple_4d()
    print(mat)
