#!/usr/bin/env python

import numpy as np


class RampInput4d:
    """
    The input ramp fit data for a 4-D image model with dimensions:
        (integrations, groups, rows, columns).
    This class contains the image data, the DQ data, the error data, and all
    processing parameters.
    """

    def __init__(self, label=None):
        """
        Initializes ramp fitting for a 4-d cub input data with default values.  The
        dimensions of the cube are (integrations, groups, rows, columns)
        """
        # Expand class variables as needed.
        self.label = label
        self.algorithm = "OLS"
        self.save_options = False

        self.multiprocessing = "None"

        self.data_4d = None
        self.groupdq = None
        self.err = None

    def __repr__(self):
        """
        Representation of the 4-D image model.
        """
        ostring = "    *** RampInput4d ***"
        if self.label is not None:
            ostring = f"{ostring}\nLabel = {self.label}"

        ostring = f"{ostring}\nRamp Fitting Algorithm = {self.algorithm}"
        if self.save_options is not None:
            ostring = f"{ostring}\nSave Options = {self.save_options}"

        ostring = f"{ostring}\nMultiprocessing = {self.multiprocessing}"
        if self.data_4d is not None:
            ostring = f"{ostring}\nShape of Data = {self.data_4d.shape}"

        if self.groupdq is not None:
            ostring = f"{ostring}\nShape of Group DQ = {self.groupdq.shape}"

        if self.err is not None:
            ostring = f"{ostring}\nShape of Errors = {self.err.shape}"

        return f"{DELIM}\n{ostring}\n{DELIM}"

    def set_algorithm(self, val):
        """
        Sets the ramp fitting algorithm.
        """
        upper_val = val.upper()
        possible_algorithms = ["OLS"]
        algos_str = ", ".join(possible_algorithms)
        if upper_val() not in possible_algorithms:
            raise ValueError(f"Ramp algorithm must be a string of type: {algos_str}")
        self.algorithm = upper_val

    def set_save_option(self, val):
        """
        Sets the ramp fitting save option.
        """
        if type(val) is not bool:
            raise TypeError("The save option must be boolean type.")
        self.save_option = val

    def set_4d_data(self, val):
        """
        Sets the ramp fitting 4-D image data.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The 4-D image array must be an ndarray.")
        self.data_4d = val

    def set_groupdq(self, val):
        """
        Sets the ramp fitting group DQ.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The group DQ must be an ndarray.")
        self.groupdq = val

    def set_err(self, val):
        """
        Sets the ramp fitting error.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The error must be an ndarray.")
        self.err = val


class RampOutput4d:
    """
    Place holder for now.
    """

    def __init__(self):
        """
        Place holder for now.
        """
        self.name = "RampOutput4d"


class RampInput3d:
    def __init__(self, label=None, ramp_input_data=None, integration=0):
        """
        This is a 3-D image cube with dimensions (groups, rows, columns).
        """
        self.label = label
        if ramp_input_data is None:
            self.algorithm = None
            self.save_options = None
            self.multiprocessing = None

            self.data_3d = None
            self.groupdq = None
            self.err = None

        if not isinstance(ramp_input_data, RampInput4d):
            raise TypeError("ramp_input_data must be an instance of RampInput4d")

        self.algorithm = ramp_input_data.algorithm
        self.save_options = ramp_input_data.save_options
        self.multiprocessing = ramp_input_data.multiprocessing

        self.data_3d = ramp_input_data.data_4d[integration, :, :, :]
        self.groupdq = ramp_input_data.groupdq[integration, :, :, :]
        self.err = ramp_input_data.err[integration, :, :, :]

    def __repr__(self):
        """
        Representation of the 4-D image model.
        """
        ostring = "    *** RampInput4d ***"
        if self.label is not None:
            ostring = f"{ostring}\nLabel = {self.label}"
        ostring = f"{ostring}\nRamp Fitting Algorithm = {self.algorithm}"
        if self.save_options is not None:
            ostring = f"{ostring}\nSave Options = {self.save_options}"
        ostring = f"{ostring}\nMultiprocessing = {self.multiprocessing}"
        if self.data_3d is not None:
            ostring = f"{ostring}\nShape of Data = {self.data_3d.shape}"
        if self.groupdq is not None:
            ostring = f"{ostring}\nShape of Group DQ = {self.groupdq.shape}"
        if self.err is not None:
            ostring = f"{ostring}\nShape of Errors = {self.err.shape}"

        return f"{DELIM}\n{ostring}\n{DELIM}"

    def set_algorithm(self, val):
        """
        Sets the ramp fitting algorithm.
        """
        upper_val = val.upper()
        possible_algorithms = ["OLS"]
        algos_str = ", ".join(possible_algorithms)
        if upper_val() not in possible_algorithms:
            raise ValueError(f"Ramp algorithm must be a string of type: {algos_str}")
        self.algorithm = upper_val

    def set_save_option(self, val):
        """
        Sets the ramp fitting save option.
        """
        if type(val) is not bool:
            raise TypeError("The save option must be boolean type.")
        self.save_option = val

    def set_3d_data(self, val):
        """
        Sets the ramp fitting 3-D image data.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The 3-D image array must be an ndarray.")
        self.data_3d = val

    def set_groupdq(self, val):
        """
        Sets the ramp fitting group DQ.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The group DQ must be an ndarray.")
        self.groupdq = val

    def set_err(self, val):
        """
        Sets the ramp fitting error.
        """
        if type(val) is not np.ndarray:
            raise TypeError("The error must be an ndarray.")
        self.err = val


class RampOutput3d:
    """
    Place holder for now.
    """

    def __init__(self):
        self.name = "RampOutput3d"
