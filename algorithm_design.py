"""
Receive 4-D RampModel with processing information:
"""

class IntermediateComputations(model):
    def __init__(self):
        shape.orig_dim = model.data.shape

def ramp_fit(model, proc_options):
    orig_dimensions = model.data.shape
    miri_answer = miri_correction(model)
