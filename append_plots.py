"""
add forgotten plots
"""

import numpy as np
from factory import PlotFactory as PF

for i in range(1134,1146):
    print i
    f = PF(i)
    f.mm_orb()


