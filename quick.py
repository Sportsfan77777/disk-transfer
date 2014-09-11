"""
test program (results)
"""

from factory import PlotFactory as PF
import plot_snaps_all_05 as ps


f = PF(228)

a = f.unbound[240] # ~ 6000 yrs
b = f.unbound[-1] # last time step

changes = [c for c,d in zip(a,b) if d not in a] # particles that are unbound at the end, but not when the n-body ends






