"""
test program (calculations)
"""


from factory import PlotFactory as PF

x = 981
n = 5

"""
for i in range(x, x+n):
    f = PF(i, snapshot_base = "snap", count = 3)
    f.simple_mask()
    f.sort()
    f.count()
    f.mm_pos()
"""

for i in range(x, x+n):
    f = PF(i, snapshot_base = "snap", count = 3)
    f.hist()

