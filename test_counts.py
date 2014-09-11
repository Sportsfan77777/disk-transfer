"""
test program (count calculations)
"""


from factory import PlotFactory as PF
import numpy as np

x = 230
n = 10

"""
for i in range(x, x+n):
    f = PF(i, snapshot_base = "snap", count = 3)
    f.simple_mask()
    f.sort()
    f.count()
    f.mm_pos()
"""

counts = []
for i in range(x, x+n):
    f = PF(i, snapshot_base = "snap", count = 3)
    final_count = len(f.passing[-3])
    counts.append(final_count)

frac_counts = [x/5.0 for x in counts]

print "Counts:", frac_counts
print "Sorted:", np.sort(frac_counts)
print "Median:", np.median(frac_counts)
print "Mean:", np.mean(frac_counts)
print "Standard Deviation:", np.std(frac_counts)
print

middle = np.sort(frac_counts)[1:-1]

print "Middle:", middle
print "Median:", np.median(middle)
print "Mean:", np.mean(middle)
print "Standard Deviation:", np.std(middle)

