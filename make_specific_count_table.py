"""
for the simple case of varying the star mass and the closest pericenter distance
"""

import numpy as np
#from factory import PlotFactory as PF
import pickle

from matplotlib import pyplot as plot
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D

mass_range = [1.0/10, 1.0/5, 1.0/4, 1.0/3, 1.0/2, 2.0/3, 1.0, 3.0/2, 2.0, 5.0/2, 3.0]  #### Note 1.0 is before 0.67

mass_range = [round(mr,2) for mr in mass_range] #### <<<---- USE THIS

m_bins = np.zeros(len(mass_range) + 1)
m_bins[0] = 0.08  # only if log plot
m_bins[-1] = 3.5
for i, m in enumerate(mass_range[1:]):
    m_bins[i+1] = round(0.5 * (mass_range[i] + mass_range[i+1]), 2) # take average of consecutive values in mass range to set up bins

m_bins = np.array(m_bins)

print m_bins

"""
m1s = mass_range * 14 # replicate once for each pericenter distance
"""

peri_range = range(40, 301, 20) #### <<<---- USE THIS
p_bins = np.array(range(30, 311, 20)) ### for the color-mesh plot

print p_bins

"""
peris = np.zeros(len(peri_range) * len(mass_range))
for i, p in enumerate(peri_range):
    peris[i*11:(i+1)*11] = p

r_ins = np.zeros(len(peri_range) * len(mass_range))

r_ins[:33] = 10 # more inner particles for disks with closer passage
r_ins[33:] = 40
"""

######################################

# set up to be a string of a certain length

fn = open("count_table.txt", 'a')

def set_dash(num):
   s = ""
   for i in xrange(num):
       s += "-"
   return s

header = '{:^7}'.format("") # centered
dash = set_dash(7)

dashes = dash

for p in peri_range:
    next = '{:^7}'.format(str(p))
    header += (" | " + next)
    dashes += ("---" + dash)

fn.write(header + "\n")
fn.write(dashes + "\n")

#######################################
## ^^^ Add this to the print out ^^^ ##


######### READ THE RESULTS!!! ########

results = np.zeros((len(peri_range), len(mass_range)))

index = 1100
for i, p in enumerate(peri_range):
    #next = '{:^6}'.format(str(m)) + "|"
    for j, m in enumerate(mass_range):
        index += 1
        if index == 1134:
           index += 1000
        print index
        print i, j

        snap_dir = "sim00" + str(index)

        pickle_file = open(snap_dir + "/results.p", "rb")
        info = pickle.load(pickle_file) # Information Dictionary
        pickle_file.close()

        if p >= 1000:
           #results[i,j] = round(0.67 * info['std_dev'] * 2.0 / 5.0, 1)
           pass
        else:
           results[i,j] = round(info['mean'], 1)
           #results[i,j] = round(info['std_dev'] * 2.0 / 5.0, 1)

np.save("count_table.npy", results)

tr_results = results.transpose() ### numpy array here ###

np.save("count_tr_table.npy", tr_results)

########## WRITE A TABLE!!! ###########

index = 100
for i, m in enumerate(mass_range):
    next = '{:^6}'.format(str(m)) + "|"
    for j, p in enumerate(peri_range):
        index += 1
        print index

        count = '{:^7}'.format(str(tr_results[i,j]))

        if j == 3:
           next += ("   " + count)
           #next += (" | " + count)
        else:
           next += ("   " + count)
    fn.write(next + "\n")
    
########### MAKE A CONTOUR PLOT!!! ###########

levels = range(0, 41, 5)

fig, ax = plot.subplots()

#cmap = plot.get_cmap('YlGnBu_r')
cmap = plot.get_cmap('RdYlBu_r')

#contour_cmap = plot.get_cmap('YlGnBu')
#contour_cmap = plot.get_cmap('winter_r')
contour_cmap = plot.get_cmap('hot')

p_map = ax.pcolormesh(p_bins, m_bins, tr_results, cmap = cmap)
p_con = ax.contour(peri_range, mass_range, tr_results, levels = levels, linewidths = 2, cmap = contour_cmap)

ax.clabel(p_con, fmt = "%1.1f")

plot.gca().invert_yaxis()

ax.set_yscale('log')

#plot.title("Transferred Particle Counts (by %)")
plot.xlabel("Impact Parameter b (in AU)", fontsize = 15)
plot.ylabel("Mass of Disk Star m (in solar masses)", fontsize = 15)

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()

ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
ax.get_yaxis().get_major_formatter().labelOnlyBase = False

ax.set_xticks(peri_range)
ax.set_yticks(mass_range)

cbar = plot.colorbar(p_map)
cbar.set_label("Transferred Particles (in %)", fontsize = 15)

plot.tight_layout()

plot.show()
#plot.savefig("count_table_cmap.png")  # this doesn't work because of use of axes 'ax'

plot.clf()

########## MAKE A 3-D PLOT!!!!!! ##########

fig = plot.figure()
ax = Axes3D(fig)

ax.set_xticks(peri_range)
ax.set_yticks(mass_range)

X, Y = np.meshgrid(peri_range, mass_range)
Z = tr_results

surface = ax.plot_surface(X, Y, Z, rstride=1, cstride = 1, cmap = cmap, linewidth = 0, antialiased = False)

#ax.set_yscale('log')

#plot.title("Transferred Particle Counts (by %)")
plot.xlabel("Impact Parameter b (in AU)", fontsize = 15)
plot.ylabel("Mass of Disk Star m (in solar masses)", fontsize = 15)

ax.xaxis.set_label_position('top')
ax.xaxis.tick_top()

ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
ax.get_yaxis().get_major_formatter().labelOnlyBase = False

ax.set_xticks(peri_range)
ax.set_yticks(mass_range)

cbar = fig.colorbar(surface)
cbar.set_label("Transferred Particles (in %)", fontsize = 15)

plot.show()




