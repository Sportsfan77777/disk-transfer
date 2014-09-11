"""
for the simple case of varying the star mass and the closest pericenter distance
"""

import numpy as np
#from factory import PlotFactory as PF
import pickle

mass_range = [1.0/10, 1.0/5, 1.0/4, 1.0/3, 1.0/2, 2.0/3, 1.0, 3.0/2, 2.0, 5.0/2, 3.0]  #### Note 1.0 is before 0.67

mass_range = [round(mr,3) for mr in mass_range] #### <<<---- USE THIS

"""
m1s = mass_range * 14 # replicate once for each pericenter distance
"""

peri_range = range(40, 141, 20) #### <<<---- USE THIS

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

results = np.zeros((len(peri_range), len(mass_range)))

index = 3100
for i, p in enumerate(peri_range):
    #next = '{:^6}'.format(str(m)) + "|"
    for j, m in enumerate(mass_range):
        index += 1
        print index
        print i, j

        #f = PF(index)
        #final = f.get_final_count(last = 3)
        #adj = round( 100.0 * final / f.info.n_disk, 1) # fractional count

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

tr_results = results.transpose()

np.save("count_tr_table.npy", tr_results)

index = 3100
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
    


