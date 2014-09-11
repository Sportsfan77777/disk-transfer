"""
test program (executions)
"""

import sys ##### <----fix this (to take in a string argument as an identifier) 
import subprocess
import random

import numpy as np

# option parser organizers
def append(arr, string):
    """ simple append structured to option parser format """
    items = string.split() # split by ' '    <<<<<------- Use shlex.split(string) instead??????
    for item in items:
        arr.append(item)

def add_option(o_str, o):
    """ add an 'optional' option """
    # careful with which 'add'
    if o is not None:
       option = "-%s %s" % (o_str, str(o))
       if len(o_str) > 1:
          add("-" + option) # "--o
       else:
          add(option) # "-option

def to_shell(command):
    string = ""
    for o in command:
        string += (o + " ")
    string = string[:-1] # remove trailing " "
    return string


options = ["sim_number", "num_format", "num_T", "nb_end", "m0", "m1", "rel_force", "br_dt", "eta", 
           "r_in", "r_out", "n_disk", "power", "peri", "ecc", "incl", "omega", "n", "seed"]

# Initialize options to None (so that they will not be called unless needed)
init_none = ""
for o in options:
    init_none += (o + " = ")
init_command = init_none + " None"
exec(init_command) # one time only <<<--- BE CAREFUL WITH EXEC!!!!!!!!!

if len(sys.argv) > 1:
   name = sys.argv[1] # 1st argument: a string name for the file containing the list of commands
else:
   # select random number (not useful in combination with the bash script)
   random_num = random.randint(0, 99999)
   name = "%05d" % (random_num)

##############################################################################################
##############################################################################################
##########################             MODIFY THIS              ##############################

x = 10101 # snap number (this format will soon be deprecated)
iterations = 154

# Parameters constant through all simulations 
# (or 'initial conditions' for zeroth simulation)
# (or array of parameters that vary)

n = 162

#num_format = 3
#snap_dir = "snap"

# to make quick
#nb_end = 0.95
#num_T = 1.5

n_disk = 500

#seeds = [11,21345,153,6136,3661836,13,36,47,983,73151,
#         22,98721,551,1234,4321362,80,25,23,542,90901,
#         55,10101,355,9919,1999135,90,65,34,562,81113,
#         45,41012,789,2468,5818528,89,72,76,781,34521]

#peris = [100, 120, 140, 150, 160, 180, 200, 225, 250, 275, 300, 325]
#peri = 200

mass_range = [1.0/10, 1.0/5, 1.0/4, 1.0/3, 1.0/2, 2.0/3, 1.0, 3.0/2, 2.0, 5.0/2, 3.0]

mass_range = [round(mr,3) for mr in mass_range]

m0s = mass_range * 14 # replicate once for each pericenter distance

peri_range = range(40, 301, 20)

peris = np.zeros(len(peri_range) * len(mass_range))
for i, p in enumerate(peri_range):
    peris[i*11:(i+1)*11] = p

r_ins = np.zeros(len(peri_range) * len(mass_range))

# Not Parabolic!!! ######
ecc = 1.5

"""
for i in range(len(r_ins)):
    m = m1s[i]
    p = peris[i]
    if p >= 180:
       r_in[i] = 40
    elif m >= 1.0:
       r_in[i] = 40
    elif m >= 0.6:
       if p >= 120:
          r_in[i] = 40          
       else:
          r_in[i] = 10
    elif m >= 0.3:
       if p >= :
          r_in[i] = 40          
       else:
          r_in[i] = 10
    elif m >= 0.2:
       if p >= :
          r_in[i] = 40
       else:
          r_in[i] = 10
       
"""

"""
r_ins[:33] = 10 # more inner particles for disks with closer passage
r_ins[33:] = 40
"""

r_ins[:] = 10


##############################################################################################
##############################################################################################
##############################################################################################

# Set up each option not dependent on i outside of the loop

#command = """amuse.sh disk_flyby_multiple_code_lestrade11.py""".split()
command = """amuse.sh disk_flyby_multiprocessing_code_lestrade11.py""".split()

# Deprecate for now, since re-written variables will not be overwritten
"""
add = lambda s1 : append(command, s1) # concern about 'late binding' problems?
for o in options:
     add_option(o, eval(o)) ####### CAREFUL WITH EVAL #######
"""

processes = []
outfiles = []
for i in range(x, x + iterations):
     # Reset command
     command_i = command[:] # copy array

     ##############################################################################################
     ##############################################################################################
     ##########################             MODIFY THIS              ##############################

     # Calculate parameters (dependent on i) if necessary [PRE]

     r_in = r_ins[i-x]
     peri = peris[i-x]
     m0 = m0s[i-x]
     
     ##############################################################################################
     ##############################################################################################
     ##############################################################################################

     # Partial Function
     add = lambda s2 : append(command_i, s2)

     # Add each option dependent on i separately
     #add("""--fout outfile%03d.hdf5""" % (i)) # string options are read differently by python and bash?
     #add("""--snap_dir snap""")

     sim_number = i

     for o in options:
         add_option(o, eval(o)) ####### CAREFUL WITH EVAL #######

     shell_command = to_shell(command_i)
     #print "Shell Command:"
     print shell_command

     # Track all processes
     processes.append(shell_command)
     
     # Track all outfiles
     outfile = """sim%06d.out""" % (i)
     outfiles.append(outfile)

     ##############################################################################################
     ##############################################################################################
     ##########################             MODIFY THIS              ##############################

     # Calculate parameters (dependent on i) if necessary [POST]

     ##############################################################################################
     ##############################################################################################
     ##############################################################################################

# Write 'commands' and 'outfiles' to separate files in "commands" directory

fn = "commands/commands_%s.txt" % (name)
f = open(fn, "w")
for p in processes:
    f.write(p + "\n")
f.close()

fn = "commands/commands_out_%s.txt" % (name)
f = open(fn, "w")
for out in outfiles:
    f.write(out + "\n")
f.close()



 


