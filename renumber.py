"""
're-number' a simulation
essentially move a simulation directory to a different number
as well as all of the contents of that directory which are labeled with that number
"""

"""
Things to change
[6](A) 'Grand' Directory <<--- Do this last, lol
[5](B) 'Grand' outfile.txt, outfile.hdf5
[4](C) 'Grand' info.txt, info.p

[3](A) sub-simulation directories (all of them) <-- glob.glob
[2](B) sub-simulation outfile.txt, outfile.hdf5
[1](C) sub-simulation info.txt, info.p
"""

from amuse.units import units,quantities
import glob
import sys
import shutil
import pickle

args = sys.argv
if not (len(args) == 3):
   raise Exception("Use two arguments!!!!")
else:
   print args
   
old_n = "%06d" % (int(args[1]))
new_n = "%06d" % (int(args[2]))

def switch(old_n, new_n):
   grand_path = "sim" + old_n
   sub_path_wild = grand_path + "/" + grand_path + "_*"

   # Find all sub-simulations + grand_path
   sub_simulations = sorted(glob.glob(sub_path_wild))
   sub_simulations.append(grand_path)

   print sub_simulations

   for sub_sim in sub_simulations:
       #[1, 4]
       #1a
       info_pickle_fn = sub_sim + "/info.p"
       info_pickle_tmp_fn = sub_sim + "/info_tmp.p"

       info_pickle = open(info_pickle_fn, "rb")
       info_pickle_tmp = open(info_pickle_tmp_fn, "wb")

       info_dict = pickle.load(info_pickle)
       info_dict.sim_number = new_n
       # also change snap_dir and outfile ###### <--- do this next

       pickle.dump(info_dict, info_pickle_tmp)
       info_pickle.close()

       shutil.move(info_pickle_tmp_fn, info_pickle_fn)

       #1b
       info_txt_fn = sub_sim + "/info.txt"
       info_txt_tmp_fn = sub_sim + "/info_tmp.txt"

       info_txt = open(sub_sim + "/info.txt", "r")
       info_txt_tmp = open(sub_sim + "/info_tmp.txt", "w")

       txt = info_txt.read()
       info_txt_tmp.write(txt.replace(old_n, new_n))
    
       info_txt.close()
       info_txt_tmp.close()
       shutil.move(info_txt_tmp_fn, info_txt_fn) ################## <---------- save fns, DON'T USE files!!!!!!!
    
       #[2, 5]
       outfiles_wild = sub_sim + "/outfile*"
       outfiles = glob.glob(outfiles_wild)
       for out_f in outfiles:
           b = out_f.rfind("/")
           end = out_f[b:]
           new_out_f = out_f[:b] + end.replace(old_n, new_n)
           shutil.move(out_f, new_out_f)
       #[3, 6]
       if i == len(sub_simulations) - 1:
          new_path = sub_sim.replace(old_n, new_n)
          shutil.move(sub_sim, new_path)
       else:
          new_path = sub_sim.replace(old_n, new_n)
          shutil.move(sub_sim, new_path)

   #### This could be combined with the sub_simulations #### (duh...) (just append grand_path to sub_simulations)
   """
   #[4]
   #4a
   info_pickle = open(grand_path + "/info.p", "rwb")

   info_dict = pickle.load(info_pickle)
   info_dict.sim_number = new_n

   pickle.dump(info_dict, info_pickle)
   info_pickle.close()

   #4b
   info_txt = open(grand_path + "/info.txt", "r")
   info_tmp = open(grand_path + "/info_tmp.txt", "w")

   txt = info_txt.read()
   info_tmp.write(txt.replace(old_n, new_n))
    
   info_txt.close()
   info_tmp.close()
   shutil.move(info_tmp, info_txt)

   #[5]
   outfiles_wild = grand_path + "/outfile*"
   outfiles = glob.glob(outfiles_wild)
   for out_f in outfiles:
       new_out_f = out_f.replace(old_n, new_n)
       shutil.move(out_f, new_out_f)

   #[6]
   new_path = grand_path.replace(old_n, new_n)
   shutil.move(grand_path, new_path)
   """

switch(old_n, new_n) #### call ####



