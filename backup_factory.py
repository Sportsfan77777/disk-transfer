"""
Sort of like a factory, this program is intended to call the functions
in a particular plotting file, which in this program, is referred to as 'ps'

Note: ps = plot_snaps methods

Shortcut Methods:

save()
restore()
possess()

mm_pos()
mm_orb()
mm_sma()
counts()
hist()
"""

# Not sure if all of these are necessary

from amuse.lab import *
from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.datamodel import Particles
from amuse.io import read_set_from_file

import subprocess
from collections import Iterable as iterable
from shutil import rmtree
import glob
from os import mkdir, symlink
import pickle

import plot_snaps_all_05 as ps

import numpy as np
from math import e

# Figure out how to make this into a class!!!!!!

class PlotFactory:  

    # Globals, yay! (class variables)
    CENTRAL_STAR = 0
    PASSING_STAR = 1

    UNBOUND_INDEX = 0
    CENTRAL_INDEX = 1
    PASSING_INDEX = 2
    BOTH_INDEX = 3
    
    # Colors
    c_STARS = 1
    c_CENTRAL = 2 # Originally had disk
    c_PASSING = 3 # Passing by star with disk
    c_BOTH = 4
    c_UNBOUND = 5

    sorted_set_fn = "/orbit_sorted_set_keys.npy"
    mask_fn = "/mask.npy"

    def __init__(self, number, snapshot_dir = None, output_loc = None, snapshot_loc = None):

        # Output File
        number = '%03d' % number

        self.f1 = "outfile" + number + ".hdf5"

        if output_loc is not None:
           self.f1 = output_loc + "/" + self.f1

        # Snapshot Directory

        if snapshot_dir is None:
           snapshot_dir = "snap" + number

        if snapshot_loc is not None:
           snapshot_dir = snapshot_loc + "/" + snapshot_dir

        self.snapshot_dir = snapshot_dir # Updates Here
        ps.set_snapshot_dir(snapshot_dir) # Updates Elsewhere

        # Bodies

        self.bodies = ps.parse_file(self.f1)

        # Simulation Info
        self.read_info_file()

        # Time

        self.times = quantities.AdaptingVectorQuantity()

        # Initialize States

        self.init_states()
        self.num_snapshots = len(self.times)
        self.save_star_separation()

        # Sorted Bodies + Keys (Change this to particle sets instead)

        self.unbound = np.zeros(len(self.times), dtype = np.ndarray) # These are arrays of Particle sets
        self.central = np.zeros(len(self.times), dtype = np.ndarray)
        self.passing = np.zeros(len(self.times), dtype = np.ndarray) 
        self.both = np.zeros(len(self.times), dtype = np.ndarray) 

        self.tmp_mask = np.zeros(len(self.times)) # for before the mask is applied
        self.mask = np.zeros(len(self.times)) # 1 indicates 'sorting' has occured at particular timestep
        # This mask needs to be incorporated into functions using sorted sets
        # For instance, the self.times array will need to be restricted to masked values
        ################# DO NOT MODIFY THE MASK DIRECTLY!!!!!! (switch to private?) ################

        # Restore Sorts

        print " *** Restoring *** "
        self.restore()

        # Misc
        self.tmp_movie_dir = "tmp_movies_for_factory"


#####################################################################################################
# 1

# This needs to be re-written

    def plot_all(self):
        """ Plot Everything (step by step) """
        self.init_orbital_parameter_sorted_sets()
        
        self.make_movies_positon()
        self.make_movies_orbital_elements()
        self.make_movies_semi_major_axes()
        self.make_count_curve()
        self.make_histograms()
        # Go on...

    def plot_all_from_files():
        """ Plot Everything! """
        ps.plot_all_snaps(None, bodies = bodies)

#####################################################################################################
# 2a

    def init_states(self):
        """ Initializes the time array, the bodies at timestep 1 (start+1), and the bodies at timestep -1 (end) """
        # Reset if called previously (or would prompting user be better?)
        self.zero_state = None
        self.initial_state = None
        self.final_state = None

        for b in self.bodies.history:
            # Mark Time
            self.times.append(b.collection_attributes.timestamp)

            # Set States
            if self.zero_state is None:
                self.zero_state = b
            elif self.initial_state is None:
                self.initial_state = b

            self.final_state = b # re-update until history ends

        # Identify limits on inner and outer disk radii (######## Eventually, replace this with reading in the info ########)
        self.min_a = 100000.0|units.AU
        self.max_a = 0.0|units.AU
        two_stars = self.initial_state[0:2]
        planetesimals = self.initial_state[2:]
        shifted = self.shift_to_star(planetesimals, two_stars, star_choice = 0)
        for b in shifted:
            a, e, T = ps.orbital_parameters(b.position, b.velocity, self.info.m0)
            if (a > self.max_a):
               self.max_a = a
            if (a < self.min_a):
               self.min_a = a

# 2b

    def read_info_file(self):
        """ loads parameters for simulation from info.txt file """
        pickle_file = open(self.snapshot_dir + "/info.p", "rb")
        self.info = pickle.load(pickle_file) # Information Dictionary
        pickle_file.close()

# 2c

    def set_mask(self, ranges):
        """
        set mask for creating and using the orbital paramter sorted sets
        e.g.
        ranges = [range(0,40,8), range(40,60,1), [99, 150], range(60,432,20)]
        bad values will be ignored

        Technically, ranges can overlap. All of the masks will be 'or'-ed together along with the existing mask
        i.e. 100100100 + 101010101 + 000101000 = 101111101

        Future Modifications: Allow for inputs to be either 'indices' or 'times'
        """

        for r in ranges:
            # First, make range nice
            r = [round(x) for x in r if x >= 0 and x < len(self.times)] # restrict values to ints in range [0, len)
            # Apply Mask to Temporary Mask. The real mask will be updated when applied to the sorting.
            self.tmp_mask[r] = 1

# 2d

    def save_star_separation(self):
        """ calculate and plot separation between two stars """
        self.star_separation = ps.plot_distances(self.times, self.bodies) 
        # the input is messed up ([0:2] does not work)

# 2e

    def q(self):
        """ shortcut for print_quartiles """
        return self.print_quartiles()

    def print_quartiles(self):
        """ return critical points in orbit """
        start = self.star_separation[0] # not yet --> .value_in(units.AU)
        start_index = 0
        start_time = self.times[0].value_in(1000 * units.yr)
        final = self.star_separation[-1].value_in(units.AU)
        final_index = len(self.times) - 1
        final_time = self.times[-1].value_in(1000 * units.yr)

        closest = min(self.star_separation) # not yet --> .value_in(units.AU)
        closest_index = ([i for i,x in enumerate(self.star_separation) if x == closest])[0]
        closest_time = self.times[closest_index].value_in(1000 * units.yr)

        critical = start - (start - closest)*(1/e) # "Scale Height-ish"
        start = start.value_in(units.AU) # now, do it
        closest = closest.value_in(units.AU) # now, do it
        
        critical_one_min = min(abs(self.star_separation[:closest_index] - critical))
        critical_two_min = min(abs(self.star_separation[closest_index:] - critical))
        critical_one_index = ([i for i,x in enumerate(abs(self.star_separation[:closest_index] - critical)) \
                                                   if x == critical_one_min])[0]
        critical_two_index = closest_index \
                           + ([i for i,x in enumerate(abs(self.star_separation[closest_index:] - critical)) \
                                                   if x == critical_two_min])[0]
        critical_one = self.star_separation[critical_one_index].value_in(units.AU)
        critical_two = self.star_separation[critical_two_index].value_in(units.AU)
        critical_one_time = self.times[critical_one_index].value_in(1000 * units.yr)
        critical_two_time = self.times[critical_two_index].value_in(1000 * units.yr)

        print "Start: ", start_time, " | Index: ", start_index,  " | Distance", start
        print "Encounter Start:", critical_one_time, " | Index: ", critical_one_index, " | Distance", critical_one
        print "Closest:", closest_time, " | Index: ", closest_index, " | Distance", closest
        print "Encounter End:", critical_two_time, " | Index: ", critical_two_index, " | Distance", critical_two
        print "End: ", final_time, " | Index: ", final_index, " | Distance", final

        return [critical_one_index, closest_index, critical_two_index, final_index]

    def simple_mask(self):
        """ sets up sparse mask that is more dense at closest encounter and end """
        quartiles = self.q()
        self.set_mask([range(quartiles[0], quartiles[2], 8), range(quartiles[1]-5, quartiles[1]+5),
                       range(quartiles[2], quartiles[3], 35), range(quartiles[3]-5, quartiles[3])])


#####################################################################################################
# 3a

    def sort(self):
        """ shortcut for sorting into orbital parameter sorted sets """
        self.init_orbital_parameter_sorted_sets()

    def init_orbital_parameter_sorted_sets(self):
        """ Initialize particle sets which sort the particles by their bounds at each timestep """
        # Update current mask
        self.mask = self.tmp_mask

        for i,b in enumerate(self.bodies.history):
            two_stars = b[0:2]
            planetesimals = b[2:]

            if (self.mask[i] and (not isinstance(self.unbound[i], iterable)) ):
                 print "Iteration ", i
                 # iterable is defined in the header
                 # only sort if mask says to sort and array does not exist yet
                 unbound, central, passing, both = ps.orbital_parameter_sorted_sets(planetesimals, two_stars)
                 (self.unbound)[i] = unbound
                 (self.central)[i] = central
                 (self.passing)[i] = passing
                 (self.both)[i] = both

        print " **** Saving Sorted Particle Sets **** "
        self.save()

        # Also, partition 'both'
        self.partition_both()
        self.save()

# 3b

    def save(self):
       """ shortcut for saving orbital parameter sorted sets """
       self.save_orbital_parameter_sorted_sets()

    def save_orbital_parameter_sorted_sets(self):
       """ Saves keys from sorted sets using Numpy 2-D array (4 x num_snapshots) """
       num_snapshots = len(self.times)
       sorted_keys = np.zeros( (4, num_snapshots), dtype = np.ndarray)

       for i in xrange(num_snapshots):
           # Only save sets that have been calculated

           keys_unbound = []
           keys_central = []
           keys_passing = []
           keys_both = []

           if (isinstance(self.unbound[i], iterable)):
              keys_unbound = np.array(self.unbound[i].key)
           if (isinstance(self.central[i], iterable)):
              keys_central = np.array(self.central[i].key)
           if (isinstance(self.passing[i], iterable)):
              keys_passing = np.array(self.passing[i].key)
           if (isinstance(self.both[i], iterable)):
              keys_both = np.array(self.both[i].key)

           sorted_keys[PlotFactory.UNBOUND_INDEX][i] = keys_unbound
           sorted_keys[PlotFactory.CENTRAL_INDEX][i] = keys_central
           sorted_keys[PlotFactory.PASSING_INDEX][i] = keys_passing
           sorted_keys[PlotFactory.BOTH_INDEX][i] = keys_both

           #sorted_keys[PlotFactory.UNBOUND_INDEX][i][1] = _unbound
           #sorted_keys[PlotFactory.CENTRAL_INDEX][i][1] = sma_central
           #sorted_keys[PlotFactory.PASSING_INDEX][i][1] = sma_passing
           #sorted_keys[PlotFactory.BOTH_INDEX][i][PlotFactory.SMA_INDEX] = _both

       fn = self.snapshot_dir + self.sorted_set_fn
       np.save(fn, sorted_keys)

       # Save the associated mask also

       fn = self.snapshot_dir + self.mask_fn
       np.save(fn, self.mask)

       

# 3c

    def restore(self):
       """ shortcut for restoring orbital parameter sorted sets """
       self.restore_orbital_parameter_sorted_sets()

    def restore_orbital_parameter_sorted_sets(self):
       """ Restores sorted_sets using Numpy 2-D array of saved keys (4 x num_snapshots) """

       mask = self.restore_mask()

       sorted_keys = self.possess_keys_for_orbital_parameter_sorted_sets()

       if sorted_keys is None:
          return # No Restore ---> End

       for i, b in enumerate(self.bodies.history):
          if ( mask[i] ):
              # Retrieve only sets that have been calculated (ignore empty non-iterables)
              (self.unbound)[i] = b.select(lambda x : x in sorted_keys[PlotFactory.UNBOUND_INDEX][i], ["key"])
              (self.central)[i] = b.select(lambda x : x in sorted_keys[PlotFactory.CENTRAL_INDEX][i], ["key"])
              (self.passing)[i] = b.select(lambda x : x in sorted_keys[PlotFactory.PASSING_INDEX][i], ["key"])
              (self.both)[i] = b.select(lambda x : x in sorted_keys[PlotFactory.BOTH_INDEX][i], ["key"])
              # Add attribute 'sm-axis'
              # Update Mask
              print i
              self.tmp_mask[i] = 1
       self.mask = self.tmp_mask

    def restore_mask(self):
       """ 
       Restores keys to (4 x num_snapshots) Numpy 2-D array so that
       only the sorts that are needed are retrieved
       """
       try:
          fn = self.snapshot_dir + self.mask_fn
          mask = np.load(fn)
          return mask
       except IOError:
          print " *** no mask to restore *** "
          return None

# 3d

    def possess(self):
       """ shortcut for possessing keys for orbital parameter sorted sets """
       return self.possess_keys_for_orbital_parameter_sorted_sets()

    def possess_keys_for_orbital_parameter_sorted_sets(self):
       """ 
       Restores keys to (4 x num_snapshots) Numpy 2-D array so that
       only the sorts that are needed are retrieved
       """
       try:
          fn = self.snapshot_dir + self.sorted_set_fn
          sorted_keys = np.load(fn)
          return sorted_keys
       except IOError:
          print " *** nothing to restore *** "
          return None

# 3e
    def compress(self, array, mask):
       """ shortcut for compress array """
       return self.compress_array(array, mask)

    def compress_array(self, array, mask):
       """ Given an array and a mask of equal length, return the array where the masked entries are selected"""
       copy = list(array)
       compressed = [ x for i,x in enumerate(copy) if mask[i] ]
       # e.g. array = [1,2,3,4,5], mask = [1,0,0,1,0] ---> compressed = [1,4]
       vector = quantities.AdaptingVectorQuantity() # formatting issues (is there a better way?)
       for x in compressed:
           vector.append(x)
       return vector

# 3f

    def partition_both(self, start_time = None, start_index = None):
       """ splits an existing array of particles sorted into both after a given timestep, clearing the 'both' sort """
       split = start_index
       if start_index is not None:
          split = start_index
       elif start_time is not None:
          times = self.compress(self.times) #### Is this needed? ####
          earliest = min(abs(times - start_time))
          print "Earliest to Partition:", earliest
          split = ([i for i,x in enumerate(times - start_time) if abs(x) == earliest])[0]
       else:
          target_separation = 1.3 * self.star_separation[0]
          earliest = min(abs(self.star_separation - target_separation))
          split = ([i for i,x in enumerate(self.star_separation - target_separation) if abs(x) == earliest])[0]

       for i, (both, bodies) in enumerate(zip(self.both, self.bodies.history)):
          if (i > split) and (self.mask[i]):
              two_stars = bodies[0:2]
              print "Iteration ", i
              for b in both:    
                  dist_c = (b.position - two_stars[0].position).lengths()
                  dist_p = (b.position - two_stars[1].position).lengths()

                  if (1.2 * dist_p < dist_c):
                      self.passing[i].add_particle(b) ### modify this to include masses (treat in terms of specific forces = potentials)
                      self.both[i].remove_particle(b)
                  elif (1.2 * dist_c < dist_p):
                      self.central[i].add_particle(b)
                      self.both[i].remove_particle(b)
                  else:
                     # not clearly bound to a specific star (using simple calculation)
                     pass
       
       
       
#####################################################################################################
# 4 (Switch Coordinate Systems!)

    def shift_to_star(self, bodies, two_stars, star_choice = 1):
       """ shift bodies with stars as an awkward, but neat input """
       # Note: the parameters can be simplified
       
       # First, account for stupid inputs
       if (star_choice != 0 and star_choice != 1):
           star_choice = 1
          
       position = two_stars[star_choice].position
       velocity = two_stars[star_choice].velocity
       return self.shift(bodies, position, velocity)
           
       
    def shift(self, bodies, pos, vel):
       """ shift bodies to new coordinate system """
       #print pos
       #print
       #print bodies.position
       c_bodies = bodies.copy() # You must make copy to avoid writing to history
       c_bodies.position -= pos
       c_bodies.velocity -= vel
       return c_bodies    

#####################################################################################################
# 5 (Plotting!)

# 5

# Set up movie prompt showing options?
    """ general make movie terminal command """
    # Doesn't Work. Blame subprocess?
    def make_movies(self, path, output, fps = 6, consec = False):
         """ Given path, make movies """
         # First make sure %03d is consecutive. (Create tmp directory of consecutive files)
         if not consec:
            try:
               mkdir(self.tmp_movie_dir)
               print " *** creating temporary directory ***"
            except OSError:
               print "\t(" , self.tmp_movie_dir, "already exists)"

            a = path.rfind('_')
            b = path.rfind('.')
            asterisk_path = path.replace(path[a:b], '*')
            
            movie_files = sorted(glob.glob(asterisk_path))
            print movie_files
            path = self.tmp_movie_dir + "/tmp_%03d.png"
            # Re-number all of the .png files (using symbolic links [ ln -s in the terminal ] )
            for i,f in enumerate(movie_files):
                print f
                lns_f = self.tmp_movie_dir + "/tmp_%03d.png" % i
                print lns_f
                symlink("../" + f, lns_f) # <--- Symbolic links are not working (Update: Now, they are!)

         # Now, we can make the movie
         command = "ffmpeg -f image2 -r %d -i %s -vcodec mpeg4 -y %s" % (fps, path, output)
         movie_pipe = subprocess.call(command.split()) # 'call' takes a list of args, e.g. ['ffmpeg', '-f', etc.]

         # Delete temporary directory
         if not consec:
            print " *** deleting temporary directory ***"
            rmtree(self.tmp_movie_dir) ############# CAREFUL!!!!!! ############## (but it works!!!!)
            pass

# 5a

    def mm_pos(self, movie = 1, consec = False):
         self.make_movies_position(movie)

    def make_movies_position(self, movie = 1, consec = False):
         # Add: First, check if the sorted_sets are not empty
         # If they are empty prompt the user
    
         for i, (b, unbound, central, passing, both) in \
                  enumerate(zip(self.bodies.history, self.unbound, self.central, self.passing, self.both)):
             # Note, the function call below will shift positions automatically
             # It is not necessary to shift the velocities. That is why the self.shift function is not used.
             if (self.mask[i]):
                two_stars = b[0:2]

                ps.make_colorcoded_position_scatter_plot(i, self.times[i], two_stars, \
                                     unbound, central, passing, both)
                   
         if movie == 1:               
            # (1) POSITION 
            pos_path = self.snapshot_dir + "/spdd/spdd_%3d.png"
            pos_movie = self.snapshot_dir + "/snaps.mp4" # Note the locations are different
            self.make_movies(pos_path, pos_movie, consec = consec)

            # (2) ZOOOOMED POSTION
            zoom_path = self.snapshot_dir + "/zoom_spdd/zoom_spdd_%3d.png"
            zoom_pos_movie = self.snapshot_dir + "/snaps_zoom.mp4" # Note the locations are different
            self.make_movies(zoom_path, zoom_pos_movie, consec = consec)

# 5b
    def mm_orb(self, movie = 1, consec = False):
        self.make_movies_orbital_elements()

    def make_movies_orbital_elements(self, movie = 1, consec = False):
         # Add: First, check if the sorted_sets are not empty
         # If they are empty prompt the user

         for i, (b, passing) in enumerate(zip(self.bodies.history, self.passing)):
            if (self.mask[i]):
               # ^^^ this could be changed to a user-selected prompt from the mask (showing the times)
               # Store sma? (or re-store sma?)

               # Note, using the mask isn't as helpful as it could be
               # The sm-axes are still re-calculated. This is bad. They should be stored in a file

               two_stars = b[0:2]
               planetesimals = b[2:]
               
               mass = self.info.m1 # This can be taken out of the loop
               bodies = passing.copy()

               if len(bodies) > 0:
                    print "Iteration ", i
                    bodies.position -= two_stars[1].position
                    bodies.velocity -= two_stars[1].velocity
                    ps.make_orbital_elements_scatter_plot(bodies, mass, i = i)

         if movie == 1:
            # (3a) ORBITAL ELEMENTS (ecc vs a)
            orb_path = self.snapshot_dir + "/orb/orb_%3d.png"
            orb_movie = self.snapshot_dir + "/orbital_elements.mp4"
            self.make_movies(orb_path, orb_movie, consec = consec)

            # (3b) ORBITAL ELEMENTS (ecc vs a)
            orb_zoom_path = self.snapshot_dir + "/orb_zoom/orb_%3d.png"
            orb_zoom_movie = self.snapshot_dir + "/orbital_elements_zoom.mp4"
            self.make_movies(orb_zoom_path, orb_zoom_movie, consec = consec)

            # (3c) ORBITAL ELEMENTS (ecc vs a)
            orb_log_path = self.snapshot_dir + "/orb_log/orb_%3d.png"
            orb_log_movie = self.snapshot_dir + "/orbital_elements_log.mp4"
            self.make_movies(orb_log_path, orb_log_movie, consec = consec)

            # (3d) ORBITAL ELEMENTS (inc vs a)
            orb_inc_path = self.snapshot_dir + "/orb_inc/orb_%3d.png"
            orb_inc_movie = self.snapshot_dir + "/orbital_elements_inc.mp4"
            self.make_movies(orb_inc_path, orb_inc_movie, consec = consec)

            # (3e) ORBITAL ELEMENTS (q vs a)
            orb_q_path = self.snapshot_dir + "/orb_pericenter_distance_q/orb_%3d.png"
            orb_q_movie = self.snapshot_dir + "/orbital_elements_pericenter_distance_q.mp4"
            self.make_movies(orb_q_path, orb_q_movie, consec = consec)

# 5c
    def mm_sma(self, movie = 1, plot_kept = False, consec = False):
         self.make_movies_semi_major_axes(plot_kept = plot_kept, consec = consec)

    def make_movies_semi_major_axes(self, movie = 1, plot_kept = False, consec = False):
         """ Plot the evolution of the semi-major axes over time """
         ############## NEEDS A LOT OF CLEAN UP, DOESN'T EVEN WORK ####################
         old_sm_axes = quantities.AdaptingVectorQuantity()
         
         init_two_stars = self.initial_state[0:2]
         planetesimals = self.initial_state[2:]
         
         ### Passing Star (Transferred) ###
                 
         for i, (bodies, passing) in enumerate(zip(self.bodies.history, self.passing)):
             if (self.mask[i]):
                old_sm_axes = quantities.AdaptingVectorQuantity()
                new_sm_axes_p = quantities.AdaptingVectorQuantity()
                eccentricities = []

                time = self.times[i]
                two_stars = bodies[0:2]

                # New Values
                
                s_planetesimals = self.shift_to_star(passing, two_stars, star_choice = 1)
                s_pls = s_planetesimals # make new name later
                for b in s_pls:
                    a, ecc, T = ps.orbital_parameters(b.position, b.velocity, self.info.m1)
                    new_sm_axes_p.append(a) # This is wrong (a is an array) <--- How to handle array?
                    eccentricities.append(ecc)
                #new_sm_axes_p = a
                print "New"
                print new_sm_axes_p
                print
                
                # Note, using the mask isn't as helpful as it could be
                # The sm-axes are still re-calculated. This is bad. They should be stored in a file

                # Old Values
                selected = planetesimals.select(lambda x : x in passing.key, ["key"])
                old_planetesimals = self.shift_to_star(selected, init_two_stars, star_choice = 0)
                old_p = old_planetesimals
                for b in old_p:
                    a, ecc, T = ps.orbital_parameters(b.position, b.velocity, self.info.m0)
                    old_sm_axes.append(a)

                print "Old"
                print old_sm_axes
                print

                if len(new_sm_axes_p) > 0:
                   print "Iteration ", i
                   ps.make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_p, colorcode = eccentricities, \
                               movie = 1, color = PlotFactory.c_PASSING, i = i, time = time) # (1) A
                   ps.make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_p - old_sm_axes, \
                               movie = 1, color = PlotFactory.c_PASSING, i = i, time = time, delta = 1) # (2) Delta A
                   # Add (3) Ecc
                   ps.make_sma_evolution_scatter_plot(old_sm_axes, eccentricities, \
                               movie = 1, color = PlotFactory.c_PASSING, i = i, time = time, eccentric = 1) # (2) Delta A
                            
         if movie == 1:
            # (4) SEMI-MAJOR AXES (transferred)
            sma_p_path = self.snapshot_dir + "/sma_transfer/sma_evolution_%3d.png"
            sma_p_movie = self.snapshot_dir + "/transferred_sma.mp4"
            self.make_movies(sma_p_path, sma_p_movie, consec = consec)

            # (5) DELTA SEMI-MAJOR AXES (transferred)
            d_sma_p_path = self.snapshot_dir + "/sma_transfer_delta/sma_evolution_%3d.png"
            d_sma_p_movie = self.snapshot_dir + "/transferred_delta_sma.mp4"
            self.make_movies(d_sma_p_path, d_sma_p_movie, consec = consec)
         
            ### Central Star (Kept) ###

         if (plot_kept):
            for i, (bodies, central) in enumerate(zip(self.bodies.history, self.central)):
               if (self.mask[i]):
                   old_sm_axes = quantities.AdaptingVectorQuantity()
                   new_sm_axes_c = quantities.AdaptingVectorQuantity()
   
                   time = self.times[i]
                   two_stars = bodies[0:2]
                   
                   s_pls = self.shift_to_star(central, two_stars, star_choice = 0) # Change Name of this later
                   for b in s_pls:
                      a, ecc, T = ps.orbital_parameters(b.position, b.velocity, self.info.m0)
                      new_sm_axes_c.append(a)
                   #new_sm_axes_c = a
   
                   selected = planetesimals.select(lambda x : x in central.key, ["key"])
                   old_planetesimals = self.shift_to_star(selected, init_two_stars, star_choice = 0)
                   old_p = old_planetesimals
                   for b in old_p:
                       a, ecc, T = ps.orbital_parameters(b.position, b.velocity, self.info.m0)
                       old_sm_axes.append(a)
               
                   if len(new_sm_axes_c) > 0:
                       print "Iteration ", i
                       ps.make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_c, \
                                  movie = 1, color = PlotFactory.c_CENTRAL, i = i, time = time) # (4) A
                       ps.make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_c - old_sm_axes, \
                                  movie = 1, color = PlotFactory.c_CENTRAL, i = i, time = time, delta = 1) # (5) Delta A
                       # Add (6) Ecc
                       ps.make_sma_evolution_scatter_plot(old_sm_axes, eccentricities, \
                               movie = 1, color = PlotFactory.c_CENTRAL, i = i, time = time, eccentric = 1)
                               
            if movie == 1:
               # (6) SEMI-MAJOR AXES (kept)
               sma_c_path = self.snapshot_dir + "/sma_kept/sma_evolution_%3d.png"
               sma_c_movie = self.snapshot_dir + "/kept_sma.mp4"
               self.make_movies(sma_c_path, sma_c_movie, consec = consec)
            
               # (7) DELTA SEMI-MAJOR AXES (kept)
               d_sma_c_path = self.snapshot_dir + "/sma_kept_delta/sma_evolution_%3d.png"
               d_sma_c_movie = self.snapshot_dir + "/kept_delta_sma.mp4"
               self.make_movies(d_sma_c_path, d_sma_c_movie, consec = consec)

# 5d

    def count(self):
         """ shortcut for initial and final count """
         self.make_end_count()

    def make_end_count(self, plot = 0):
         """ print the numbers of the final bounds (along with anything else that has been calculated)"""
         self.tmp_mask[1] = 1
         self.tmp_mask[-1] = 1
         self.sort()

         self.simple_count(plot = plot)

    def simple_count(self, plot = 0):
         """ shortcut for make_simple_count """
         self.make_simple_count(plot = plot)

    def make_simple_count(self, plot = 0):
         """ print the numbers of the selected bounds wih no plot """

         print "Applying Existing Mask. To further modify the mask, use set_mask(ranges) and sort()"

         ############### FIXING THIS ###############

         counts_unbound = [len(x) for x in self.unbound if isinstance(x, iterable)]
         counts_central = [len(x) for x in self.central if isinstance(x, iterable)]
         counts_passing = [len(x) for x in self.passing if isinstance(x, iterable)]
         counts_both = [len(x) for x in self.both if isinstance(x, iterable)]

         times = self.compress(self.times, self.mask)

         ps.make_count_curve(times, counts_unbound, counts_central, \
                                      counts_passing, counts_both, plot = plot)

    def counts(self):
         self.sort()
         self.make_count_curve()

    def make_count_curve(self):
         """ plot number of each type of planetesimal over time """
         # Change this to make plot anyway
         counts_unbound = [len(x) for x in self.unbound if isinstance(x, iterable)]
         counts_central = [len(x) for x in self.central if isinstance(x, iterable)]
         counts_passing = [len(x) for x in self.passing if isinstance(x, iterable)]
         counts_both = [len(x) for x in self.both if isinstance(x, iterable)]

         times = self.compress(self.times, self.mask)
         
         ps.make_count_curve(times, counts_unbound, counts_central, \
                                      counts_passing, counts_both)
         """
         else:
             print "cannot plot until everything has been sorted"
             print "instead, doing make_simple_count"
             self.make_simple_count()
         """

# 5e

    ################ THIS WORKS #################

    def hist(self):
         self.make_histograms()

    def make_histograms(self):
         """ 
         plots histograms of the final sma results with respect 
         to the initial sma values of the particles
         
         sma = semi-major axis
         """

         ###### RE-FORMAT THIS ######

         two_stars = self.initial_state[0:2]
         
         transferred_sma = quantities.AdaptingVectorQuantity()
         keys_passing = self.passing[-1].key
         transferred = self.initial_state.select(lambda x : x in keys_passing , ["key"])

         for b in transferred:
             s_b = self.shift_to_star(b, two_stars, star_choice = 0)
             a, ecc, T = ps.orbital_parameters(s_b.position, s_b.velocity, self.info.m0)
             transferred_sma.append(a)
             
         ejected_sma = quantities.AdaptingVectorQuantity()
         keys_unbound = self.unbound[-1].key
         ejected = self.initial_state.select(lambda x : x in keys_unbound , ["key"])

         for b in ejected:
             s_b = self.shift_to_star(b, two_stars, star_choice = 0)
             a, ecc, T = ps.orbital_parameters(s_b.position, s_b.velocity, self.info.m0)
             ejected_sma.append(a)

         ###### RE-FORMAT THIS ######

         global_max = max( len(transferred_sma), len(ejected_sma) )
         
         name = "/plot_transferred_semimajor_axes_histogram.png"
         ps.make_sm_axis_histogram(transferred_sma, name, PlotFactory.c_PASSING,
                                         step_size = 5, max_count = global_max)

         name = "/plot_transferred_semimajor_axes_histogram_cum.png"
         ps.make_sm_axis_histogram(transferred_sma, name, PlotFactory.c_PASSING,
                                         cum = True, max_count = global_max)

         name = "/plot_ejected_semimajor_axes_histogram.png"
         ps.make_sm_axis_histogram(ejected_sma, name, PlotFactory.c_UNBOUND,
                                         step_size = 5, max_count = global_max)

         name = "/plot_ejected_semimajor_axes_histogram_cum.png"
         ps.make_sm_axis_histogram(ejected_sma, name, PlotFactory.c_UNBOUND,
                                         cum = True, max_count = global_max)
         

#####################################################################################################
# 6 (Compare to another factory?)


