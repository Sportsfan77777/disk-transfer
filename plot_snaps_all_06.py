"""
Plots some .hdf5 files
Format? Use .history from read_set_from_file

Centered on Passing Star
***Really depends on input, but the size of the plot here makes more sense if center on passing star***

Example Plots:
x-y scatter plots of multiple snapshots of two stars, one with a disk
Each snapshot is saved to a different formatted_number .png file


#######################################################################
#######################################################################
################   Corrects a format with an error   ##################
#######################################################################
#######################################################################

"""

""" 
Note: The whole thing has only been tested with the disk star as the passing star 
It's not really set up to do the reverse, and the whole thing is kind of messy
"""

# Should I make this into a class? probably...

import os
import subprocess
import numpy

import matplotlib
from optparse import OptionParser
from matplotlib import pyplot,gridspec

from amuse.units import units
from amuse.units import constants
from amuse.units import quantities
from amuse.datamodel import Particles
from amuse.io import read_set_from_file

c0 = 'fffffff' # Either white or black (not sure)
c1 = '#e41a1c' # Red
c2 = '#377eb8' # Light Blue
c3 = '#228b22' # Forest Green
c4 = '#8a2be2' # Purple
c5 = '#ffa500' # Orange
colors = [c0, c1, c2, c3, c4, c5] # Note: colors[0] is a filler
c_STARS = 1
c_CENTRAL = 2 # Originally had disk
c_PASSING = 3 # Passing by disk
c_BOTH = 4
c_UNBOUND = 5

snapshot_dir = None

######################################################################################################
# -9 (FIXED)

def set_snapshot_dir(snap_dir):
  """
  Updates the snapshot_dir
  Warning: This is dangerous and stupid. A class structure would probably be better.
  """
  global snapshot_dir
  snapshot_dir = snap_dir

  try:
     os.mkdir(snapshot_dir)
  except:
     print "\t(" , snapshot_dir, "already exists)"

######################################################################################################
# -8 (FIXED)

def parse_file(f1):
  """
  Parses two files, each using a different coordinate system (one for each star)
  Note: Data stored in a different format will require a different function

  f1 = coordinates with respect to center of mass
  snapshot_dir = directory of choice
  """
  #print " ** plotting: "

  bodies = read_set_from_file(snapshot_dir + "/" + f1, 'hdf5')

  return bodies

######################################################################################################
# -7 (Problems Unknown)

def init_fig():
    """
    Initializes a basic figure
    Note: I fear that this is not general enough (The histograms look weird)
    """
    fig = pyplot.figure(figsize=(6,6))

    """
    lm = 0.12
    rm = 0.99
    tm = 0.95
    bm = 0.1
    gs1 = gridspec.GridSpec(1, 1)
    gs1.update(left=lm, right=rm, top=tm, bottom=bm)
    xy_plane = pyplot.subplot(gs1[:, :])
    
    return fig, xy_plane
    """
    return fig

######################################################################################################
# -5 (Helper)

def orbital_parameters(rel_position, rel_velocity, total_mass):
    """
    Calculates three basic orbital parameters: (a = semimajor axis, e = eccentricity, T = period)
    """
    separation = rel_position.lengths()
    speed_squared = rel_velocity.lengths_squared()

    # a
    semimajor_axis = (constants.G * total_mass * separation / (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)

    # ecc
    eccentricity = numpy.sqrt( 1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=-1) / (constants.G * total_mass * semimajor_axis))

    # T
    period = (2 * numpy.pi * semimajor_axis**1.5 / (numpy.sqrt(constants.G * total_mass))).as_quantity_in(units.yr)
    
    return semimajor_axis, eccentricity, period

######################################################################################################
# -5 (FIXED)

def orbital_parameter_sorted_sets(planetesimals, two_stars):
    """
    Returns four sets of particles, sorted by their 'bounds' to each of the two stars
    """
    unbound = Particles(0)
    one = Particles(0)
    two = Particles(0)
    both = Particles(0)

    # two_stars can't be three stars or one star (hence the name)
    if len(two_stars) != 2:
        print "This is not really a Type Error, but whatever..."
        raise TypeError  # Replace with error message or something

    coor = []
    
    copy0 = planetesimals.copy_to_memory()
    copy1 = planetesimals.copy_to_memory()
    coor.append(copy0)
    coor.append(copy1)

    coor[0].position -= two_stars[0].position
    coor[1].position -= two_stars[1].position

    coor[0].velocity -= two_stars[0].velocity
    coor[1].velocity -= two_stars[1].velocity

    for pl, pl_0, pl_1 in zip(planetesimals, coor[0], coor[1]):
        rel_pos = []
        rel_pos.append(pl_0.position - two_stars[0].position) # 505, 507 fix
        rel_pos.append(pl_1.position - two_stars[0].position)

        rel_vel = []
        rel_vel.append(pl_0.velocity - two_stars[0].velocity)
        rel_vel.append(pl_1.velocity - two_stars[0].velocity)

        a_a, ecc_a, T_a = orbital_parameters(rel_pos[0], rel_vel[0], two_stars[0].mass)
        a_b, ecc_b, T_b = orbital_parameters(rel_pos[1], rel_vel[1], two_stars[1].mass)

        if (0 <= ecc_a < 1 and 0 <= ecc_b < 1):
             both.add_particle(pl)
        elif (0 <= ecc_a < 1):
             one.add_particle(pl)
        elif (0 <= ecc_b < 1):
             two.add_particle(pl)
        else:
             unbound.add_particle(pl)

    return unbound, one, two, both

######################################################################################################
# -4 (FIXED)

def make_colorcoded_position_scatter_plot(i, time, two_stars, unbound, central, passing, both,
                                             min_a = 40, max_a = 100, star_choice = 1):
    """
    Make Movies!
    Particles are color-coded based on how they are bound

    Input:
    ?????????????????

    Ouput:
    return lengths of each array (this is not needed)
    """

    try:
      os.mkdir(snapshot_dir + "/spdd")
    except:
      print "\t(" , snapshot_dir, "already exists)"

    try:
      os.mkdir(snapshot_dir + "/zoom_spdd")
    except:
      print "\t(" , snapshot_dir, "already exists)"

    print i

    fig = init_fig()

    pyplot.axes().set_aspect('equal') # Will this come back to haunt me later? Not sure
    
    lim=450
    c = 2.8
    pyplot.xlim(-c*lim, c*lim)
    pyplot.ylim(-c*lim, c*lim)

    # Set Coordinate System

    # First, account for stupid inputs
    if (star_choice != 0 and star_choice != 1):
        star_choice = 1

    c_unbound = unbound.copy() # You must copy to avoid trying to write to history
    c_central = central.copy() # However, the 'c_' prefix is probably unnecessary (investigate later...)
    c_passing = passing.copy()
    c_both = both.copy()
    c_two_stars = two_stars.copy()

    # Then, set the star_choice to the center
    # And, plot planetesimals with star_choice at the center

    if (len(unbound) > 0):
        c_unbound.position -= two_stars[star_choice].position + two_stars[0].position
        xd_unbound = c_unbound.x.value_in(units.AU) # array of x-cor for disk particles -- unbound (ORANGE)
        yd_unbound = c_unbound.y.value_in(units.AU)
        pyplot.scatter(xd_unbound, yd_unbound, c=colors[c_UNBOUND], edgecolors = "none", alpha = 0.4, zorder=2)

    if (len(central) > 0):
        c_central.position -= two_stars[star_choice].position + two_stars[0].position
        xd_central = c_central.x.value_in(units.AU) # array of x-cor-- bound to central (moving plot) (BLUE)
        yd_central = c_central.y.value_in(units.AU)
        pyplot.scatter(xd_central, yd_central, c=colors[c_CENTRAL], edgecolors = "none", alpha = 0.4, zorder=1)

    if (len(passing) > 0):
        c_passing.position -= two_stars[star_choice].position + two_stars[0].position
        xd_passing = c_passing.x.value_in(units.AU) # array of x-cor -- bound to passing (center plot) (GREEN)
        yd_passing = c_passing.y.value_in(units.AU)
        pyplot.scatter(xd_passing, yd_passing, c=colors[c_PASSING], edgecolors = "none", alpha = 0.4, zorder=3)

    if (len(both) > 0):
        c_both.position    -= two_stars[star_choice].position + two_stars[0].position
        xd_both = c_both.x.value_in(units.AU) # array of x-cor for disk particles -- bound to both (PURPLE)
        yd_both = c_both.y.value_in(units.AU)
        pyplot.scatter(xd_both, yd_both, c=colors[c_BOTH], edgecolors = "none", alpha = 0.4, zorder=4)

    time_yr_str = "{0:.1f} yr".format(time.value_in(units.yr))
    
    plot_i = snapshot_dir+"/spdd/spdd_{0:03d}.png".format(i)
    plot_j = snapshot_dir+"/zoom_spdd/zoom_spdd_{0:03d}.png".format(i)
    print " \t", time_yr_str, " -> ", plot_i
    print " \t", time_yr_str, " -> ", plot_j
    print "Unbound:", len(unbound), "| Transferred:", len(passing), "| Kept:", len(central), "| Both:", len(both)

    # Now, plot stars with star_choice at the center

    c_two_stars[1 - star_choice].position -= two_stars[star_choice].position
    c_two_stars[star_choice].position     -= two_stars[star_choice].position # careful not to overwrite too early

    xb = c_two_stars.x.value_in(units.AU)
    yb = c_two_stars.y.value_in(units.AU)
    
    pyplot.scatter(xb, yb, c=colors[c_STARS], lw=0.5, zorder=0)
    pyplot.xlabel('x [AU]')
    pyplot.ylabel('y [AU]')
    

    #pyplot.text(0.5, 0.98, time_yr_str, 
    #            horizontalalignment='center', verticalalignment='bottom', 
    #            transform=xy_plane.transAxes)
    pyplot.title(time_yr_str)

    # Inner and Outer Radii of Initial Disks
    c_1 = pyplot.Circle((xb[0],yb[0]), radius = min_a, ec='black', fc = "none", zorder = 10)
    c_2 = pyplot.Circle((xb[0],yb[0]), radius = max_a, ec='black', fc = "none", zorder = 10)

    c_3 = pyplot.Circle((xb[1],yb[1]), radius = min_a, ec='black', fc = "none", zorder = 10)
    c_4 = pyplot.Circle((xb[1],yb[1]), radius = max_a, ec='black', fc = "none", zorder = 10)

    pyplot.axes().add_artist(c_1)
    pyplot.axes().add_artist(c_2)
    pyplot.axes().add_artist(c_3)
    pyplot.axes().add_artist(c_4)
    
    #pyplot.tight_layout()
    pyplot.savefig(plot_i)

    c = 0.8
    pyplot.xlim(-c*lim, c*lim)
    pyplot.ylim(-c*lim, c*lim)

    pyplot.savefig(plot_j)

    pyplot.cla() # Clears Figure

    pyplot.close() # Closes Figure (both are probably not necessary)

    return len(unbound), len(central), len(passing), len(both)

######################################################################################################
# -3 (FIXED)

def make_sma_evolution_scatter_plot(initial_sm_axes, final_sm_axes,
                                      delta = 0, movie = 0, color = 3, i = -1, time = None):
   """
   Makes a scatter plot of the initial 'a' to the transferred 'a' around the new star
   Options:
   color = 3 : transferred
   color = 2 : kept

   Note: This should be modified to taken in additional parameters (sm-axis of inner and outermost particles)
   """

   if color == 3:
      directory = snapshot_dir + "/sma_transfer"
   elif color == 2:
      directory = snapshot_dir + "/sma_kept"

   if delta == 1:
      directory += "_delta"

   try:
      os.mkdir(directory)
   except:
      print "\t(" , directory, "already exists)"
 
   # Add disk range parameter initialized from different arrays that contain all the particles
   # Change function name

   fig = init_fig()

   pyplot.scatter(initial_sm_axes.value_in(units.AU), final_sm_axes.value_in(units.AU), \
                            c=colors[3], lw=0.5, zorder=0)

   min_x = 40
   mx_x = 100
   if movie == 1:
      mx_y = 1000
   else:
      mx_y = max(final_sm_axes.value_in(units.AU))
   pyplot.xlim(min_x / 1.02, mx_x * 1.02)
   pyplot.ylim(0, mx_y * 1.05)

   pyplot.xlabel('initial sm-axis [AU]')
   pyplot.ylabel('final sm-axis [AU]')

   # Plot 'y' vs 'x'
   pyplot.plot([min_x / 1.02, mx_x * 1.02], [min_x / 1.02, mx_x * 1.02], linestyle='--', color ='black', zorder = 1)

   pyplot.plot([40, 40], [0, 1], linestyle='--', color = 'black', zorder = 1)
   pyplot.plot([100, 100], [0, 1], linestyle='--', color = 'black', zorder = 1)

   time_yr_str = ""
   if (time is not None):
       time_yr_str = "{0:.1f} yr".format(time.value_in(units.yr))

   title_str = "Evolution of Semi-Major Axes for Transferred Particles\n" + time_yr_str
   #pyplot.text(0.5, 0.999, title_str, 
   #             horizontalalignment='center', verticalalignment='bottom', 
   #             transform=xy_plane.transAxes)
   pyplot.title(title_str)

   plot_sma = snapshot_dir+"/sma_evolution_{0:03d}.png".format(i)
   pyplot.savefig(plot_sma)

   pyplot.cla()
   
######################################################################################################
# -2 (FIXED)

def make_orbital_elements_scatter_plot(bodies, mass, i = -1, min_a = 40, max_a = 100):
   """
   Given a set of particles and the mass of the appropriate, user-selected central object
   plot the eccentricities against the semi-major axes ('ecc' vs 'a')
   """

   try:
      os.mkdir(snapshot_dir + "/orb")
   except:
      print "\t(" , snapshot_dir, "already exists)"

   semimajor_axes = []
   eccentricities = []

   for b in bodies:
      a, ecc, T = orbital_parameters(b.position, b.velocity, mass)
      semimajor_axes.append(a.value_in(units.AU))
      eccentricities.append(ecc)

   fig = init_fig()

   pyplot.scatter(semimajor_axes, eccentricities, c=colors[3], lw=0.5, zorder=0)

   min_x = min(semimajor_axes)
   mx_x = max(semimajor_axes)

   pyplot.xlim(min([min_a / 1.05, min_x / 1.05]), max([mx_x * 1.02, 1000])) # sm-axes
   pyplot.ylim(0, 1) # ecc

   # For comparison, plot initial disk limits for the other star
   pyplot.plot([min_a, min_a], [0, 1], linestyle='--', color = 'black', zorder = 1)
   pyplot.plot([max_a, max_a], [0, 1], linestyle='--', color = 'black', zorder = 1)

   pyplot.xlabel('final sm-axis [AU]')
   pyplot.ylabel('final eccentricity')

   title_str = "Orbital Elements for Transferred Particles"
   # There is something wrong with the placement (unknown)
   #pyplot.text(0.5, 0.999, title_str, 
   #             horizontalalignment='center', verticalalignment='bottom', 
   #             transform=xy_plane.transAxes)
   pyplot.title(title_str)

   plot_ecc = ""
   if (not (i == -1)):
      plot_ecc = snapshot_dir+"/orb/orb_{0:03d}.png".format(i)
   else:
      plot_ecc = snapshot_dir+"/plot_transferred_new_orbital_elements.png"

   pyplot.savefig(plot_ecc)
   pyplot.cla()

######################################################################################################
# -1 (FIXED)

def make_sm_axis_histogram(sm_axes, filename, color, \
               cum = False, min_a = 40, max_a = 100, max_count = None, step_size = 1):
   """
   Given a set of sm_axes, make a histogram

   Specific Options:
   color = 3 : Transferred Particles
   color = 5 : Unbound Particles
   """

   fig = init_fig()
   axes = fig.add_subplot(111)

   title_str = ""
   if (color == 3):
      # 3, Green, Transferred
      title_str = "Initial Semimajor-Axes of Transferred Particles"
   elif (color == 5):
      # 5, Orange, Unbound
      title_str = "Initial Semimajor-Axes of Unbound Particles"

   if (not (0 <= color <= 5)):
      color = 0

   min_x = min_a
   max_x = max_a
   step = step_size
   bin_array = range(min_x, max_x + step, step)

   # the histogram of the data
   n, bins, patches = axes.hist(sm_axes.value_in(units.AU), bin_array, facecolor = colors[color], \
                                        alpha=0.75, cumulative = cum)

   # There is something wrong with the placement (unknown)
   #pyplot.text(0.5, 0.999, title_str, 
   #             horizontalalignment='center', verticalalignment='bottom', 
   #             transform=xy_plane.transAxes)
   pyplot.title(title_str)

   pyplot.xlabel('initial sm-axis [AU] bins')
   pyplot.ylabel('count')

   pyplot.xlim(min_x, max_x)
   if max_count is not None:
      pyplot.ylim(0, max_count)
   
   pyplot.savefig(snapshot_dir + filename)
   
   pyplot.cla()

######################################################################################################
# 0 (FIXED)

# Consider changing input? Or deal with it in factory? Both would work, but the latter is probably better

def make_count_curve(times, counts_unbound, counts_central, counts_passing, counts_both, plot = 1):

   t = times.value_in(1e3 * units.yr)

   count_file = snapshot_dir + "/particle_final_count.txt"

   if plot == 1:

      pyplot.xlabel('time (in kyr)')
      pyplot.ylabel('counts')

      pyplot.axes().set_yscale('log') # FIX THIS

      pyplot.xlim(0, max(t))
      pyplot.ylim(0, 1000)

      pyplot.plot(t, counts_unbound, color = colors[c_UNBOUND])
      pyplot.plot(t, counts_passing, color = colors[c_PASSING])
      pyplot.plot(t, counts_central, color = colors[c_CENTRAL])
      pyplot.plot(t, counts_both, color = colors[c_BOTH])

      # Add title?
      title_str = "Particle Counts Over Time"
      pyplot.title(title_str)
   
      plot_count_curve = snapshot_dir + "/particle_counts_over_time.png"
      pyplot.savefig(plot_count_curve)
      pyplot.cla()

      # Write to File
      count_file = snapshot_dir + "/particle_count_table.txt"

   f = open(count_file, 'a')

   # set up to be a string of a certain length
   t_str = '{:^18}'.format("Time (kyr)") # centered
   u_str = '{:^18}'.format("Count (Unbound)")
   k_str = '{:^18}'.format("Count (Kept)") # Central
   tr_str = '{:^18}'.format("Count (Transfer)") # Passing
   b_str = '{:^18}'.format("Count (Both)")

   # Header
   f.write(t_str + "|" + u_str + "|" + k_str + "|" + tr_str + "|" + b_str + "\n")
   
   ts = t
   for (t, a, b, c, d) in zip(ts, counts_unbound, counts_central, counts_passing, counts_both):
      t_str = '{:^18}'.format("%.1f" % t) # centered
      u_str = '{:^18}'.format(str(a))
      k_str = '{:^18}'.format(str(b)) # Central
      tr_str = '{:^18}'.format(str(c)) # Passing
      b_str = '{:^18}'.format(str(d))
      f.write(t_str + "|" + u_str + "|" + k_str + "|" + tr_str + "|" + b_str + "\n")

   f.close()

######################################################################################################
# 1 (FIXED)

def plot_distances(time, bodies):
   """ star distances vs time """
   distances = quantities.AdaptingVectorQuantity()

   fig = init_fig()

   for b in bodies.history:
      two_stars = b[0:2]
      rel_position = two_stars[0].position - two_stars[1].position
      separation = rel_position.lengths()
      distances.append(separation)

   pyplot.xlabel('time (in kyr)')
   pyplot.ylabel('distance (in AU)')

   pyplot.plot(time.value_in(1000 * units.yr), distances.value_in(units.AU), c = colors[c_STARS])

   pyplot.title("Star Separation Over Time")

   fn = snapshot_dir + "/star_separation.png"
   pyplot.savefig(fn)
   pyplot.cla()

   return distances


#########################################################################################################
#########################################################################################################
######################################## Make all of the plots ##########################################
#########################################################################################################
#########################################################################################################

# Set up two input options?
# Deprecate this?
# Create a function in the factory class to do this (only)?

def plot_all_snaps(f1, bodies = None):
  """ Deprecated """
  if bodies is None and bodies2 is None:
      bodies, bodies2 = parse_files(f1, snapshot_dir)

  fig = init_fig()

  ############################

  # Save initial and final conditions
  bi_i = 0
  bi_f = 0

  for i, (bi, b2i) in enumerate(zip(bodies.history, bodies2.history)):
    # Save initial conditions (Note: this is the 1st timestep, not the 0th timestep)
    if (i == 1):
       bi_i = b2i
       
    # Save final conditions
    bi_f = bi

  ############################

  mass_zero = bi[0].mass # central mass (w/ disk)
  mass_one = bi[1].mass # passing mass (w/ no disk)

  times = quantities.AdaptingVectorQuantity()

  counts_unbound = []
  counts_passing = []
  counts_central = []
  counts_both = []
  
  # Call .history to get data at each time step
  for i, (bi, b2i) in enumerate(zip(bodies.history, bodies2.history)):
    # Mark time
    time =  bi.collection_attributes.timestamp
    times.append(time)

    # Sort particles by bounds (through orbital element calculations)
    unbound, passing, central, both, k_unbound, k_passing, k_central, k_both = \
              init_orbital_parameter_arrays(bi, b2i)

    # Make Position Plots ('y' vs 'x')
    count_unbound, count_passing, count_central, count_both = \
              make_colorcoded_position_scatter_plot(i, time, unbound, passing, central, both, snapshot_dir)

    counts_unbound.append(count_unbound)
    counts_passing.append(count_passing)
    counts_central.append(count_central)
    counts_both.append(count_both)

    # Make 4 "SM-axis" Scatter Plot(s)
    old_sm_axes = quantities.AdaptingVectorQuantity()
    new_sm_axes_t = quantities.AdaptingVectorQuantity()
    new_sm_axes_k = quantities.AdaptingVectorQuantity()

    for (keyp, keyc) in zip(k_passing, k_central):
       b_old = [x for x in bi_i if x.key == key]
       b_new_t = [x for x in bi if x.keyp == keyp] # t = transfer
       b_new_k = [x for x in bi if x.keyc == keyc] # k = kept

       if (len(b_old) > 1 or len(b_new_t) > 1 or len(b_new_k) > 1):
          print "Is this key copied?", keyp, keyc

       if (len(b_old) and len(b_new) >= 1):
          b_old = b_old[0]
          b_new = b_new[0]

          a1, e1, T1 = orbital_parameters(b_old.position, b_old.velocity, mass_zero)
          a2, e2, T2 = orbital_parameters(b_new.position, b_new.velocity, mass_one)

          old_sm_axes.append(a1)
          new_sm_axes.append(a2)

    make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_t, movie = 1, color = 3, i = i, time = time) # (1) A
    make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_t - old_sm_axes, \
                                                       movie = 1, color = 3, i = i, time = time) # (2) Delta A

    make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_k, movie = 1, color = 2, i = i, time = time) # (3) A
    make_sma_evolution_scatter_plot(old_sm_axes, new_sm_axes_k - old_sm_axes, \
                                                       movie = 1, color = 2, i = i, time = time) # (4) Delta A

    # Make Orbital Elements Plot
    b_new = [x for x in bi if x.key in (set(bi.key) & set(k_passing))]
    make_orbital_elements_scatter_plot(b_new, mass_one, i)

    
  # make movies!

  fps = 6

  # (1) POSITION 
  pos_path = snapshot_dir + "/spdd/spdd_%3d.png"
  pos_movie = snapshot_dir + "/spdd/snaps.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + pos_path + "-vcodec mpeg4 -y " + pos_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (2) ZOOOOMED POSTION
  zoom_path = snapshot_dir + "/zoom_spdd/zoom_spdd_%3d.png"
  zoom_pos_movie = snapshot_dir + "/spdd/snaps_zoom.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + zoom_path + "-vcodec mpeg4 -y " + zoom_pos_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (3) ORBITAL ELEMENTS (ecc vs a)
  orb_path = snapshot_dir + "/orb/orb_%3d.png"
  orb_movie = snapshot_dir + "/orb/orbital_elements.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + orb_path + "-vcodec mpeg4 -y " + orb_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (4) SEMI-MAJOR AXES (transferred)
  sma_t_path = snapshot_dir + "/sma_transfer/sma_evolution_%3d.png"
  sma_t_movie = snapshot_dir + "/sma_transfer/transferred_sma.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + sma_t_path + "-vcodec mpeg4 -y " + sma_t_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (5) SEMI-MAJOR AXES (kept)
  sma_k_path = snapshot_dir + "/sma_kept/sma_evolution_%3d.png"
  sma_k_movie = snapshot_dir + "/sma_kept/kept_sma.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + sma_k_path + "-vcodec mpeg4 -y " + sma_k_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (6) DELTA SEMI-MAJOR AXES (transferred)
  sma_t_path = snapshot_dir + "/sma_transfer_delta/sma_evolution_%3d.png"
  sma_t_movie = snapshot_dir + "/sma_transfer_delta/transferred_sma.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + sma_t_path + "-vcodec mpeg4 -y " + sma_t_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  # (7) DELTA SEMI-MAJOR AXES (kept)
  sma_k_path = snapshot_dir + "/sma_kept_delta/sma_evolution_%3d.png"
  sma_k_movie = snapshot_dir + "/sma_kept_delta/kept_sma.mp4"
  command = "ffmpeg -f image2 -r " + fps + " -i " + sma_k_path + "-vcodec mpeg4 -y " + sma_k_movie
  movie_pipe = subprocess.Popen(command, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

  ###################

  # make plots (that aren't movies)

  # (P1) BOUND COUNTS (all)

  make_count_curve(times, counts_unbound, counts_passing, counts_central, counts_both, snapshot_dir)

  # (P2) HISTOGRAMS (4 of them)
  # Make "Ejected SM-axis" + "Transferred SM-axis" Histograms (normal + cumulative)

  # Note: collapse this into a separate function
  transferred_sma = quantities.AdaptingVectorQuantity()
  for key in k_passing:
    b_old = [x for x in bi_i if x.key == key]

    if (len(b_old) > 1 or len(b_new) > 1):
       print "Is this key copied?", key

    if (len(b_old) >= 1):
       b_old = b_old[0]

       a, ecc, T = orbital_parameters(b_old.position, b_old.velocity, mass_zero)
       transferred_sma.append(a)

  ejected_sma = quantities.AdaptingVectorQuantity()
  for key in k_unbound:
    b_old = [x for x in bi_i if x.key == key]

    mass_zero = bi[0].mass # central mass

    if (len(b_old) > 1):
       print "Is this key copied?", key

    if (len(b_old) >= 1):
       b_old = b_old[0]

       a, ecc, T = orbital_parameters(b_old.position, b_old.velocity, mass_zero)
       ejected_sma.append(a)

  # add cumulative histograms also (mostly for total count)

  name = "/plot_transferred_semimajor_axes_histogram.png"
  make_sm_axis_histogram(transferred_sma, name, c_PASSING) # 2nd parameter is color (an identifier)

  name = "/plot_transferred_semimajor_axes_histogram_cum.png"
  make_sm_axis_histogram(transferred_sma, name, c_PASSING, cum = True) # 2nd parameter is color (an identifier)

  name = "/plot_ejected_semimajor_axes_histogram.png"
  make_sm_axis_histogram(ejected_sma, name, c_UNBOUND)

  name = "/plot_ejected_semimajor_axes_histogram_cum.png"
  make_sm_axis_histogram(ejected_sma, name, c_UNBOUND, cum = True)

  # write things to files

  # start with pretty much any plot data (either (A) make it neat or (B) use .hdf5 format or (C) both)
    
  return
