import numpy 
import os
import time

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.community.huayno.interface import Huayno
from amuse.community.kepler6.interface import Kepler
#from kepler6.interface import Kepler
from amuse.io import write_set_to_file
from amuse.couple import bridge

import pickle

from disk_ini import get_planetesimals_disk
from orbit_ini_01 import get_orbit_ini
import plot_snaps_all_05 as ps

from factory import PlotFactory as PF ######### All Data Analysis #########

import disk_flyby_nbody_02 as nbody

pi_180 = numpy.pi/180.0

def rm_file(file_name):
  """
  delete the file 'file_name', if it exists
  """  
  if(os.path.isfile(file_name)):
    os.system('rm %s' % (file_name))

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()

class planet_gravity_for_disk(object):
    """
    copy = copycat(base_class, grav_instance, converter)
    derived system, returns copy of grav instance with
    get_gravity_at_point, get_potential_at_point reimplemented in 
    base_class
    -- transforms from the coordinate system of the star--planets system to the
       coordinate system centered at the star (sun) and gives the gravity from the
       planet in this system
       + corrects for the non-inetrtial system of the star (sun) due to planets
    """
    def __init__(self, baseclass, sun_and_planets, converter, center = 0):
        self.center = center

        self.baseclass = baseclass
        self.converter = converter
        self.sun_and_planets = sun_and_planets
        
        sun = sun_and_planets[center:center+1]
        not_sun = sun_and_planets - sun
        
        # initializing the gravity code
        self.instance=self.baseclass(self.converter)
        self.instance.initialize_code()
        self.instance.particles.add_particles(not_sun) # add particles that are not the center particle
        
        #self.channel_from_planets_to_instance = planets.new_channel_to(self.instance.particles, ["x", "y", "z", "vx", "vy", "vz"])
          
    def get_gravity_at_point(self,radius,x,y,z):
        
        # set planets heliocentric
        local_sun_and_planets = self.sun_and_planets.copy()
        sun = local_sun_and_planets[self.center:self.center + 1]
        planets = local_sun_and_planets - sun 
        self.instance.particles.position = planets.position - sun.position
        self.instance.particles.velocity = planets.velocity - sun.velocity
        
        # acceleration from the planets
        ax,ay,az = self.instance.get_gravity_at_point(radius,x,y,z)
        
        # acceleration of the central particle (sun) due to planets
        ax0,ay0,az0 = self.instance.get_gravity_at_point(radius,(0.0|units.AU),(0.0|units.AU),(0.0|units.AU))
        # ^^^^ sun not at (0,0,0)?

        """
        print "Sun:", sun.position
        print "Planets:", planets.position
        print "ax", ax[0:3]
        print "ax0", ax0[0:3]
        print "ax_diff", str(ax[0:3] - ax0[0:3])
        print "ay", ay[0:3]
        print "ay0", ay0[0:3]
        print "ay_diff", str(ay[0:3] - ay0[0:3])
        print
        """        

        return ax-ax0,ay-ay0,az-az0

    def get_potential_at_point(self,radius,x,y,z):
        local_sun_and_planets = self.sun_and_planets.copy()
        sun = local_sun_and_planets[self.center:self.center + 1]
        planets = local_sun_and_planets - sun 
        self.instance.particles.position = planets.position - sun.position
        self.instance.particles.velocity = planets.velocity - sun.velocity
        phi = self.instance.get_potential_at_point(radius,x,y,z)
        phi0 = self.instance.get_potential_at_point(radius,(0.0|units.AU),(0.0|units.AU),(0.0|units.AU))
        # ^^^^ sun not at (0,0,0)?

        """
        print "Sun:", sun.position
        print "Planets:", planets.position
        print "Phi:", phi
        print "Phi0:", phi0
        """

        return phi-phi0
   
    def stop(self):
        self.instance.stop()

def initialize_kepler(planetesimals, converter, file_redir, center = 0):
  if file_redir is None:
    planetesimals_gravity = Kepler(converter,channel_type="sockets")
    #print planetesimals_gravity
  elif file_redir=="0":
    planetesimals_gravity = Kepler(converter,channel_type="sockets", redirection="none")
  else:
    planetesimals_gravity = Kepler(converter,channel_type="sockets", redirection="file", redirect_file=file_redir)
    
  sun = stars[center:center+1] # weird format: because input must be array

  planetesimals_gravity.central_particle.add_particles(sun) ####### <<<---- Set central particle here #######
  planetesimals_gravity.central_particle.position = (0.0, 0.0, 0.0) | units.AU
  planetesimals_gravity.central_particle.velocity = (0.0, 0.0, 0.0) | units.kms
  planetesimals_gravity.orbiters.add_particles(planetesimals)
  planetesimals_gravity.commit_particles()
  planetesimals_gravity.particles = planetesimals_gravity.orbiters # make sure the kicks are done on orbiters only

  return planetesimals_gravity

def initialize_huayno(bodies, converter, huayno_eta):
  stars_gravity = Huayno(converter,channel_type="sockets")
  stars_gravity.particles.add_particles(bodies) # either bodies or just stars
  stars_gravity.commit_particles()
  #stars_gravity.set_timestep_parameter(huayno_eta)
  stars_gravity.parameters.timestep_parameter = huayno_eta
  stars_gravity.set_inttype_parameter(12) # CC_KEPLER

  return stars_gravity
        
def integrate_disk_flyby(stars, planetesimals, t_end, t_peri, nb_end, n_steps,
                         r_step, snap_dir, file_out, file_redir, bridge_dt, huayno_eta):
  """
  stars --- Star bodies
  planetesimals --- Non-star bodies
  t_end --- Total Integration Time
  n_steps --- Number of Snapshots
  r_step ---
  snap_dir --- Snapshot Directory
  file_out --- Output File (cannot be None)
  file_redir --- Redirection File (can be None)
  bridge_dt --- Bridge Timestep
  huayno_eta --- Huayno Timestep
  """

  converter=nbody_system.nbody_to_si(1|units.MSun, 1|units.AU) # set nbody scales
  
  # Initialize Huayno (N-body Solver)
  stars_gravity = initialize_huayno(stars, converter, huayno_eta)
  
  #print file_redir
  
  t0 = time.time()

  """
  Notes: 
  (1) Careful with coordinate system
  (2) Different timesteps
  (3) 
  """

  # Calculate t_end_a, t_end_b

  rm_file(snap_dir + "/" + file_out)

  # (Part A)
  # planetesimals_gravity = Kepler_0 [0, t_end_a]
  planetesimals_gravity = initialize_kepler(planetesimals, converter, file_redir, center = 0)

  # Calculate number of steps for each integration
  t_0 = 0.0 | units.yr
  dt = (t_end - t_0) / float(n_steps)

  k_end = 0.84
  t_end_a = k_end * abs(t_peri)
  n_steps_a = round(t_end_a / dt) # number of steps for (Part A)

  t_end_b = t_end_a + (nb_end - k_end) * abs(t_peri) #### <---- fix this later
  n_steps_b = round ( (t_end_b - t_end_a) / dt) # number of steps for (Part B)

  n_steps_c = n_steps - (n_steps_a + n_steps_b)

  print "t_peri", abs(t_peri.value_in(units.yr)), "t_end_a", t_end_a, "n_steps_a", n_steps_a, "n_steps_c", n_steps_c, "dt", dt.value_in(units.yr)
  print

  stars, planetesimals = evolve_disk_flyby(stars, planetesimals, stars_gravity, planetesimals_gravity, 
                             t_0, t_end_a, n_steps_a, converter, snap_dir, file_out, r_step, bridge_dt, center = 0)

  print stars
  print
  print planetesimals
  print

  # (Part B)
  # planetesimals_gravity = Huayno [t_end_a, t_end_b]
  bodies = ParticlesSuperset([stars, planetesimals])
  gravity = initialize_huayno(bodies, converter, huayno_eta)
  stars, planetesimals = evolve_disk_flyby_together(stars, planetesimals, gravity,
                             t_end_a, t_end_b, n_steps_b, converter, snap_dir, file_out)

  planetesimals.position -= stars[1].position
  planetesimals.velocity -= stars[1].velocity

  print stars
  print
  print planetesimals
  print

  # (Part C)
  # planetesimals_gravity = Kepler_1 [t_end_b, t_end]
  stars_gravity = initialize_huayno(stars, converter, huayno_eta)
  planetesimals_gravity = initialize_kepler(planetesimals, converter, file_redir, center = 1)

  stars, planetesimals = evolve_disk_flyby(stars, planetesimals, stars_gravity, planetesimals_gravity,
                           t_end_b, t_end, n_steps_c, converter, snap_dir, file_out, r_step, bridge_dt, center = 1)

  print stars
  print
  print planetesimals
  print

  t1= time.time()
  dt = t1-t0
  print " ** Performace data: N =", len(stars)+len(planetesimals), "dt=", dt, "s =", dt/60.0, "min"
  
  return

def evolve_disk_flyby_together(stars, planetesimals, gravity, 
                      t_start, t_end, n_steps, converter, snap_dir, file_out):
  bodies = ParticlesSuperset([stars, planetesimals])
  
  channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)
  
  rm_file(file_out)
    
  Etot_init = gravity.kinetic_energy + gravity.potential_energy
  Etot = Etot_init
  
  ps.mkdir(snap_dir)
  
  duration = t_end - t_start
  dt = duration / float(n_steps)
  time = 0.0 | units.yr

  print "Duration:", duration, "t_start", t_start, "t_end", t_end, "dt:", dt
  
  print " ** evolving: t_start = ", t_start.value_in(1000 * units.yr), "t_end = ", t_end.value_in(1000 * units.yr), ", dt = ", dt
  print " \t", "time", "\t\t\t", "E", "\t\t", "dE"

  stdout = (file_out.split('.'))[0]
  stdout += '.txt'
  f = open(snap_dir + "/" + stdout, 'a')


  while time<=duration:
    gravity.evolve_model(time)
    channel_from_gr_to_framework.copy()
    
    bodies.collection_attributes.timestamp = time + t_start

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot

    nb_E = converter.to_nbody(Etot)
    #nb_J = converter.nbody_length ** 2 * units.nbody_mass * units.nbody_time ** -2 # not supposed to work

    # A formatted string would work better than tabs. (Tabs never work)
    line = " \t" + str((time + t_start).value_in(units.yr)) + "\t" + str(nb_E) + "\t" + str(dE/Etot_init)
    print line
    
    f.write(line + "\n")

    # Write coordinates in Center of Mass frame

    # Move Stars to CoM (initially in CoM)
    #### The stars are already in CoM coordinates ####

    # Move planetesimals to CoM (initially in CoM)
    #### The stars are already in CoM coordinates #### (Does this mean the other method is wrong?)
    
    write_set_to_file(bodies, snap_dir + "/" + file_out, "hdf5")
    
    time += dt
  
  gravity.stop()
  
  return stars, planetesimals
  
def evolve_disk_flyby(stars, planetesimals, stars_gravity, planetesimals_gravity, 
                    t_start, t_end, n_steps, converter, snap_dir, file_out, r_step, bridge_dt, center = 0):

  """
  Input Parameters:
  stars --- Star bodies
  planetesimals --- Non-star bodies
  stars_gravity --- 
  planetesimals_gravity --- 
  t_end --- Total Integration Time
  n_steps --- Number of Snapshots
  converter --- 
  snap_dir --- 
  file_out --- 
  r_step --- 
  bridge_dt --- Bridge Timestep
  """

  bodies = ParticlesSuperset([stars, planetesimals])
  
  channel_from_stars_to_framework = stars_gravity.particles.new_channel_to(stars)
  channel_from_planetesimals_to_framework = planetesimals_gravity.particles.new_channel_to(planetesimals)
  
  pl_gravity = planet_gravity_for_disk(Huayno, stars_gravity.particles, converter, center = center)
  gravity = bridge.Bridge(use_threading=False)
  gravity.add_system(planetesimals_gravity, (pl_gravity,) ) # stars work on planetesimals
  gravity.add_system(stars_gravity, () ) # stars are self aware (evolves itself)
  
  # timestep for BRIDGE relative to the period at r_step (or r_min, if no r_step is given)
  #   BRIDGE_timestep = bridge_dt * P(r_step)
  if r_step is None:
    Porb_min = orbital_period((bodies[2:].position.lengths()).min(), bodies[0].mass)
    print " ** P_min_disk = ", Porb_min.in_(units.yr)
    gravity.timestep = bridge_dt*Porb_min
  else:
    Porb_min = orbital_period(r_step, bodies[0].mass)
    print " ** P_rin_disk = ", Porb_min.in_(units.yr)
    gravity.timestep = bridge_dt*Porb_min
    
  time_step = stars_gravity.get_timestep_parameter()
  print ' ** timesteps: \t stars gravity timestep parameter =', time_step 
  print '\t\t bridge =', gravity.timestep
  
  Etot_init = stars_gravity.kinetic_energy + stars_gravity.potential_energy
  Etot = Etot_init

  ps.mkdir(snap_dir)
  
  duration = t_end - t_start
  dt = duration / float(n_steps)
  time = 0.0 | units.yr

  print "Duration:", duration, "t_start", t_start, "t_end", t_end, "dt:", dt
  
  print " ** evolving: t_start = ", t_start.value_in(1000 * units.yr), "t_end = ", t_end.value_in(1000 * units.yr), ", dt = ", dt
  print " \t", "time", "\t\t\t", "E", "\t\t", "dE"

  stdout = (file_out.split('.'))[0]
  stdout += '.txt'
  f = open(snap_dir + "/" + stdout, 'a')

  # Joules
  J = units.m ** 2 * units.kg * units.s ** -2

  while time<=duration:
    gravity.evolve_model(time)
    channel_from_stars_to_framework.copy()
    channel_from_planetesimals_to_framework.copy()
    
    bodies.collection_attributes.timestamp = time + t_start

    Ekin = stars_gravity.kinetic_energy 
    Epot = stars_gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot

    nb_E = converter.to_nbody(Etot)
    #nb_J = converter.nbody_length ** 2 * units.nbody_mass * units.nbody_time ** -2 # not supposed to work

    # A formatted string would work better than tabs. (Tabs never work)
    line = " \t" + str((time + t_start).value_in(units.yr)) + "\t" + str(nb_E) + "\t" + str(dE/Etot_init)
    print line
    
    f.write(line + "\n")

    # Write coordinates in Center of Mass frame

    # Move Stars to CoM (initially in CoM)
    #### The stars are already in CoM coordinates ####

    # Move planetesimals to CoM (initially w/ respect to star zero)
    planetesimals.velocity += stars[center].velocity
    planetesimals.position += stars[center].position
    
    write_set_to_file(bodies, snap_dir + "/" + file_out, "hdf5")
    
    time += dt

  f.close() # stdout
     
  gravity.stop()
  stars_gravity.stop()
  planetesimals_gravity.stop()
  pl_gravity.stop()
  
  # retrieval?
  return stars, planetesimals

def new_option_parser():
  result = OptionParser()
  result.add_option("--sim_number", 
                    dest="sim_number", type="int", default = 9999,
                    help="simulation id number [%default]")
  result.add_option("--num_format", 
                    dest="num_format", type="int", default = 6,
                    help="number of digits in simulation id [%default]")
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 32,
                    help="number of steps [%default]")
  result.add_option("--fout", 
                    dest="fout", default="outfile",
                    help="output file base name[%default]")
  result.add_option("--snap_dir", 
                    dest="snap_dir", default="sim",
                    help="snapshot directory base name [%default]")
  result.add_option("--m0", unit=units.MSun,
                    dest="m0", type="float", default = 1.0|units.MSun,
                    help="mass of the disk-central star (with disk) in MSun [%default]")
  result.add_option("--m1", unit=units.MSun,
                    dest="m1", type="float", default = 1.0|units.MSun,
                    help="mass of the passing star (with no disk) in MSun [%default]")
  result.add_option("--peri", unit=units.AU,
                    dest="peri", type="float", default = 100|units.AU,
                    help="pericenter of the orbit in AU  [%default]")
  result.add_option("--ecc",
                    dest="ecc", type="float", default = 1.0,
                    help="eccentricity of the orbit  [%default]")
  result.add_option("--incl",
                    dest="incl", type="float", default = 0.0,
                    help="inclination of the orbit in deg [%default]")
  result.add_option("--omega",
                    dest="omega", type="float", default = 90.0,
                    help="argument of periapsis [%default]") # was 90.0 before
  ###### EXPERIMENT #####
  result.add_option("--rel_force",
                    dest="rel_force", type="float", default = 0.1,
                    help="initial relative force from m1 and m0 at the r_disk_out [%default]")
  ###### EXPERIMENT #####
  result.add_option("--n_disk", 
                    dest="n_disk", type="int", default=1000,
                    help="number of disk particles [%default]")
  result.add_option("--r_in", 
                    unit=units.AU, dest="r_in", type="float", default=40|units.AU,
                    help="inner radius of the disk in AU [%default]")
  result.add_option("--r_out", 
                    unit=units.AU, dest="r_out", type="float", default=100|units.AU,
                    help="outer radius of the disk in AU [%default]")
  result.add_option("--power", 
                    dest="power", type="float", default=1.5, 
                    help="negative index of the power-law surface number density of the disk particles [%default]") ### Careful with this one!! (and it is negative) ###
  result.add_option("--eta",
                    dest="eta", type="float", default=0.001,
                    help="Huayno eta parameter (~timestep) [%default]") # Huayno timestep here
  result.add_option("--fredir",
                    dest="fredir", type="string", default=None,
                    help="redirection file [%default]")
  ###### EXPERIMENT #####
  result.add_option("--br_dt",
                    dest="br_dt", type="float", default=0.01,
                    help="bridge timestep -- fraction of inner disk period [%default]") # Bridge timestep here
  ###### EXPERIMENT #####
  result.add_option("--num_T",
                    dest="num_T", type="float", default=4.0,
                    help="total simulation time (in 'pericenter arc' units) [%default]")
  result.add_option("--nb_end",
                    dest="nb_end", type="float", default=1.5,
                    help="end n-body time (in 'pericenter arc' units) [%default]")
  result.add_option("--seed",
                    dest="seed", type="int", default=42,
                    help="random seed for disk generator [%default]")

  return result
  
if __name__ in ('__main__', '__plot__'):
  """
  to reproduce Fig.1. in Lestrade+2011
  """
  
  # Option Parser Dictionary (This isn't a dictionary?!)
  o, arguments  = new_option_parser().parse_args()
  if o.num_format == 3:
     o.sim_number = "%03d" % o.sim_number
  else:
     o.sim_number = "%06d" % o.sim_number

  snap_base = o.snap_dir

  o.snap_dir = o.snap_dir + o.sim_number
  o.fout = o.fout + o.sim_number + ".hdf5"
  
  stars, time_peri = get_orbit_ini(o.m0, o.m1, o.peri, o.ecc, o.incl*pi_180, o.omega*pi_180, 
                     o.rel_force, o.r_out)
  
  planetesimals = get_planetesimals_disk(o.n_disk, o.r_in, o.r_out, o.m0, alpha = o.power, seed = o.seed)
  
  #print stars
  #print planetesimals

  ps.mkdir(o.snap_dir, safety = True)

  # Write Info File
  f = open(o.snap_dir + "/info.txt", 'a')
  f.write("Simulation ID Number: " + str(o.sim_number) + "\n")
  f.write("Disk Random Seed: " + str(o.seed) + "\n")
  f.write("Number of Pericenter Arcs: " + str(o.num_T) + "\n")
  f.write("End of N-Body Time: " + str(o.nb_end) + "\n")
  f.write("Mass (Disk Star): " + str(o.m0) + "\n")
  f.write("Mass (Other Star): " + str(o.m1) + "\n")
  f.write("Initial Relative Force: " + str(o.rel_force) + "\n")
  f.write("Bridge Timestep: " + str(o.br_dt) + "\n")
  f.write("Huayno Timestep (eta): " + str(o.eta) + "\n")
  f.write("Inner Disk Radius: " + str(o.r_in) + "\n")
  f.write("Outer Disk Radius: " + str(o.r_out) + "\n")
  f.write("Number of Disk Particles: " + str(o.n_disk) + "\n")
  f.write("Disk Particle Power Law: " + str(o.power) + "\n")
  f.write("Orbit Pericenter: " + str(o.peri) + "\n")
  f.write("Orbit Eccentricity: " + str(o.ecc) + "\n")
  f.write("Orbit Inclination: " + str(o.incl) + "\n")
  f.write("Orbit Argument of Periapsis: " + str(o.omega) + "\n")
  f.write("Number of Steps: " + str(o.n_steps) + "\n")
  f.close()

  # Write 'o' dictionary into file (using Pickle)
  pickle_f = open(o.snap_dir + "/info.p", "wb")
  pickle.dump(o, pickle_f)
  pickle_f.close()

  r_step = o.r_in ### set timestep based on period of innermost disk particles ###
  t_end = o.num_T * abs(time_peri) # Current Target Period: 4.0 T
  #t_end = 1300.0 | units.yr
  
  integrate_disk_flyby(stars, planetesimals, t_end, time_peri, o.nb_end, o.n_steps,
                       r_step, o.snap_dir, o.fout, o.fredir, o.br_dt, o.eta)

  factory = PF(int(o.sim_number), snapshot_base = snap_base, count = o.num_format)

  # Process snapshots ---> Make Plots
  factory.plot_all() # <<<--- 
  
