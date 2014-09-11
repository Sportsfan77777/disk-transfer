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

from disk_ini import get_planetesimals_disk
from orbit_ini_01 import get_orbit_ini

pi_180 = numpy.pi/180.0

def rm_file(file_name):
  """
  delete the file 'file_name', if it exists
  """  
  if(os.path.isfile(file_name)):
    os.system('rm %s' % (file_name))

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()
        
def integrate_disk_flyby(stars, planetesimals, t_end, n_steps,
                         snap_dir, file_out, file_redir, huayno_eta):
  converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  planetesimals.position += stars[0].position
  planetesimals.velocity += stars[0].velocity
  
  bodies = ParticlesSuperset([stars, planetesimals])
  
  if file_redir is None:
    gravity = Huayno(converter, channel_type="sockets",
                       mode="openmp") # Add this (specify a number of cores)
  elif file_redir=="0":
    gravity = Huayno(converter, channel_type="sockets", redirection="none",
                       mode="openmp")
  else:
    gravity = Huayno(converter, channel_type="sockets", redirection="file", redirect_file=file_redir,
                       mode="openmp")
                       
  gravity.particles.add_particles(bodies)
  gravity.commit_particles()
  gravity.parameters.timestep_parameter = huayno_eta
  time_step = gravity.get_timestep_parameter()
  #gravity.set_inttype_parameter(12)
  #gravity.set_inttype_parameter(8)
  gravity.set_eps2_parameter(0.001*0.001) # Softening Parameter
  print ' ** timestep: ', gravity.get_timestep_parameter()
  print ' ** inttype: ', gravity.get_inttype_parameter()
  print ' ** eps2: ', gravity.get_eps2_parameter(), numpy.sqrt(gravity.get_eps2_parameter())
  
  t0 = time.time()
  evolve_disk_flyby(bodies, gravity, 
                    t_end, n_steps, converter, snap_dir, file_out)
  t1= time.time()
  dt = t1-t0
  print "Performace data: N =", len(bodies), "dt=", dt, "s =", dt/60.0, "min"
  
  return
  
def evolve_disk_flyby(bodies, gravity, 
                      t_end, n_steps, converter, snap_dir, file_out):
  bodies = ParticlesSuperset([stars, planetesimals])
  
  channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)
  
  rm_file(file_out)
    
  Etot_init = gravity.kinetic_energy + gravity.potential_energy
  Etot = Etot_init
  
  dt = t_end / float(n_steps)
  time = 0.0 | units.yr

  stdout = (file_out.split('.'))[0]
  stdout += '.txt'
  f = open(snap_dir + "/" + stdout, 'a')
  
  print " ** evolving: t_end = ", t_end, ", dt = ", dt.in_(units.yr)
  print " \t\t", "time", "\t\t", "dE"
  while time<=t_end:
    gravity.evolve_model(time)
    channel_from_gr_to_framework.copy()
    
    bodies.collection_attributes.timestamp = time

    Ekin = gravity.kinetic_energy 
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot

    nb_E = converter.to_nbody(Etot)
    #nb_J = converter.nbody_length ** 2 * units.nbody_mass * units.nbody_time ** -2 # not supposed to work

    # A formatted string would work better than tabs. (Tabs never work)
    line = " \t" + str(time.value_in(units.yr)) + "\t" + str(nb_E) + "\t" + str(dE/Etot_init)
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
  
  return

def new_option_parser():
  result = OptionParser()
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 10,
                    help="number of steps [%default]")
  result.add_option("--fout", 
                    dest="fout", default="disk_flyby_nbody.hdf5",
                    help="output file [%default]")
  result.add_option("--snap_dir", 
                    dest="snap_dir", default="disk_flyby_lestrade11",
                    help="output file [%default]")
  result.add_option("--m0", unit=units.MSun,
                    dest="m0", type="float", default = 1.0|units.MSun,
                    help="mass of the disk-central star in MSun [%default]")
  result.add_option("--m1", unit=units.MSun,
                    dest="m1", type="float", default = 1.0|units.MSun,
                    help="mass of the passing star in MSun [%default]")
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
                    help="argument of periapsis [%default]")
  result.add_option("--rel_force",
                    dest="rel_force", type="float", default = 0.01,
                    help="initial relative force from m1 and m0 at the r_disk_out [%default]")
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
                    dest="power", type="float", default=None, 
                    help="index of the power-law surface number density of the disk particles [%default]")
  result.add_option("--eta",
                    dest="eta", type="float", default=0.0001,
                    help="Huayno eta parameter (~timestep) [%default]")
  result.add_option("--fredir",
                    dest="fredir", type="string", default=None,
                    help="redirection file [%default]")
  result.add_option("--br_dt",
                    dest="br_dt", type="float", default=0.1,
                    help="bridge timestep -- fraction of inner disk period [%default]") # Bridge timestep here
  return result
  
if __name__ in ('__main__', '__plot__'):
  o, arguments  = new_option_parser().parse_args()
  
  stars, time_peri = get_orbit_ini(o.m0, o.m1, o.peri, o.ecc, o.incl*pi_180, o.omega*pi_180, 
                     o.rel_force, o.r_out)
  planetesimals = get_planetesimals_disk(o.n_disk, o.r_in, o.r_out, o.m0)
  
  #print stars
  #print planetesimals

  # Number of Periods
  num_T = 10.0

  try:
     os.mkdir(o.snap_dir)
  except:
     print "\t(" , o.snap_dir, "already exists)"

  # Write Info Files
  f = open(o.snap_dir + "/info.txt", 'a')
  f.write("Number of Periods: " + str(num_T) + "\n")
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
  
  r_step = o.r_in
  t_end = num_T * abs(time_peri.value_in(units.yr)) | units.yr
  
  integrate_disk_flyby(stars, planetesimals, t_end, o.n_steps, o.snap_dir, o.fout, o.fredir, o.eta)
  
