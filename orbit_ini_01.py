"""
binary on orbit with given parameters
"""

import numpy 
import os

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.community.huayno.interface import Huayno
from amuse.community.kepler.interface import Kepler as Kepler_twobody
from amuse.io import write_set_to_file

pi_180 = numpy.pi/180.0

def orbital_parameters(rel_position, rel_velocity, total_mass):
    separation = rel_position.lengths()
    speed_squared = rel_velocity.lengths_squared()

    semimajor_axis = (constants.G * total_mass * separation / (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)
    eccentricity = numpy.sqrt( 1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=0) / (constants.G * total_mass * semimajor_axis))
    period = (2 * numpy.pi * semimajor_axis**1.5 / (numpy.sqrt(constants.G * total_mass))).as_quantity_in(units.yr)
    
    return semimajor_axis, eccentricity, period
  
def get_orbit_ini(m0, m1, peri, ecc, incl, omega, 
                  rel_force=0.01, r_disk=50|units.AU):
  converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  # semi-major axis
  if ecc!=1.0:
    semi = peri/(1.0-ecc)
  else:
    semi = 1.0e10 | units.AU
  
  # relative position and velocity vectors at the pericenter using kepler
  kepler = Kepler_twobody(converter)
  kepler.initialize_code()
  kepler.initialize_from_elements(mass=(m0+m1), semi=semi, ecc=ecc, periastron=peri) # at pericenter
  
  # moving particle backwards to radius r where: F_m1(r) = rel_force*F_m0(r_disk)
  #r_disk = peri
  r_ini = r_disk*(1.0 + numpy.sqrt(m1/m0)/rel_force)
  kepler.return_to_radius(radius=r_ini)
  
  rl = kepler.get_separation_vector()
  r = [rl[0].value_in(units.AU), rl[1].value_in(units.AU), rl[2].value_in(units.AU)] | units.AU
  vl = kepler.get_velocity_vector()
  v = [vl[0].value_in(units.kms), vl[1].value_in(units.kms), vl[2].value_in(units.kms)] | units.kms
  period_kepler = kepler.get_period()
  time_peri = kepler.get_time()
  
  kepler.stop()
  
  # rotation of the orbital plane by inclination and argument of periapsis
  a1 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a2 = ([numpy.cos(omega), -numpy.sin(omega), 0.0], [numpy.sin(omega), numpy.cos(omega), 0.0], [0.0, 0.0, 1.0])
  rot = numpy.dot(a1,a2)
  r_au = numpy.reshape(r.value_in(units.AU), 3, 1)
  v_kms = numpy.reshape(v.value_in(units.kms), 3, 1)
  r_rot = numpy.dot(rot, r_au) | units.AU
  v_rot = numpy.dot(rot, v_kms) | units.kms
  
  bodies = Particles(2)
  bodies[0].mass = m0
  bodies[0].radius = 1.0|units.RSun
  bodies[0].position = (0,0,0) | units.AU
  bodies[0].velocity = (0,0,0) | units.kms
  bodies[1].mass = m1
  bodies[1].radius = 1.0|units.RSun
  bodies[1].x = r_rot[0]
  bodies[1].y = r_rot[1]
  bodies[1].z = r_rot[2]
  bodies[1].vx = v_rot[0]
  bodies[1].vy = v_rot[1]
  bodies[1].vz = v_rot[2]
  
  bodies.age = 0.0 | units.yr
  bodies.move_to_center()
  
  print "\t r_rel_ini  = ", r_rot.in_(units.AU)
  print "\t v_rel_ini  = ", v_rot.in_(units.kms)
  print "\t time since peri = ", time_peri.in_(units.yr)
  a_orbit, e_orbit, p_orbit = orbital_parameters(r_rot, v_rot, (m0+m1))
  print "\t a = ", a_orbit.in_(units.AU), "\t e = ", e_orbit, "\t period = ", p_orbit.in_(units.yr)
  
  return bodies, time_peri
  
def orbit_evolve(bodies, time_peri, eta, n_steps=100):
  converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  gravity = Huayno(converter,channel_type="sockets")
  gravity.particles.add_particles(bodies)
  gravity.commit_particles()
  gravity.parameters.timestep_parameter = eta
  
  channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
  
  Etot_init = gravity.kinetic_energy + gravity.potential_energy
  Etot = Etot_init
  
  file_snap = "orbit_ini.hdf5"
  
  t_end = 2.0*abs(time_peri.value_in(units.yr)) | units.yr
  dt = t_end / float(n_steps)
  time = 0.0 | units.yr
  
  while time <= t_end:
    
    gravity.evolve_model(time)
    channel_from_gravity_to_framework.copy()
    bodies.age = time
    
    # checking energy conservation
    Ekin = gravity.kinetic_energy
    Epot = gravity.potential_energy
    Etot = Ekin + Epot
    dE = Etot_init-Etot
    
    # output centerd on the bodies[0] star
    bodies.position -= bodies[0].position
    bodies.velocity -= bodies[0].velocity
    write_set_to_file(bodies, file_snap, "hdf5")
    
    rel_r = (bodies[1].position - bodies[0].position).lengths()
    print " \t\t", time, "\t", dE/Etot_init, rel_r.in_(units.AU)
    time += dt
  
  gravity.stop()
  bodies.position -= bodies[0].position
  bodies.velocity -= bodies[0].velocity
  print bodies
  
  return

def new_option_parser():
  result = OptionParser()
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 10,
                    help="number of steps [%default]")
  result.add_option("--m0", unit=units.MSun,
                    dest="m0", type="float", default = 1.0|units.MSun,
                    help="mass of the disk-central star in MSun [%default]")
  result.add_option("--m1", unit=units.MSun,
                    dest="m1", type="float", default = 1.0|units.MSun,
                    help="mass of the passing star in MSun [%default]")
  result.add_option("--peri", unit=units.AU,
                    dest="peri", type="float", default = 200|units.AU,
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
  result.add_option("--eta",
                    dest="eta", type="float", default=0.001,
                    help="Huayno eta parameter (~timestep) [%default]")
  return result


if __name__ in ('__main__', '__plot__'):
  o, arguments  = new_option_parser().parse_args()
  
  bodies, time_peri = get_orbit_ini(o.m0, o.m1, o.peri, o.ecc, o.incl*pi_180, o.omega*pi_180)
  orbit_evolve(bodies, time_peri, o.eta, o.n_steps)

