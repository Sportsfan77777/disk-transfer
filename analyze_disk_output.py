"""
Analyze results of stellar flyby trying to steal particles from another star's debris disk
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


def parse_output(f1):
    """
    f1 = output file (in .hdf5 format)
    """
    bodies = read_set_from_file(f1, 'hdf5')


def new_option_parser():
  result = OptionParser()
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 13,
                    help="number of steps [%default]")
  result.add_option("--fout", 
                    dest="fout", default="disk_flyby.hdf5",
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
                    dest="peri", type="float", default = 200|units.AU,
                    help="pericenter of the orbit in AU  [%default]")
  result.add_option("--ecc",
                    dest="ecc", type="float", default = 1.0,
                    help="eccentricity of the orbit  [%default]")
  result.add_option("--incl",
                    dest="incl", type="float", default = 0.0,
                    help="inclination of the orbit in deg [%default]")
  result.add_option("--n_disk", 
                    dest="n_disk", type="int", default=1000,
                    help="number of disk particles [%default]")
  result.add_option("--r_in", 
                    unit=units.AU, dest="r_in", type="float", default=40|units.AU,
                    help="inner radius of the disk in AU [%default]")
  result.add_option("--r_out", 
                    unit=units.AU, dest="r_out", type="float", default=100|units.AU,
                    help="outer radius of the disk in AU [%default]")
  result.add_option("--fredir",
                    dest="fredir", type="string", default=None,
                    help="redirection file [%default]")
  result.add_option("--center",
                    dest="center", type="int", default=0,
                    help="central star in plot [%default]")
  return result


if __name__ in ('__main__', '__plot__'):
  """
  determine # of particles in stable orbits around 'passing star'
  and properties regarding their new orbital parameters

  re-plot so that "stable particles" are green and "unstable particles" are orange (stars are already red)?
  stars = "Red"
  bound to center = "Blue"
  bound to passing = "Green"
  bound to both? = "Purple" <--- Is this possible? (probably?)
  unbound = "Orange"
  """
  
  o, arguments  = new_option_parser().parse_args()
  
  stars, time_peri = get_orbit_ini(o.m0, o.m1, o.peri, o.ecc, o.incl*pi_180, o.omega*pi_180, 
                     o.rel_force, o.r_out)
  
  planetesimals = get_planetesimals_disk(o.n_disk, o.r_in, o.r_out, o.m0)
  
  #print stars
  #print planetesimals
  
  r_step = o.r_in
  t_end = 25.0*abs(time_peri)
  #t_end = 1300.0 | units.yr
  
  integrate_disk_flyby(stars, planetesimals, t_end, o.n_steps,
                       r_step, o.fout, o.fredir, o.br_dt, o.eta, o.center)
                        
  plot_all_snaps(o.fout, o.snap_dir)
  
