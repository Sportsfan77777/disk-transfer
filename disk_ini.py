"""
Intializes particles of planetesimal disk using built-in ProtoPlanetaryDisk function and input parameters
-- n_disk particles
-- r_in<r<r_out
-- circular orbits around mass m_star
-- in the plane z=0
-- possibily to get power-law surface number density of the disk particles (using alpha)
"""

import numpy 
import os

from amuse.units import units
from amuse.units import constants
from amuse.datamodel import Particles
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.ext.protodisk import ProtoPlanetaryDisk

def get_planetesimals_disk(n_disk, r_in=20.0|units.AU, r_out=50.0|units.AU, m_star=1.0|units.MSun, 
                           alpha=None, m_disk=1.0e-15|units.MSun, seed = 42, disk_num = 1):
    """
    returns particles of planetesimal disk with given parameters
    """
    numpy.random.seed(seed) # Mess with random seed

    for i in xrange(int(disk_num) + 4):
        planetesimals = Particles(n_disk)
    print "Seed:", seed
    print planetesimals.key[:10]
    planetesimals.mass = 0.0|units.MJupiter
    planetesimals.radius = 100.0|units.km
    planetesimals.collection_attributes.timestamp = 0.0 | units.yr
    
    if alpha is not None:
      converter = nbody_system.nbody_to_si(m_disk, r_in)
      power_disk = ProtoPlanetaryDisk(n_disk, convert_nbody=converter, densitypower=alpha, 
                                      Rmin=1.0, Rmax=1.0*r_out/r_in, q_out=0.0, discfraction=1.0).result
      x = power_disk.x
      y = power_disk.y
      z = power_disk.z # <--- Mystery error?

      print "X"
      print x.value_in(units.AU)
      print "Y"
      print y.value_in(units.AU)
      print "Z"
      print z.value_in(units.AU)
      #z = 0

      print "MASS"
      print power_disk.mass

      #power_disk.mass = 0.0 * power_disk.mass ###### THIS WORKS!!!! (if you want to switch to this later)

      print power_disk.mass
      
      a = (x**2 + y**2)**0.5
      print "SM-AXIS"
      print a.value_in(units.AU)
      
      phi = numpy.arctan2(y.value_in(units.AU), x.value_in(units.AU))
      vc = (constants.G*m_star/a)**0.5
      vx = - vc * numpy.sin(phi)
      vy = vc * numpy.cos(phi)
      vz = 0.0 * vc
      # vz = - vc * numpy.sin(phi) # ???????????????????????????????????????????????????????????????? #

      print "VX"
      print vx.value_in(units.km / units.s)
      print "VY"
      print vy.value_in(units.km / units.s)
      print "VZ"
      print vz.value_in(units.km / units.s)

      print "PLANAR VELOCITY VECTOR"
      print ((vx**2 + vy**2)**(0.5)).value_in(units.km / units.s)

      #vx = power_disk.vx
      #vy = power_disk.vy
      #vz = power_disk.vz

      #print "POWER DISK VX"
      #print vx.value_in(units.km / units.s)
      #print "POWER DISK VY"
      #print vy.value_in(units.km / units.s)
      #print "POWER DISK VZ"
      #print vz.value_in(units.km / units.s)
      
    else:
      a = r_in + (r_out-r_in)*numpy.random.rand(n_disk)
      phi_rand = 2.0 * numpy.pi * numpy.random.rand(n_disk)
   
      x = a * numpy.cos(phi_rand)
      y = a * numpy.sin(phi_rand)
      z = 0.0 * a
   
      vc = (constants.G*m_star/a)**0.5
      vx = - vc * numpy.sin(phi_rand)
      vy = vc * numpy.cos(phi_rand)
      vz = 0.0 * vc
      
    planetesimals.x = x
    planetesimals.y = y
    planetesimals.z = z
    
    planetesimals.vx = vx
    planetesimals.vy = vy
    planetesimals.vz = vz
    
    return planetesimals
    
    
