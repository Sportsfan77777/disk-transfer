"""
close encounters in star cluster
  some of used equtions adopted from Lestrade+11
"""

import numpy
from scipy import integrate
from amuse.units import units, constants

def number_density_t(t, tcl, pop_fraction, n0, n1=0.1|units.parsec**(-3)):
  """
  star number density [pc^-3], eq. (3) in Lestrade+11
  linerly decreasing with time -- expanding cluster
  t -- time
  tcl -- cluster lifetime
  n0 -- initial density
  n1 -- final density (after dissolution, e.i. the field stars density)
  pop_fraction -- fraction of given popula
  tion (in given mass range)
  """
  den = n0*pop_fraction - (n0-n1)*pop_fraction*(t/tcl)
  return den

def kroupa_imf(m):
  """
  IMF, Kroupa 2002
  m [MSun]
  """
  if ((0.08<=m<0.1)):
    pm = m**-0.3
  elif ((0.1<=m<0.5)):
    pm = m**-1.3
  elif ((0.5<=m<120.)):
    pm = m**-2.3
  else:
    pm = 0.
  
  return pm

def pop_fraction(m_min, m_max):
  """
  number fraction of given mass bin m_min--m_max [MSun]
  using Kroupa IMF
  """
  # auxiliary constant 
  #  total probability of Kroupa IMF
  imf_tot = integrate.quad(kroupa_imf,0.08,0.1)[0] + \
    integrate.quad(kroupa_imf,0.1,0.5)[0] + \
    integrate.quad(kroupa_imf,0.5,120.)[0]
  
  pop_frac = integrate.quad(kroupa_imf,m_min,m_max)[0]/imf_tot
  
  return pop_frac

def encounter_rate(dens, sigma, impact, mc, mp):
  """
  encounter rate, eq. (1)
  """
  a_const, b_const = t_per_n_const(sigma, mc, mp)
  rate = a_const * dens * impact**2 + \
    b_const * dens * impact
  return rate
  
def t_per_n_const(sigma, mc, mp):
  # auxiliary constants
  a_const = 4.*numpy.sqrt(numpy.pi)*sigma
  b_const = 2.*numpy.sqrt(numpy.pi)*constants.G*(mc+mp) / sigma
  
  return a_const, b_const
  
def encounter_d_for_n(n_enc, sigma, mc, mp, pop_frac, n0, n1, tcl):
  """
  solving eq. (4) for d
  """
  
  # auxiliary constants
  a_const, b_const = t_per_n_const(sigma, mc, mp)
  c_const = 0.5*pop_frac*tcl*(n0+n1)
  
  # quadratic equation coeff.
  a = a_const*c_const
  b = b_const*c_const
  c = -1.*n_enc

  d = b*b - 4.*a*c
  x1 = (-b + numpy.sqrt(d)) / (2.*a)
  x2 = (-b - numpy.sqrt(d)) / (2.*a)
  
  return x1, x2

def wrong_d_for_n(n_enc, sigma, mc, mp, pop_frac, n0, n1, tcl):
  """
  solving eq. (4) for d
  usinf Lestrade's coefficients in eq. (2) -- these are WRONG
  """
  
  # auxiliary constants
  a_const = (1.9e-8) *sigma/(1.0|units.kms) / (100.0|units.AU)**2
  b_const = (8.8e-9) *(mc+mp)/(1.0|units.MSun) / (sigma/(1.0|units.kms)) / (100.0|units.AU)
  c_const = 0.5*pop_frac*tcl/(1.0|units.yr)*(n0+n1)/(1000|units.parsec**-3)
  
  # quadratic equation coeff.
  a = a_const*c_const
  b = b_const*c_const
  c = -1.*n_enc

  d = b*b - 4.*a*c
  x1 = (-b + numpy.sqrt(d)) / (2.*a)
  x2 = (-b - numpy.sqrt(d)) / (2.*a)
  
  return x1, x2

def encounter_n_per_t_d(d, sigma, mc, mp, pop_frac, tcl, n0, n1=0.1|units.parsec**-3):
  """
  int_0^tcl (1/t_enc) dt
 
  """
  # auxiliary constants
  a_const, b_const = t_per_n_const(sigma, mc, mp)
  c_const = 0.5*pop_frac*tcl*(n0+n1)
  n_enc = c_const*(a_const*d**2 + b_const*d)
  return n_enc

  
def enc_rate(d, sigma, mc, mp, den):
  """
  Binney and Tremaine 1987, eq. (8-122), 
    where 2m=(mc+mp) since we do not have equal mass case
  """
  rate = 4.0*numpy.sqrt(numpy.pi)*den*sigma*d**2 +\
    2.0*numpy.sqrt(numpy.pi)*constants.G*(mp+mc)*den*d / sigma
  return rate
  
def get_encountarers_for_bins(dens=3000.|units.parsec**-3, 
                              sigma=5.|units.kms, 
                              mc=1.|units.MSun,
                              tcl=100.|units.Myr):
  """
  get encounter rates -- # of encounters per time
  for table of different masses and impact parameters
  cluster given by
    dens -- number density of the cluster [m^-3]
    sigma -- velocity dispersion [kms]
    mc -- mass of the 'central' star [kg]
    tcl -- lifetime [s]
  """
  
  # bins in Mp [MSun] / b [AU] Michael's table
  mcbins = numpy.array([0.1, 0.2, 0.25, 0.333, 0.5, 0.667, 1., 1.5, 2., 2.5, 3.])
  dcbins = numpy.linspace(40.,300.,14)
  
  #for full mass and impact spectrum
  #mbins = numpy.linspace(0.001,121.,20)
  #dbins = numpy.linspace(0.,200000.,10) * (1.|units.AU)
  
  # mass bins
  mbins = numpy.zeros(numpy.size(mcbins)+1)
  mbins[0] = 0.  # lower mass limit
  mbins[-1] = 10. # upper mass limit
  for i,mi in enumerate(mcbins[:-1]):
    mbins[i+1] = 0.5*(mcbins[i]+mcbins[i+1])

  # the population fractions for mass bins
  pop_fraction_mbins = numpy.zeros_like(mcbins)
  for i,mi in enumerate(mcbins[:]):
    pop_fraction_mbins[i] = pop_fraction(mbins[i],mbins[i+1])
  #print pop_fraction_mbins
  
  # impact factor bins
  dbins = numpy.zeros(numpy.size(dcbins)+1)
  dbins[0] = 30.   # lower impact limit
  dbins[-1] = 310. # upper impact limit
  for i,mi in enumerate(dcbins[:-1]):
    dbins[i+1] = 0.5*(dcbins[i]+dcbins[i+1])

  dbins[0] = 0.0
  
  # get rates per 1myr assuming constant denstity dens
  enc_myr = numpy.zeros([numpy.size(mcbins), numpy.size(dcbins)])
  for i,di in enumerate(dcbins[:]):
    for j,mj in enumerate(mcbins[:]):
      enc_myr_d0 = encounter_rate(dens, sigma, dbins[i]*(1.|units.AU), mc, mj*(1.|units.MSun)) \
        * pop_fraction_mbins[j] * (1.|units.Myr)
      enc_myr_d1 = encounter_rate(dens, sigma, dbins[i+1]*(1.|units.AU), mc, mj*(1.|units.MSun)) \
        * pop_fraction_mbins[j] * (1.|units.Myr)
      enc_myr[j,i] = enc_myr_d1 - enc_myr_d0
      
  # get rates for cluster with linearly decreasing density over t_cl
  enc_tcl = numpy.zeros([numpy.size(mcbins), numpy.size(dcbins)])
  for i,di in enumerate(dcbins[:]):
    for j,mj in enumerate(mcbins[:]):
      enc_myr_d0 = encounter_n_per_t_d(dbins[i]*(1.|units.AU), sigma, 
                                       mc, mj*(1.|units.MSun), pop_fraction_mbins[j], tcl, dens)
      enc_myr_d1 = encounter_n_per_t_d(dbins[i+1]*(1.|units.AU), sigma, 
                                       mc, mj*(1.|units.MSun), pop_fraction_mbins[j], tcl, dens)
      enc_tcl[j,i] = enc_myr_d1 - enc_myr_d0
  
  #print enc_myr
  #print enc_tcl
  
  return enc_tcl, enc_myr
  
def sun_example():
  
  # read Michael's table
  counts = numpy.loadtxt('count_table.dat') / 100.
  
  # number of encounters for the table
  # and Sun's birth cluster from Lestrade+
  ref = Lestrade()
  enc_tcl, enc_myr = get_encountarers_for_bins(ref.n0, ref.sigma, 1.|units.MSun, ref.tcl)
  # transferred rate * encounter rate
  tot = counts * enc_tcl
  print ' ** Lestrade+2011', sum(sum(tot))
  
  # Sun's birth cluster from Adams2010
  ref = Adams2010()
  enc_tcl, enc_myr = get_encountarers_for_bins(ref.n0, ref.sigma, 1.|units.MSun, ref.tcl)

  print enc_myr.shape
  print enc_myr

  peri = range(40, 301, 20) #### <<<---- USE THIS

  mass_range = [1.0/10, 1.0/5, 1.0/4, 1.0/3, 1.0/2, 2.0/3, 1.0, 3.0/2, 2.0, 5.0/2, 3.0]  #### Note 1.0 is before 0.67
  mass = [round(mr,2) for mr in mass_range] #### <<<---- USE THIS

  print "For each impact parameter"
  total = 0
  total -= 0.00633953383867
  #total = 0
  for i in range(14):
      print peri[i]
      count = sum(enc_myr[:,i])
      print count
      total += count
      print 1.0 / count
      total += count
      print total
      print 1.0 / total
      print

  """
  print "For each disk star mass"
  for i in range(11):
      print mass[i]
      print sum(enc_myr[i,:])
      print 1.0 / sum(enc_myr[i,:])
      print
  """

  tot = counts * enc_tcl
  print ' ** Adams2012', sum(sum(tot))
  
  return
  

class Lestrade:
  """ 
  Sun's birth cluster from Lestrade+11
  """
  def __init__(self):
    self.n0 = 3000. |units.parsec**-3
    self.n1 = 0.1 |units.parsec**-3
    self.sigma = 5. | units.kms
    self.tcl = 100. | units.Myr
    self.mc = 1. | units.MSun
    self.mass_bins_l11 = [0.1,0.21], [0.21,0.47], [0.47,0.8], [0.8,1.7], [1.7,3.2], [3.2,6.5]
    self.pop_frac_l11 = 0.43, 0.315, 0.124, 0.085, 0.027, 0.013
    
class Adams2010:
  """ 
  Sun's birth cluster from Adams2010
  see Sec. 4.1
  """
  def __init__(self):
    self.n0 = 300. |units.parsec**-3
    self.n1 = 0.1 |units.parsec**-3
    self.sigma = 3. | units.kms
    self.tcl = 10. | units.Myr
    
if __name__ in ('__main__', '__plot__'):
  sun_example()


  
