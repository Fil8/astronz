#!/usr/bin/env python

# ASTRONZ: astro calc utilities


import sys,string,os,math,numpy

RAD2DEG=180./math.pi
HI=1.42040575177e+09 #GHz
TSPIN=100 #K
MSUN=1.98855e33
MHI=1.6749E-24
CHI=2.36E5
PC=3.08567758E18 #cm

JANSKY=1e-23    #erg/scm2Hz

C=2.99792458E10    #cm/s
G=6.6742E-08       #cm3kg-1s-1      
MP=1.67492728E-24  #g
SIGMAT=6.66524E-25 #cm2


class Cosmo:
  
  def __init__(self, h0=69.6, omega_l=0.714, omega_m=0.286, ztime=False):


    #define constants
    # H0
    self.h0=h0

    # Cosmological constant
    self.cosmo_c=omega_l/(3.*self.h0**2.)

    # q0
    if omega_l == 0:
      self.q0=omega_m/2.
    else:
      self.q0=(3.*omega_m/2.) - 1.

    self.age_uni = self.compute_age(0)

    # z-time curve
    if ztime:
      self.z_time_compute()

    self.jy=JANSKY

  def compute_age(self, z):

    # Returns time at redshift z (from GALAXEV)

    a = lambda z: (math.sqrt(1. + ((2. * self.q0) * z)) / (1. - (2. * self.q0))) / (1. + z)
    c = lambda z: ((1. - (self.q0 * (1. - z))) / self.q0) / (1. + z)
    d = lambda x, omegainv: math.sqrt(omegainv) / (x*math.sqrt((omegainv-1.)+(1./(x**3.))))

    hh0 = self.h0 * 0.001022 # in (billion years)**(-1)

    if self.cosmo_c:

      omega0 = (2. * (self.q0 + 1.)) / 3.
      aa = 0
      bb = 1. / (1. + z)

      ok=0
      s0=1.e-10
      npts=0
      while not ok:
        npts=npts+1

        if npts==1:
          s=(bb-aa)*d(0.5*(aa+bb),1/omega0)
        else:
          it=3**(npts-2)
          tnm=it
          dd=(bb-aa)/(3.*tnm)
          ddel=dd+dd
          x=aa+0.5*dd
          sum=0.
          for j in range (1, it+1):
            sum=sum+d(x,1/omega0)
            x=x+ddel
            sum=sum+d(x,1/omega0)
            x=x+dd
          s=(s+(bb-aa)*sum/tnm)/3.

        epsr=math.fabs(s-s0)/s0
        if epsr < 1.0e-4:
          ok=True
        else:
          s0=s

      t=s

    elif self.q0==0:
      t = 1. / (1. + z)

    elif self.q0==0.5:
      t = (2. / 3.) / ((1. + z) ** 1.5)

    else:

      b = self.q0 / (math.fabs((2. * self.q0) - 1.) ** 1.5) 

      if self.q0<0.5:
        t = a(z) - (b * math.cosh(c(z))) 

      else:
        t = a(z) + (b * math.cos(c(z)))


    t = t / hh0

    return t


  def lum_dist(self, z):

    # Computes luminosity distance corresponding to a redshift z.
    # Uses Mattig formulae for qo both 0 and non 0
    # Ho in km/sec/Mpc
    # DL is in cm
    #
    # from GALAXEV

    e = lambda x, omegainv: 1. / math.sqrt(((x ** 3.) + omegainv) - 1.)

    if z <=0:
      # 10 pc
      return 1.e-5

    if self.q0 == 0:
      dl = ((3.e5 * z) * (1 + (z / 2.))) / self.h0
    elif self.q0 > 0:
      d1 = (self.q0 * z) + ((self.q0 - 1.) * (math.sqrt(1. + ((2. * self.q0) * z)) - 1.))
      d2 = ((self.h0 * self.q0) * self.q0) / 3.e5
      dl = d1 / d2
    elif self.q0 < 0:
      omega0 = (2. * (self.q0 + 1.)) / 3.
      aa = 1.
      bb = 1. + z
      ok=None
      s0=1.e-10
      npts=0
      while not ok:
        npts=npts+1

        if npts==1:
          s=(bb-aa)*e(0.5*(aa+bb),1/omega0)
        else:
          it=3**(npts-2)
          tnm=it
          dd=(bb-aa)/(3.*tnm)
          ddel=dd+dd
          x=aa+0.5*dd
          sum=0.
          for j in range (1, it+1):
            sum=sum+e(x,1/omega0)
            x=x+ddel
            sum=sum+e(x,1/omega0)
            x=x+dd
          s=(s+(bb-aa)*sum/tnm)/3.

        epsr=abs(s-s0)/s0
        if epsr < 1.e-4:
          ok=True
        else:
          s0=s
      dd1=s
      dd2 = (3.e5 * (1. + z)) / (self.h0 * math.sqrt(omega0))
      dl = dd1 * dd2

    dl=dl*3.085678e24 

    return dl

  def z_time_compute(self):

    # create the time/z curve
    self.z_time=[[],[]]
    for z in numpy.arange(0, 30, 0.01):

      self.z_time[0].append(self.compute_age(z))
      self.z_time[1].append(z)


  def lin2ang(self, r, z): # r in Mpc

    dl = self.lum_dist(z)/3.085678e24 # Mpc

    ang = RAD2DEG * 3600. * r * (1.+z)**2 / dl # arcsec

    return ang


  def ang2lin(self, ang, z): # r in arcsec

    dl = self.lum_dist(z)/3.085678e24 # Mpc
    r = ang * dl / (RAD2DEG * 3600 * (1+z)**2) # Mpc

    return r
    
  def luminosity(self,z, flux):
    
    f=flux*self.jy
    #print f
    
    dl=self.lum_dist(z)
    #print dl
    lum=f*4*math.pi*pow(dl,2)
    
    return lum
    
  def main(self):

    in1="\t... Here's what I can do: ...\n"
    in2='''\n\t - luminosity distance (DL)
\t - size from arcsec (LA)
\t - size from Mpc (LR)
\t - age at redshift (CA)
\t - luminosity (L)\n
    '''
    
    inp=str(raw_input(in1+in2))
    
    if (inp=='dl' or inp=='lumdist' or inp=='DL'):
        in1="\nz =\t "
        dl=self.lum_dist(float(raw_input(in1)))/3.085678e24
        print("D_L =\t %g\t Mpc"% (dl))
    elif(inp=='la' or inp=='scalearcsec' or inp=='LA'):
        in1="\nz= "
        in2="R(kpc)= "
        z=float(raw_input(in1))
        r=float(raw_input(in2))/1e3   #since the module wants the radius in Mpc
        print("R =\t %g\t arcsec"% (self.lin2ang(r,z)))
    elif(inp=='lr' or inp=='scaleMpc' or inp=='LR'):
        in1="\nz =\t "
        in2="R [arcsec] =\t "
        z=float(raw_input(in1))
        r=float(raw_input(in2))
        print("R = %g kpc"% (self.ang2lin(r,z)*1e3))
    elif(inp=='ca' or inp=='scalearcsec' or inp=='CA'):
        in1="\nz =\t "
        dl=self.compute_age(float(raw_input(in1)))
        print("Age =\t %g\t Gyears"% (dl))
    elif(inp=='l' or inp=='luminosity' or inp=='L'):
        in1="\nz= "
        in2="F(Jy)= "
        z=float(raw_input(in1))
        f=float(raw_input(in2))
        a=(self.luminosity(z,f))
        print("L =\t %g\t erg/s\n"% a)
        b=a*1e-7
        print("L =\t %g\t W"% b)
        c=a*HI
        print("L =\t %g\t W"% c)
    else:
        print ('\n\t ... you have not entered correct function ... \n')
        print ('\t************* --- COSMO : ERROR --- **************\n')
        sys.exit(0)

    print ('\n\t************* --- COSMO : OVER --- **************\n')

          