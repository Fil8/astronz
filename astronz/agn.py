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

class AGN:
      
  def __init__(self,args):



    self.c=C
    self.G=G
    self.mp=MP
    self.sigmaT=SIGMAT
    self.msun=MSUN
    self.myr=(3600*365*24)/MSUN
    self.yr=(3600*365*24)
    self.ledd=1.3e38
    self.rad=7.3e36
    self.oiii=3500
    self.effb=1
    self.eff=0.1
    self.p0=1.0e44   #erg/s
    self.l0=7e29     #erg/sHz
    self.ev= 1.602177e-12 #from ev to erg
    self.pc=PC
    
  def eddington(self,mbh):
    
    #set mbh in g:
    mbhg=mbh*self.msun
    a=4*math.pi*self.G*self.mp*self.c
    ledd=a*mbhg/self.sigmaT
    medd=ledd/(self.eff*(self.c)**2)*self.myr

    ledd1=1.3e38*mbh
    
    return medd, ledd, ledd1
    
    
  def mechanical(self,z,f,mbh):
    
    lum=c.luminosity(z,f)*1e-7
    lmec=self.rad*pow(lum/1e24,0.71) *1e7       #relation in W result in erg/s
    #a=a.eddington(mbh)[1]
    
    lmecedd=lmec/a.eddington(mbh)[1]        #in erg/s
    
    return lmec,lmecedd
  
  def radiative(self, z,foiii,mbh):
    
    dl=c.lum_dist(z)
    lum=(foiii*4*math.pi*(dl)**2)


    lrad=self.oiii*lum
    
    lradedd=lrad/a.eddington(mbh)[1]        #in erg/s
    
    return lrad, lradedd
    
  def lamb(self,z,f,foiii,mbh):
    
    lr=self.radiative(z,foiii,mbh)[1]
    lm=self.mechanical(z,f,mbh)[1]
    ledd=a.eddington(mbh)[1]

  
    lamb=(lm+lr)
  
  
  
    loglam=math.log10(lamb)
    
    return lamb, loglam
    

  def bondi (self,z,f,mbh):
    
    lr=c.luminosity(z,f)
    
    pjet=lr/self.l0
    pjet=(pjet**0.71)
    print lr

    pjet=1e44*pjet

    print pjet
    #pjet=1.8e44/1e43


    # from Balmaverde et al. 2008 Pjet => Pbondi = Mbondi*c^2
                                # Pjet= 0.012 Mbondi * c^2
    
    pbondi=(math.log10(pjet/1e43)+1.91)/1.1   #the formula is in 1e43 units
    pbondi=pow(10,pbondi)
    pbondi=pbondi*1e43

    mbondi=pbondi/((self.c**2)*1)*self.myr

    mbondijet=pjet/((self.c**2)*0.012)*self.myr
   
    # from Allen et al. 2006 Pjet => Pbondi = 0.1 * Mbondi *c^2
    pallen=(math.log10(pjet/1e43)*0.77)+0.65
    pallen=pow(10,pallen)
    pallen=pallen*1e43
    
    mbondia=pallen/((self.c**2)*0.1)*self.myr

    
    #bondi and mass accretion rate from BH mass (Willett 2010)
  
    #determine c_s
    kbt=0.7e3*self.ev
    cs=math.sqrt((5/3*kbt)/(0.62*self.mp))
    
    #determine the density of the gas
    ne=1
    rho=ne*1.13*self.mp

  
    #pbondi=4*pi lambda G^2 M_Bh^2 rho eff c^2 cs^-3
    
    #pbondi = 2* G M_Bh cs^-2
    
    ra=2*self.G*mbh*self.msun/(cs**2)
  
    lambdaa=0.25
    eff=0.1

    bondi=lambdaa*math.pi*cs*ra**2*rho*self.myr
    
    pbonditop=0.017*bondi/self.myr*self.c**2

    #bondiratio=bondi/self.myr*(self.c**2)*self.eff
    #print bondiratio

    #general relativity
    #you confront the mass accretion rate with the eddington luminosity
    
    ledd=a.eddington(mbh)[1]
    ledduni=ledd/self.c**2*self.myr
    
    accr=bondi/ledduni
    print accr
    
    ##about PKS1718-649
    #infall
    dist=100*self.pc*1e-5
    vcloud=70
    mcloud=1e4
    
    #pjet allen
    dl=c.lum_dist(z)

    minorax=0.4*1e-3*dl  #cm
    majorax=1.*1e-3*dl    #cm
    
    gamma=5./3./(5./3.-1.)
    T=1. #kev
    
    press=kbt*T*rho/(0.62*self.mp)
    print press
    
    V=math.pi*minorax*majorax**2
    print V

    tage=1e2*self.yr
    print tage
    
    pjetta=gamma*press*V/(tage*self.yr)

    t=dist/vcloud
    accrcloud=mcloud*self.msun/t*self.myr
    
    # from Allen et al. 2006 Pjet => Pbondi = 0.1 * Mbondi *c^2
    pallen2=(math.log10(pjetta/1e43)*0.77)+0.65
    pallen2=pow(10,pallen2)
    pallen2=pallen2*1e43
    
    mbondia2=pallen2/((self.c**2)*1)*self.myr

    
    

    return pjet,pbondi,mbondi,mbondijet,mbondia,bondi,accr,accrcloud,pjetta,mbondia2,pbonditop
    
    
  def Xaccretion (self,z,f,mbh):
    
    ledd=a.eddington(mbh)[1]
    lx=c.luminosity(z,f)
    lx=1e41
    lxratio=math.log10(lx/ledd)

    return lxratio
