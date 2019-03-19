#!/usr/bin/env python

# ASTRONZ: astro calc utilities


import sys,string,os
import numpy as np

import cosmo


class AGN:
      
  def __init__(self):

    self.cc = cosmo.Cosmo()


    self.rad2deg = 180./np.pi
    self.C=2.99792458E10    #cm/s
    self.G=6.6742E-08       #cm3kg-1s-1
    self.mp=1.67492728E-24  #g
    self.sigmaT=6.66524E-25 #cm2
    self.msun=1.98855e33
    self.myr=(3600*365*24)/self.msun
    self.yr=(3600*365*24)
    self.ledd=1.3e38
    self.rad=7.3e36
    self.oiii=3500
    self.effb=1
    self.eff=0.1
    self.p0=1.0e44   #erg/s
    self.l0=7e29     #erg/sHz
    self.ev= 1.602177e-12 #from ev to erg
    self.pc=3.08567758E18 #cm
    
  def eddington(self,mbh):
    
    #set mbh in g:
    mbhg=mbh*self.msun
    a=4*np.pi*self.G*self.mp*self.C
    ledd=a*mbhg/self.sigmaT
    medd=ledd/(self.eff*(self.C)**2)*self.myr

    ledd1=1.3e38*mbh
    
    return medd, ledd, ledd1
    
    
  def mechanical(self,z,f,mbh):
    
    lum=self.cc.luminosity(z,f)*1e-7
    lmec=self.rad*pow(lum/1e24,0.71) *1e7       #relation in W result in erg/s
    #a=self.eddington(mbh)[1]
    
    lmecedd=lmec/self.eddington(mbh)[1]        #in erg/s
    
    return lmec,lmecedd
  
  def radiative(self, z,foiii,mbh):
    
    dl=self.cc.c.lum_dist(z)
    lum=(foiii*4*np.pi*(dl)**2)


    lrad=self.oiii*lum
    
    lradedd=lrad/self.eddington(mbh)[1]        #in erg/s
    
    return lrad, lradedd
    
  def lamb(self,z,f,foiii,mbh):
    
    lr=self.radiative(z,foiii,mbh)[1]
    lm=self.mechanical(z,f,mbh)[1]
    ledd=self.eddington(mbh)[1]

  
    lamb=(lm+lr)
  
  
  
    loglam=np.log10(lamb)
    
    return lamb, loglam
    

  def bondi (self,z,f,mbh):
    
    lr=self.cc.luminosity(z,f)
    
    pjet=lr/self.l0
    pjet=(pjet**0.71)
    print lr

    pjet=1e44*pjet

    print pjet
    #pjet=1.8e44/1e43


    # from Balmaverde et al. 2008 Pjet => Pbondi = Mbondi*c^2
                                # Pjet= 0.012 Mbondi * c^2
    
    pbondi=(np.log10(pjet/1e43)+1.91)/1.1   #the formula is in 1e43 units
    pbondi=pow(10,pbondi)
    pbondi=pbondi*1e43

    mbondi=pbondi/((self.C**2)*1)*self.myr

    mbondijet=pjet/((self.C**2)*0.012)*self.myr
   
    # from Allen et al. 2006 Pjet => Pbondi = 0.1 * Mbondi *c^2
    pallen=(np.log10(pjet/1e43)*0.77)+0.65
    pallen=pow(10,pallen)
    pallen=pallen*1e43
    
    mbondia=pallen/((self.C**2)*0.1)*self.myr

    
    #bondi and mass accretion rate from BH mass (Willett 2010)
  
    #determine c_s
    kbt=0.7e3*self.ev
    cs=np.sqrt((5/3*kbt)/(0.62*self.mp))
    
    #determine the density of the gas
    ne=1
    rho=ne*1.13*self.mp

  
    #pbondi=4*pi lambda G^2 M_Bh^2 rho eff c^2 cs^-3
    
    #pbondi = 2* G M_Bh cs^-2
    
    ra=2*self.G*mbh*self.msun/(cs**2)
  
    lambdaa=0.25
    eff=0.1

    bondi=lambdaa*np.pi*cs*ra**2*rho*self.myr
    
    pbonditop=0.017*bondi/self.myr*self.C**2

    #bondiratio=bondi/self.myr*(self.c**2)*self.eff
    #print bondiratio

    #general relativity
    #you confront the mass accretion rate with the eddington luminosity
    
    ledd=self.eddington(mbh)[1]
    ledduni=ledd/self.C**2*self.myr
    
    accr=bondi/ledduni
    print accr
    
    ##about PKS1718-649
    #infall
    dist=100*self.pc*1e-5
    vcloud=70
    mcloud=1e4
    
    #pjet allen
    dl=self.cc.lum_dist(z)

    minorax=0.4*1e-3*dl  #cm
    majorax=1.*1e-3*dl    #cm
    
    gamma=5./3./(5./3.-1.)
    T=1. #kev
    
    press=kbt*T*rho/(0.62*self.mp)
    print press
    
    V=np.pi*minorax*majorax**2
    print V

    tage=1e2*self.yr
    print tage
    
    pjetta=gamma*press*V/(tage*self.yr)

    t=dist/vcloud
    accrcloud=mcloud*self.msun/t*self.myr
    
    # from Allen et al. 2006 Pjet => Pbondi = 0.1 * Mbondi *c^2
    pallen2=(np.log10(pjetta/1e43)*0.77)+0.65
    pallen2=np.power(10,pallen2)
    pallen2=pallen2*1e43
    
    mbondia2=pallen2/((self.c**2)*1)*self.myr

    return pjet,pbondi,mbondi,mbondijet,mbondia,bondi,accr,accrcloud,pjetta,mbondia2,pbonditop
    
    
  def Xaccretion (self,z,f,mbh):
    
    ledd=self.eddington(mbh)[1]
    lx=self.cc.luminosity(z,f)
    lx=1e41
    lxratio=np.log10(lx/ledd)

    return lxratio

  def main(self):

    in1='''\t... Here's what I can do: ...\n'''
    in2='''\n\t - Eddington\t->\tE
\t - Mech. Lum.\t->\tLM
\t - Rad. Lum.\t->\tLR
\t - Lambda\t->\tLA
\t - Bondi\t->\tB
\t - X-accretion\t->\tX\n
      '''

    inp=str(raw_input(in1+in2))
  
    if (inp=='e' or inp=='E' or inp=='eddington'):
        in1='''\nM_BH [M_sun] = '''
        m=self.eddington(float(raw_input(in1)))
        print('''\nMedd =\t %g\t M_sun/yr'''% (m[0]))
        print('''Ledd =\t %g\t erg/s'''% (m[1]))
        print('''Ledd =\t %g\t M_sun/yr'''% (m[2]))
    elif (inp=='lm' or inp=='Mechanical Luminosity' or inp=='LM'):
        in1='''\nz =\t '''
        in2='''F [Jy] =\t '''
        in3='''M_BH [M_sun] =\t'''
        z=float(raw_input(in1))
        f=float(raw_input(in2))
        mbh=float(raw_input(in3))
        lmech=self.mechanical(z,f,mbh)
        print('''\nL_mech =\t %g\t erg/s\n'''% (lmech[0]))
        print('''L_mech/Ledd=\t %g\t '''%(lmech[1]))
    elif (inp=='lr' or inp=='Radiative Luminosity' or inp=='LR'):
        in1='''\nz =\t '''
        in2='''F_OIII [erg/scm2] =\t '''
        in3='''M_BH (M_sun) =\t '''
        z=float(raw_input(in1))
        f=float(raw_input(in2))
        mbh=float(raw_input(in3))
        lmech=self.radiative(z,f,mbh)
        print('''\nL_rad =\t %g\t erg/s\n'''% (lmech[0]))
        print('''L_rad/Ledd =\t %g '''%(lmech[1]))
    elif (inp=='LA' or inp=='lambda' or inp=='la'):
        in1='''\nz =\t '''
        in2='''F_OIII [erg/scm2] =\t '''
        in3='''F_1.4 [Jy] =\t '''
        in4='''M_BH [M_sun]=\t '''
        z=float(raw_input(in1))
        foiii=float(raw_input(in2))
        f=float(raw_input(in3))
        mbh=float(raw_input(in4))
        lmech=self.lamb(z,f,foiii,mbh)
        print('''\nLambda =\t %g \n'''% (lmech[0]))
        print('''LOGLambda =\t %g'''% (lmech[1]))
    elif (inp=='b' or inp=='Bondi' or inp=='B'):
        in1='''\nz =\t '''
        in2='''F [Jy] =\t '''
        in3='''M_BH [M_sun] =\t '''
        z=float(raw_input(in1))
        f=float(raw_input(in2))
        mbh=float(raw_input(in3))
        pj=self.bondi(z,f,mbh)
        print('''\tP_jet =\t %g\t erg/s\n'''% (pj[0]))
        print('''P_Bondi =\t %g\t erg/s\n'''% (pj[1]))
        print('''M_Bondi(B) =\t %g\t M_sun/yr\n'''% (pj[2]))
        print('''M_Bondi_jet(B) =\t %g\t M_sun/yr\n'''% (pj[3]))
        print('''M_Bondi(A) =\t %g\t M_sun/yr\n'''% (pj[4]))
        print('''M_bondi(BH) =\t %g\n'''% (pj[5]))
        print('''eff=\t %g\n'''% (pj[6]))
        print('''cloud=\t %g\n'''% (pj[7]))
        print('''pjet=\t %g\n'''% (pj[8]))
        print('''pallent=\t %g\n'''% (pj[9]))
        print('''pbonditop=\t %g\n'''% (pj[10]))
    elif (inp=='x' or inp=='Xaccretion' or inp=='X-ray'):
        in1='''\nz =\t '''
        in2='''F [Jy] =\t '''
        in3='''M_BH [M_sun] =\t '''
        z=float(raw_input(in1))
        f=float(raw_input(in2))
        mbh=float(raw_input(in3))
        pj=self.Xaccretion(z,f,mbh)
        print '''Xratio =\t %g \n'''% (pj)
    else:
        print ('\n\t ... you have not entered correct function ... \n')
        print ('\t************* --- AGN : ERROR --- **************\n')
        sys.exit(0)

    print ('\n\t************* --- AGN : DONE --- **************\n')