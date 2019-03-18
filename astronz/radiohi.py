#!/usr/bin/env python

# ASTRONZ: astro calc utilities


import sys,string,os,math,numpy

RAD2DEG=180./math.pi
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

    
class radioHI():
  
  def __init__(self):

    print '*************ASTRONZ**************'
    in1="What can I do for you?\n"
    in2="CosmoCalc (CC)\tRadioCalc(RC)\tAGN(A)\n"
    
    inp=str(raw_input(in1+in2))
    self.hi=1.42040575177e+09 
    self.m26=-23.33
    self.s_tully=-9.64
    self.nhi=1.8216E18
    self.nhiem=3.1E17
    self.T=TSPIN
    self.mhi=MHI
    self.msun=MSUN
    self.chi=CHI
    self.mp=MP

  def RadPow(self,z,scont):

    dl=c.lum_dist(z)*1e-2
    print dl
    power=math.log10(pow(10,-26)*scont*4*math.pi*(dl**2))
    


    return power

  def vel_res(self,delta_lambda,wavelength):

    delta_v = delta_lambda/wavelength*C*1e-2

    delta_v /=1e3 #in km/s

    return delta_v 

  def Velo(self,z):

    freq=self.hi/(1+z)
    velocity=C*((self.hi-freq)/freq)/1e5

      
    return velocity

  def hiatz(self,z):
    hobs = self.hi/(1+z)/1e06   #MHz
    
    return hobs
    
  def tully(self,M,z):
    
    vflat= 10**((M-self.m26)/self.s_tully+2.6)/2
    
    return vflat
    
  def tully_app(self,m,z):
    
    M=r.abs_mag(m,z)
    T=r.tully(M,z)
    
    return T

  def abs_mag(self,m,z):
    dl=c.lum_dist(z)
    
    M=m-5*(math.log10(dl/3.085678e18)-1)
    
    return M

  def tau_abs(self,scont,sabs):
    
    tau=math.log(1-(sabs/scont))
    
    return tau
    
  def column(self,tau,dv):
    
    nhi=self.nhi*self.T*tau*dv
    
    return nhi
        
  def column_em(self,s,bx,by,dv):
    
    bxa=bx/60
    bya=by/60   #convert beam in arcmin
    nhi=self.nhiem*dv*s/(bxa*bya)
    
    return nhi
    
  def mass_hi(self,scont,sabs,a,dv):
    
    tau=r.tau_abs(scont,sabs)
    nhi=r.column(tau,dv)
    print nhi
    area=((a*PC)**2)
    
    mhi=area*self.mp*nhi/self.msun
  
    
    return mhi
  
  def mass_hi_1(self,s,bx,by,pix,z):
    
    dl=c.lum_dist(z)/3.085678e24
    beamcorr=2.*math.pi*(bx*by/(2.35482**2))/(math.pow(pix,2))
    bla=s/beamcorr

    mhi=self.chi*(dl**2)*bla
    
    return mhi
    
  def mimimum_mass(self, nhi, bx,by,z):
    
    bxa=c.ang2lin(bx,z)*1e6
    bya=c.ang2lin(by,z)*1e6

    beamarea=(bxa*bya)*(PC**2)
    min_mass=nhi*self.mhi*beamarea/MSUN
    
    return min_mass
    