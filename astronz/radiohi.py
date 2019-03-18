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

# elif (inp=='RC' or inp=='rc' or inp=='RadioCalc'):
#   in1="\nHere's what I can do:\n"
#   in2="""*HI at z(HI)\t*Tully-Fisher(TF)\n*TF-apparent(TFM)\t*Absolute Magnitude(M)\n*Optical depth(TAU)\t*ABSColumn Density(NHI)\n*fluxColumn Density(FHI)\t*EMColumn Density(EHI)\n*Mass HI(MHI)\t *Mass HI(CHI)\t*Mimimum_mass(MM)\n*Radio Power(RP)\t*VEL at z(HI)\t*Incl deVAB\t\n*Velocity resolution (VRES)\n"""
  
#   inp=str(raw_input(in1+in2))
    
#   if (inp=='hi' or inp=='HI at z' or inp=='HI'):
#       in1="\nz= "
#       dl=r.hiatz(float(raw_input(in1)))
#       print "HI = %g MHz"% (dl)
#   elif(inp=='tf' or inp=='Tully-Fisher' or inp=='TF'):
#       in1="\nM(K20)= "
#       in2="z= "
#       M=float(raw_input(in1))
#       z=float(raw_input(in2))
#       dl=r.tully(M,z)
#       print "V_flat = %g km/s"% (dl)
#   elif(inp=='tfm' or inp=='TF-apparent' or inp=='TFM'):
#       in1="\nm(K20)= "
#       in2="z= "
#       m=float(raw_input(in1))
#       z=float(raw_input(in2))
#       dl=r.tully_app(m,z)
#       print "V_flat = %g km/s"% (dl)
#   elif(inp=='rp' or inp=='RadPow' or inp=='RP'):
#       in1="\nScont (Jy/beam)= "
#       in2="z= "
#       scont=float(raw_input(in1))
#       z=float(raw_input(in2))
#       p=r.RadPow(z,scont)
#       print "Radio_Power = %g W/Hz"% (p)
#   elif(inp=='vel' or inp=='VEL' or inp=='V'):
#       in1="z= "
#       z=float(raw_input(in1))
#       p=r.Velo(z)
#       print "velocity = %g km/s"% (p)  
#   elif(inp=='deVAB' or inp=='VAB' or inp=='deV'):
#       in1="AB= "
#       ratio=float(raw_input(in1))
#       result=numpy.arccos(ratio)/3.14159265*180.
#       print "inclination = %g degrees"% (result)
#   elif(inp=='M' or inp=='absmag' or inp=='m'):
#       in1="\nm= "
#       in2="z= "
#       m=float(raw_input(in1))
#       z=float(raw_input(in2))
#       M=r.abs_mag(m,z)
#       print "M = %g "% (M)
#   elif (inp=='tau' or inp=='Optical Depth' or inp=='TAU'):
#       in1="\nScont (Jy/beam)= "
#       in2="Sabs (-Jy/beam)= "
#       scont=float(raw_input(in1))
#       sabs=float(raw_input(in2))
#       dl=r.tau_abs(scont,sabs)
#       print "tau = %g "% (dl)
#   elif(inp=='nhi' or inp=='ABSColumn Density' or inp=='NHI'):
#       in1="\ntau= "
#       in3="FWHM (km/s)= "
#       n=float(raw_input(in1))
#       nn=float(raw_input(in3))
#       dl=r.column(n,nn)
#       print "N(HI) = %g cm^-2"% (dl)
#   elif(inp=='ehi' or inp=='EMColumn Density' or inp=='EHI'):
#       in1="\nSem (mJy/beam)= "
#       in2="Beamx (arcsec)= "
#       in3="Beamy (arcsec)= "
#       in4="DeltaV (km/s)= "
#       s=float(raw_input(in1))
#       bx=float(raw_input(in2))
#       by=float(raw_input(in3))
#       dv=float(raw_input(in4))
#       dl=r.column_em(s,bx,by,dv)
#       print "N(HI) = %g cm^-2"% (dl)
#   elif (inp=='mhi' or inp=='Mass HI' or inp=='MHI'):
#       in1="\nScont (Jy/beam)= "
#       in2="Sabs (-Jy/beam)= "
#       in3="FWHM (km/s)= "
#       in4="R radiosource (pc)= "
#       scont=float(raw_input(in1))
#       sabs=float(raw_input(in2))
#       dv=float(raw_input(in3))
#       a=float(raw_input(in4))
#       dl=r.mass_hi(scont,sabs,a,dv)
#       print "M = %g Msun"% (dl)
#   elif (inp=='chi' or inp=='Mass HI' or inp=='CHI'):
#       in1="\nSint (Jy/beam*km/s imstat)= "
#       in2="Beamx (arcsec)= "
#       in3="Beamy (arcsec)= "
#       in4="pix (arcsec)= "
#       in5="z = "
#       sint=float(raw_input(in1))
#       bx=float(raw_input(in2))
#       by=float(raw_input(in3))
#       pix=float(raw_input(in4))
#       z=float(raw_input(in5))
#       dl=r.mass_hi_1(sint,bx,by,pix,z)
#       print "M = %g Msun"% (dl)
#   elif (inp=='fhi' or inp=='fluxColumn Density' or inp=='FHI'):
#       in1="\nScont (Jy/beam)= "
#       in2="Sabs (-Jy/beam)= "
#       in3="FWHM (km/s)= "
#       scont=float(raw_input(in1))
#       sabs=float(raw_input(in2))
#       dv=float(raw_input(in3))
#       tau=r.tau_abs(scont,sabs)
#       nhi=r.column(tau,dv)
#       print "tau = %g \n"% (tau)
#       print "N(HI) = %g cm^-2"% (nhi)
#   elif (inp=='mm' or inp=='Minimum mass' or inp=='MM'):
#       in1="\nNHI (cm^-2)= "
#       in2="Beamx (arcsec)= "
#       in3="Beamy (arcsec)= "
#       in4="z ="
#       nhi=float(raw_input(in1))
#       bx=float(raw_input(in2))
#       by=float(raw_input(in3))
#       z=float(raw_input(in4))
#       dl=r.mimimum_mass(nhi,bx,by,z)
#       print "M = %g Msun"% (dl)
#   elif (inp=='vres' or inp=='vel_res' or inp=='VRES'):
#       in1="\nDeltaLambda (Hz)= "
#       in2="Lambda (Hz) = "
#       delta_lambda=float(raw_input(in1))
#       lamb=float(raw_input(in2))
#       deltav=r.vel_res(delta_lambda,lamb)
#       print "Vel Res = %g km/s"% (deltav)
#   else:
#       print 'you have not entered correct function\n EXIT FAILURE MOROAN'

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
    