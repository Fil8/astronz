#!/usr/bin/env python

import sys,string,os,math,numpy
import cosmo

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

    self.cc = cosmo.Cosmo()


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

    dl=self.cc.lum_dist(z)*1e-2
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
    
    M=self.abs_mag(m,z)
    T=self.tully(M,z)
    
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
    
    tau=self.tau_abs(scont,sabs)
    nhi=self.column(tau,dv)
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

  def main(self):

    in1='''\t... Here's what I can do: ...\n'''
    in2='''\n\t- HI at z\t->\tHI
\t- Tully-Fisher\t->\tTF
\t- TF-apparent\t->\tTFM
\t- Absolute Magnitude\t->\tM
\t- Optical depth\t->\tTAU
\t- Abs. Column Density\t->\tNHI
\t- fluxColumn Density\t->\tFHI
\t- EMColumn Density\t->\tEHI
\t- Mass HI\t->\tMHI
\t- Mass HI\t->\tCHI
\t- Mimimum_mass\t->\tMM
\t- Radio Power\t->\tRP
\t- Vel. at z(HI)\t->\tVEL
\t- Incl. deVAB\t->deVAB
\t- Vel. res.\t->\tVRES
        '''
  
    inp=str(raw_input(in1+in2))
    
    if (inp=='hi' or inp=='HI at z' or inp=='HI'):
        in1='''\nz =\t '''
        dl=self.hiatz(float(raw_input(in1)))
        print('''HI =\t %g\t MHz'''% (dl))
    elif(inp=='tf' or inp=='Tully-Fisher' or inp=='TF'):
        in1='''\nM [K20] =\t '''
        in2='''z =\t '''
        M=float(raw_input(in1))
        z=float(raw_input(in2))
        dl=self.tully(M,z)
        print('''\nV_flat =\t %g\t km/s'''% (dl))
    elif(inp=='tfm' or inp=='TF-apparent' or inp=='TFM'):
        in1='''\nm [K20] =\t'''
        in2='''z =\t '''
        m=float(raw_input(in1))
        z=float(raw_input(in2))
        dl=self.tully_app(m,z)
        print('''\nV_flat =\t %g\t km/s'''% (dl))
    elif(inp=='rp' or inp=='RadPow' or inp=='RP'):
        in1='''\nScont [Jy/beam] =\t'''
        in2='''z =\t '''
        scont=float(raw_input(in1))
        z=float(raw_input(in2))
        p=self.RadPow(z,scont)
        print('''\nRadio_Power =\t %g\t W/Hz'''% (p))
    elif(inp=='vel' or inp=='VEL' or inp=='V'):
        in1='''\nz =\t '''
        z=float(raw_input(in1))
        p=self.Velo(z)
        print('''velocity = %g km/s'''% (p) )
    elif(inp=='deVAB' or inp=='VAB' or inp=='deV'):
        in1='''\nAB =\t '''
        ratio=float(raw_input(in1))
        result=numpy.arccos(ratio)/3.14159265*180.
        print '''\ninclination =\t %g\t degrees'''% (result)
    elif(inp=='M' or inp=='absmag' or inp=='m'):
        in1='''\nm =\t '''
        in2='''z =\t '''
        m=float(raw_input(in1))
        z=float(raw_input(in2))
        M=self.abs_mag(m,z)
        print('''M =\t %g '''% (M))
    elif (inp=='tau' or inp=='Optical Depth' or inp=='TAU'):
        in1='''\nScont [Jy/beam] =\t '''
        in2='''Sabs [-Jy/beam] = '''
        scont=float(raw_input(in1))
        sabs=float(raw_input(in2))
        dl=self.tau_abs(scont,sabs)
        print '''tau =\t %g '''% (dl)
    elif(inp=='nhi' or inp=='ABSColumn Density' or inp=='NHI'):
        in1='''\ntau =\t '''
        in3='''FWHM [km/s] =\t '''
        n=float(raw_input(in1))
        nn=float(raw_input(in3))
        dl=self.column(n,nn)
        print('''N(HI) =\t %g\t cm^-2'''% (dl))
    elif(inp=='ehi' or inp=='EMColumn Density' or inp=='EHI'):
        in1='''\nSem [mJy/beam] =\t '''
        in2='''Beamx [arcsec] =\t '''
        in3='''Beamy [arcsec] =\t '''
        in4='''DeltaV [km/s] =\t '''
        s=float(raw_input(in1))
        bx=float(raw_input(in2))
        by=float(raw_input(in3))
        dv=float(raw_input(in4))
        dl=self.column_em(s,bx,by,dv)
        print('''N(HI) =\t %g\t cm^-2'''% (dl))
    elif (inp=='mhi' or inp=='Mass HI' or inp=='MHI'):
        in1='''\nScont [Jy/beam] =\t '''
        in2='''Sabs [-Jy/beam] =\t '''
        in3='''FWHM [km/s] =\t '''
        in4='''R radiosource [pc] =\t '''
        scont=float(raw_input(in1))
        sabs=float(raw_input(in2))
        dv=float(raw_input(in3))
        a=float(raw_input(in4))
        dl=self.mass_hi(scont,sabs,a,dv)
        print('''M =\t %g\t Msun'''% (dl))
    elif (inp=='chi' or inp=='Mass HI' or inp=='CHI'):
        in1='''\nSint (Jy/beam*km/s imstat] =\t '''
        in2='''Beamx [arcsec] =\t '''
        in3='''Beamy [arcsec] =\t '''
        in4='''pix [arcsec] =\t '''
        in5='''z =\t '''
        sint=float(raw_input(in1))
        bx=float(raw_input(in2))
        by=float(raw_input(in3))
        pix=float(raw_input(in4))
        z=float(raw_input(in5))
        dl=self.mass_hi_1(sint,bx,by,pix,z)
        print('''M = %g Msun'''% (dl))
    elif (inp=='fhi' or inp=='fluxColumn Density' or inp=='FHI'):
        in1='''\nScont [Jy/beam]=\t '''
        in2='''Sabs [-Jy/beam]=\t '''
        in3='''FWHM [km/s]=\t '''
        scont=float(raw_input(in1))
        sabs=float(raw_input(in2))
        dv=float(raw_input(in3))
        tau=self.tau_abs(scont,sabs)
        nhi=self.column(tau,dv)
        print('''tau =\t %g \n'''% (tau))
        print('''N(HI) =\t %g\t cm^-2'''% (nhi))
    elif (inp=='mm' or inp=='Minimum mass' or inp=='MM'):
        in1='''\nNHI [cm^-2] =\t '''
        in2='''Beamx [arcsec] =\t '''
        in3='''Beamy [arcsec] =\t '''
        in4='''z =\t'''
        nhi=float(raw_input(in1))
        bx=float(raw_input(in2))
        by=float(raw_input(in3))
        z=float(raw_input(in4))
        dl=self.mimimum_mass(nhi,bx,by,z)
        print('''M =\t %g\t Msun'''% (dl))
    elif (inp=='vres' or inp=='vel_res' or inp=='VRES'):
        in1='''\nDeltaLambda [Hz] =\t '''
        in2='''Lambda [Hz] =\t '''
        delta_lambda=float(raw_input(in1))
        lamb=float(raw_input(in2))
        deltav=self.vel_res(delta_lambda,lamb)
        print('''Vel Res =\t %g\t km/s'''% (deltav))
    else:
        print ('\n\t ... you have not entered correct function ... \n')
        print ('\t************* --- HI : ERROR --- **************\n')
        sys.exit(0)

    print ('\n\t************* --- HI : DONE --- **************\n')