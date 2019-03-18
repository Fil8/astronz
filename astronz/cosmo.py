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


    print '*************ASTRONZ**************'
    in1="What can I do for you?\n"
    in2="CosmoCalc (CC)\tRadioCalc(RC)\tAGN(A)\n"

# if (inp=='CC' or inp=='cc' or inp=='CosmoCalc'):
#   in1="\nHere's what I can do:\n"
#   in2="*lumdist(DL)\t*scalearcsec(LA)\t*scaleMpc(LR)\t*ageatz(CA)\tluminosity(L)\n"
  
#   inp=str(raw_input(in1+in2))
  
#   if (inp=='dl' or inp=='lumdist' or inp=='DL'):
#       in1="\nz= "
#       dl=c.lum_dist(float(raw_input(in1)))/3.085678e24
#       print "D_L = %g Mpc"% (dl)
#   elif(inp=='la' or inp=='scalearcsec' or inp=='LA'):
#       in1="\nz= "
#       in2="R(kpc)= "
#       z=float(raw_input(in1))
#       r=float(raw_input(in2))/1e3   #since the module wants the radius in Mpc
#       print "R = %g arcsec"% (c.lin2ang(r,z))
#   elif(inp=='lr' or inp=='scaleMpc' or inp=='LR'):
#       in1="\nz= "
#       in2="R(arcsec)= "
#       z=float(raw_input(in1))
#       r=float(raw_input(in2))
#       print "R = %g kpc"% (c.ang2lin(r,z)*1e3)
#   elif(inp=='ca' or inp=='scalearcsec' or inp=='CA'):
#       in1="\nz= "
#       dl=c.compute_age(float(raw_input(in1)))
#       print "Age = %g Gyears"% (dl)
#   elif(inp=='l' or inp=='luminosity' or inp=='L'):
#       in1="\nz= "
#       in2="F(Jy)= "
#       z=float(raw_input(in1))
#       f=float(raw_input(in2))
#       a=(c.luminosity(z,f))
#       print "L = %g erg/s\n"% a
#       b=a*1e-7
#       print "L = %g W"% b
#       c=a*HI
#       print "L = %g W"% c
#   else:
#       print 'you have not entered correct function\n EXIT FAILURE MOROAN'


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
    
class radioHI():
  def __init__(self):
    self.hi=HI
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
    
if __name__=='__main__':
    
    c=Cosmo()
    r=radioHI()
    a=AGN()


    print '*************ASTRONZ**************'
    in1="What can I do for you?\n"
    in2="CosmoCalc (CC)\tRadioCalc(RC)\tAGN(A)\n"
    
    inp=str(raw_input(in1+in2))
    
    if (inp=='CC' or inp=='cc' or inp=='CosmoCalc'):
      in1="\nHere's what I can do:\n"
      in2="*lumdist(DL)\t*scalearcsec(LA)\t*scaleMpc(LR)\t*ageatz(CA)\tluminosity(L)\n"
      
      inp=str(raw_input(in1+in2))
      
      if (inp=='dl' or inp=='lumdist' or inp=='DL'):
          in1="\nz= "
          dl=c.lum_dist(float(raw_input(in1)))/3.085678e24
          print "D_L = %g Mpc"% (dl)
      elif(inp=='la' or inp=='scalearcsec' or inp=='LA'):
          in1="\nz= "
          in2="R(kpc)= "
          z=float(raw_input(in1))
          r=float(raw_input(in2))/1e3   #since the module wants the radius in Mpc
          print "R = %g arcsec"% (c.lin2ang(r,z))
      elif(inp=='lr' or inp=='scaleMpc' or inp=='LR'):
          in1="\nz= "
          in2="R(arcsec)= "
          z=float(raw_input(in1))
          r=float(raw_input(in2))
          print "R = %g kpc"% (c.ang2lin(r,z)*1e3)
      elif(inp=='ca' or inp=='scalearcsec' or inp=='CA'):
          in1="\nz= "
          dl=c.compute_age(float(raw_input(in1)))
          print "Age = %g Gyears"% (dl)
      elif(inp=='l' or inp=='luminosity' or inp=='L'):
          in1="\nz= "
          in2="F(Jy)= "
          z=float(raw_input(in1))
          f=float(raw_input(in2))
          a=(c.luminosity(z,f))
          print "L = %g erg/s\n"% a
          b=a*1e-7
          print "L = %g W"% b
          c=a*HI
          print "L = %g W"% c
      else:
          print 'you have not entered correct function\n EXIT FAILURE MOROAN'

          
    elif (inp=='RC' or inp=='rc' or inp=='RadioCalc'):
      in1="\nHere's what I can do:\n"
      in2="""*HI at z(HI)\t*Tully-Fisher(TF)\n*TF-apparent(TFM)\t*Absolute Magnitude(M)\n*Optical depth(TAU)\t*ABSColumn Density(NHI)\n*fluxColumn Density(FHI)\t*EMColumn Density(EHI)\n*Mass HI(MHI)\t *Mass HI(CHI)\t*Mimimum_mass(MM)\n*Radio Power(RP)\t*VEL at z(HI)\t*Incl deVAB\t\n*Velocity resolution (VRES)\n"""
      
      inp=str(raw_input(in1+in2))
        
      if (inp=='hi' or inp=='HI at z' or inp=='HI'):
          in1="\nz= "
          dl=r.hiatz(float(raw_input(in1)))
          print "HI = %g MHz"% (dl)
      elif(inp=='tf' or inp=='Tully-Fisher' or inp=='TF'):
          in1="\nM(K20)= "
          in2="z= "
          M=float(raw_input(in1))
          z=float(raw_input(in2))
          dl=r.tully(M,z)
          print "V_flat = %g km/s"% (dl)
      elif(inp=='tfm' or inp=='TF-apparent' or inp=='TFM'):
          in1="\nm(K20)= "
          in2="z= "
          m=float(raw_input(in1))
          z=float(raw_input(in2))
          dl=r.tully_app(m,z)
          print "V_flat = %g km/s"% (dl)
      elif(inp=='rp' or inp=='RadPow' or inp=='RP'):
          in1="\nScont (Jy/beam)= "
          in2="z= "
          scont=float(raw_input(in1))
          z=float(raw_input(in2))
          p=r.RadPow(z,scont)
          print "Radio_Power = %g W/Hz"% (p)
      elif(inp=='vel' or inp=='VEL' or inp=='V'):
          in1="z= "
          z=float(raw_input(in1))
          p=r.Velo(z)
          print "velocity = %g km/s"% (p)  
      elif(inp=='deVAB' or inp=='VAB' or inp=='deV'):
          in1="AB= "
          ratio=float(raw_input(in1))
          result=numpy.arccos(ratio)/3.14159265*180.
          print "inclination = %g degrees"% (result)
      elif(inp=='M' or inp=='absmag' or inp=='m'):
          in1="\nm= "
          in2="z= "
          m=float(raw_input(in1))
          z=float(raw_input(in2))
          M=r.abs_mag(m,z)
          print "M = %g "% (M)
      elif (inp=='tau' or inp=='Optical Depth' or inp=='TAU'):
          in1="\nScont (Jy/beam)= "
          in2="Sabs (-Jy/beam)= "
          scont=float(raw_input(in1))
          sabs=float(raw_input(in2))
          dl=r.tau_abs(scont,sabs)
          print "tau = %g "% (dl)
      elif(inp=='nhi' or inp=='ABSColumn Density' or inp=='NHI'):
          in1="\ntau= "
          in3="FWHM (km/s)= "
          n=float(raw_input(in1))
          nn=float(raw_input(in3))
          dl=r.column(n,nn)
          print "N(HI) = %g cm^-2"% (dl)
      elif(inp=='ehi' or inp=='EMColumn Density' or inp=='EHI'):
          in1="\nSem (mJy/beam)= "
          in2="Beamx (arcsec)= "
          in3="Beamy (arcsec)= "
          in4="DeltaV (km/s)= "
          s=float(raw_input(in1))
          bx=float(raw_input(in2))
          by=float(raw_input(in3))
          dv=float(raw_input(in4))
          dl=r.column_em(s,bx,by,dv)
          print "N(HI) = %g cm^-2"% (dl)
      elif (inp=='mhi' or inp=='Mass HI' or inp=='MHI'):
          in1="\nScont (Jy/beam)= "
          in2="Sabs (-Jy/beam)= "
          in3="FWHM (km/s)= "
          in4="R radiosource (pc)= "
          scont=float(raw_input(in1))
          sabs=float(raw_input(in2))
          dv=float(raw_input(in3))
          a=float(raw_input(in4))
          dl=r.mass_hi(scont,sabs,a,dv)
          print "M = %g Msun"% (dl)
      elif (inp=='chi' or inp=='Mass HI' or inp=='CHI'):
          in1="\nSint (Jy/beam*km/s imstat)= "
          in2="Beamx (arcsec)= "
          in3="Beamy (arcsec)= "
          in4="pix (arcsec)= "
          in5="z = "
          sint=float(raw_input(in1))
          bx=float(raw_input(in2))
          by=float(raw_input(in3))
          pix=float(raw_input(in4))
          z=float(raw_input(in5))
          dl=r.mass_hi_1(sint,bx,by,pix,z)
          print "M = %g Msun"% (dl)
      elif (inp=='fhi' or inp=='fluxColumn Density' or inp=='FHI'):
          in1="\nScont (Jy/beam)= "
          in2="Sabs (-Jy/beam)= "
          in3="FWHM (km/s)= "
          scont=float(raw_input(in1))
          sabs=float(raw_input(in2))
          dv=float(raw_input(in3))
          tau=r.tau_abs(scont,sabs)
          nhi=r.column(tau,dv)
          print "tau = %g \n"% (tau)
          print "N(HI) = %g cm^-2"% (nhi)
      elif (inp=='mm' or inp=='Minimum mass' or inp=='MM'):
          in1="\nNHI (cm^-2)= "
          in2="Beamx (arcsec)= "
          in3="Beamy (arcsec)= "
          in4="z ="
          nhi=float(raw_input(in1))
          bx=float(raw_input(in2))
          by=float(raw_input(in3))
          z=float(raw_input(in4))
          dl=r.mimimum_mass(nhi,bx,by,z)
          print "M = %g Msun"% (dl)
      elif (inp=='vres' or inp=='vel_res' or inp=='VRES'):
          in1="\nDeltaLambda (Hz)= "
          in2="Lambda (Hz) = "
          delta_lambda=float(raw_input(in1))
          lamb=float(raw_input(in2))
          deltav=r.vel_res(delta_lambda,lamb)
          print "Vel Res = %g km/s"% (deltav)
      else:
          print 'you have not entered correct function\n EXIT FAILURE MOROAN'
          
    elif (inp=='A' or inp=='a' or inp=='AGN'):
      in1="\nHere's what I can do:\n"
      in2="*Eddington(e)\t*Mechanical Luminosity (LM)\t*Radiative Luminosity(LR)\n*Lambda(LA)\t*Bondi(B)\t*X-accretion\n"
      
      inp=str(raw_input(in1+in2))
      
      if (inp=='e' or inp=='E' or inp=='eddington'):
          in1="\nM_BH(M_sun)= "
          m=a.eddington(float(raw_input(in1)))
          print "Medd = %g M_sun/yr\n"% (m[0])
          print "Ledd = %g erg/s"% (m[1])
          print "Ledd1 = %g M_sun/yr"% (m[2])
      elif (inp=='lm' or inp=='Mechanical Luminosity' or inp=='LM'):
          in1="\nz = "
          in2="F (Jy)= "
          in3="M_BH (M_sun)= "
          z=float(raw_input(in1))
          f=float(raw_input(in2))
          mbh=float(raw_input(in3))
          lmech=a.mechanical(z,f,mbh)
          print "L_mech = %g erg/s\n"% (lmech[0])
          print "L_mech/Ledd= %g "%(lmech[1])
      elif (inp=='lr' or inp=='Radiative Luminosity' or inp=='LR'):
          in1="\nz = "
          in2="F_OIII (erg/scm2)= "
          in3="M_BH (M_sun)= "
          z=float(raw_input(in1))
          f=float(raw_input(in2))
          mbh=float(raw_input(in3))
          lmech=a.radiative(z,f,mbh)
          print "L_rad = %g erg/s\n"% (lmech[0])
          print "L_rad/Ledd= %g "%(lmech[1])
      elif (inp=='LA' or inp=='lambda' or inp=='la'):
          in1="\nz = "
          in2="F_OIII (erg/scm2)= "
          in3="F_1.4 (Jy)= "
          in4="M_BH (M_sun)= "
          z=float(raw_input(in1))
          foiii=float(raw_input(in2))
          f=float(raw_input(in3))
          mbh=float(raw_input(in4))
          lmech=a.lamb(z,f,foiii,mbh)
          print "Lambda = %g \n"% (lmech[0])
          print "LOGLambda = %g \n"% (lmech[1])
      elif (inp=='b' or inp=='Bondi' or inp=='B'):
          in1="\nz = "
          in2="F (Jy)= "
          in3="M_BH (M_sun)= "
          z=float(raw_input(in1))
          f=float(raw_input(in2))
          mbh=float(raw_input(in3))
          pj=a.bondi(z,f,mbh)
          print "P_jet = %g erg/s\n"% (pj[0])
          print "P_Bondi = %g erg/s\n"% (pj[1])
          print "M_Bondi(B) = %g M_sun/yr\n"% (pj[2])
          print "M_Bondi_jet(B) = %g M_sun/yr\n"% (pj[3])
          print "M_Bondi(A) = %g M_sun/yr\n"% (pj[4])
          print "M_bondi(BH) = %g\n"% (pj[5])
          print "eff= %g\n"% (pj[6])
          print "cloud= %g\n"% (pj[7])
          print "pjet= %g\n"% (pj[8])
          print "pallent= %g\n"% (pj[9])
          print "pbonditop= %g\n"% (pj[10])
      elif (inp=='x' or inp=='Xaccretion' or inp=='X-ray'):
          in1="\nz = "
          in2="F (Jy)= "
          in3="M_BH (M_sun)= "
          z=float(raw_input(in1))
          f=float(raw_input(in2))
          mbh=float(raw_input(in3))
          pj=a.Xaccretion(z,f,mbh)
          print "Xratio = %g \n"% (pj)
    else:
      print 'you have not entered correct function\n EXIT FAILURE MOROAN'
