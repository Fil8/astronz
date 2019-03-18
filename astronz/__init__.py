# Import modules
import os
import sys
import string
import numpy as np

import logging
import warnings

import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap

from astropy.io import fits, ascii
from astropy import units as u
#from astropy.time import Time, TimeDelta
#from astropy.table import Table, Column, MaskedColumn


# get rfinder install directory
ASTRONZ_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASTRONZ_DIR = ASTRONZ_PATH+'/astronz/'

sys.path.append(os.path.join(ASTRONZ_PATH, 'astronz'))

import agn
import cosmo
import radiohi

#c= cosmo.Cosmo()
#r= hi.radioHI()
#a= agn.AGN()

#DEFAULT_CONFIG = 'rfinder_default.yml'

#if not sys.warnoptions:
#    warnings.simplefilter("ignore")


import pkg_resources

try:
    __version__ = pkg_resources.require("astronz")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"

####################################################################################################

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

#    '''

#    Class to investigate the RFI behaviour during observations

#    '''

#C = 2.99792458e5  # km/s
#HI = 1.420405751e9  # Hz

#def __init__(self):
#    '''

#    Set self.logger for spectrum extraction
#    Find config file
#    If not specified by user load rfinder_default.yml

#    '''

    # set self.logger
#    logging.basicConfig(level=logging.INFO)
#    self.logger = logging.getLogger(__name__)






# print '*************ASTRONZ**************'
# in1="What can I do for you?\n"
# in2="CosmoCalc (CC)\tRadioCalc(RC)\tAGN(A)\n"

# inp=str(raw_input(in1+in2))

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
      
# elif (inp=='A' or inp=='a' or inp=='AGN'):
#   in1="\nHere's what I can do:\n"
#   in2="*Eddington(e)\t*Mechanical Luminosity (LM)\t*Radiative Luminosity(LR)\n*Lambda(LA)\t*Bondi(B)\t*X-accretion\n"
  
#   inp=str(raw_input(in1+in2))
  
#   if (inp=='e' or inp=='E' or inp=='eddington'):
#       in1="\nM_BH(M_sun)= "
#       m=a.eddington(float(raw_input(in1)))
#       print "Medd = %g M_sun/yr\n"% (m[0])
#       print "Ledd = %g erg/s"% (m[1])
#       print "Ledd1 = %g M_sun/yr"% (m[2])
#   elif (inp=='lm' or inp=='Mechanical Luminosity' or inp=='LM'):
#       in1="\nz = "
#       in2="F (Jy)= "
#       in3="M_BH (M_sun)= "
#       z=float(raw_input(in1))
#       f=float(raw_input(in2))
#       mbh=float(raw_input(in3))
#       lmech=a.mechanical(z,f,mbh)
#       print "L_mech = %g erg/s\n"% (lmech[0])
#       print "L_mech/Ledd= %g "%(lmech[1])
#   elif (inp=='lr' or inp=='Radiative Luminosity' or inp=='LR'):
#       in1="\nz = "
#       in2="F_OIII (erg/scm2)= "
#       in3="M_BH (M_sun)= "
#       z=float(raw_input(in1))
#       f=float(raw_input(in2))
#       mbh=float(raw_input(in3))
#       lmech=a.radiative(z,f,mbh)
#       print "L_rad = %g erg/s\n"% (lmech[0])
#       print "L_rad/Ledd= %g "%(lmech[1])
#   elif (inp=='LA' or inp=='lambda' or inp=='la'):
#       in1="\nz = "
#       in2="F_OIII (erg/scm2)= "
#       in3="F_1.4 (Jy)= "
#       in4="M_BH (M_sun)= "
#       z=float(raw_input(in1))
#       foiii=float(raw_input(in2))
#       f=float(raw_input(in3))
#       mbh=float(raw_input(in4))
#       lmech=a.lamb(z,f,foiii,mbh)
#       print "Lambda = %g \n"% (lmech[0])
#       print "LOGLambda = %g \n"% (lmech[1])
#   elif (inp=='b' or inp=='Bondi' or inp=='B'):
#       in1="\nz = "
#       in2="F (Jy)= "
#       in3="M_BH (M_sun)= "
#       z=float(raw_input(in1))
#       f=float(raw_input(in2))
#       mbh=float(raw_input(in3))
#       pj=a.bondi(z,f,mbh)
#       print "P_jet = %g erg/s\n"% (pj[0])
#       print "P_Bondi = %g erg/s\n"% (pj[1])
#       print "M_Bondi(B) = %g M_sun/yr\n"% (pj[2])
#       print "M_Bondi_jet(B) = %g M_sun/yr\n"% (pj[3])
#       print "M_Bondi(A) = %g M_sun/yr\n"% (pj[4])
#       print "M_bondi(BH) = %g\n"% (pj[5])
#       print "eff= %g\n"% (pj[6])
#       print "cloud= %g\n"% (pj[7])
#       print "pjet= %g\n"% (pj[8])
#       print "pallent= %g\n"% (pj[9])
#       print "pbonditop= %g\n"% (pj[10])
#   elif (inp=='x' or inp=='Xaccretion' or inp=='X-ray'):
#       in1="\nz = "
#       in2="F (Jy)= "
#       in3="M_BH (M_sun)= "
#       z=float(raw_input(in1))
#       f=float(raw_input(in2))
#       mbh=float(raw_input(in3))
#       pj=a.Xaccretion(z,f,mbh)
#       print "Xratio = %g \n"% (pj)
# else:
#   print 'you have not entered correct function\n EXIT FAILURE MOROAN'


def main (argv):

    print '\t*************ASTRONZ**************\n'

    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    parser = ArgumentParser(description='ASTRONZ: tools to analyse astronomical data'
                            '|n version {:s} |n install path {:s} |n '
                            'Filippo Maccagni <filippo.maccagni@gmial.com>'.format(__version__,
                                                                               os.path.dirname(__file__)),
                            formatter_class=MultilineFormatter,
                            add_help=False)

    add = parser.add_argument

    add("-h", "--help",  action="store_true",
            help="Print help message and exit")

    add("-v","--version", action='version',
            version='{:s} version {:s}'.format(parser.prog, __version__))

    add('-c', '--cosmo',
        #type= str,
        default = False,
        action='store_true',
        help= 'select set of cosmological tools')

    add('-a', '--agn',
        #type=str,
        default=False,
        action='store_true',
        help='select set of tools for AGN science')

    add('-hi', '--radioHI',
        #type=bool,
        #default=False,
        action='store_true',
        help='''select set of tools for neutral hydrogen science''')

    args = parser.parse_args(argv)

    if args.help:  #rfinder -h 
        print ('\t... help: called for help ...\n')
        parser.print_help()

        print ("""\nRun a command. This can be:\n
astronz
astronz -cosmo (-c)
astronz -radioHI (-hi)
astronz -agn (-a)
            """)

        sys.exit(0)

    elif args.cosmo:    #rfinder -c config_file.yml         
        print ('\t... cosmo: Cosmological Tools ... \n')
        c = cosmo.Cosmo()

    elif args.radioHI:
        print ('\t... radioHI: Neutral Hydrogen Tools ... \n')
        hi = radiohi.radioHI()

    elif args.agn:
        print ('\t... AGN: AGN science Tools ... \n')
        a = agn.AGN()

