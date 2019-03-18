# Import modules
import os
import sys
import string
import numpy as np

import logging
import warnings

from astropy.io import fits, ascii
from astropy import units as u
#from astropy.time import Time, TimeDelta
#from astropy.table import Table, Column, MaskedColumn


# get rfinder install directory
RFINDER_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RFINDER_DIR = RFINDER_PATH+'/astronz/'

sys.path.append(os.path.join(RFINDER_PATH, 'astronz'))

import agn
import cosmo
import hi

c= cosmo.Cosmo()
r= hi.radioHI()
a= agn.AGN()

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


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file '%s' does not exist!" % arg)

    return arg

class rfinder:
    '''

    Class to investigate the RFI behaviour during observations

    '''

    C = 2.99792458e5  # km/s
    HI = 1.420405751e9  # Hz

    def __init__(self, file=None):
        '''

        Set self.logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml

        '''

        # set self.logger
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)


    def set_dirs(self):
        '''
     
        Sets directory strucure and filenames
        Creates directory abs/ in basedir+beam and subdirectories spec/ and plot/

        OUTPUT:
            tab : table of catalog
            flux : flux of continuum sources abouve set threshold
     
        '''

        key = 'general'

        self.workdir  = self.cfg_par[key].get('workdir', None)
        self.msfile = self.workdir + self.cfg_par[key].get('msname', None)[0]
        self.cfg_par[key]['msfullpath'] = self.msfile      

        self.outdir  = self.cfg_par[key].get('outdir', None)
        self.rfidir  = self.outdir+'rfi_'+self.cfg_par['rfi']['polarization']+'/'
        self.cfg_par[key]['rfidir'] = self.rfidir
        self.rfifile = self.rfidir+'rfi_flagged_vis.MS'
        self.rfi_freq_base = self.rfidir+'freq_base.fits'
        self.rfimsfile = self.rfidir+'rfi_flagged.MS'
        self.tabledir = self.rfidir+'tables/'
        self.cfg_par[key]['tabledir'] = self.tabledir
        self.rfi_table = self.tabledir+'rfi_table.fits'
    
        self.rfiplotdir = self.rfidir+'plots/'
        self.cfg_par[key]['plotdir'] = self.rfiplotdir 

 
        self.moviedir = self.rfidir+'plots/movies/'
        self.cfg_par[key]['moviedir'] = self.moviedir        

        if os.path.exists(self.moviedir) == False:
             os.makedirs(self.moviedir)

        if os.path.exists(self.rfidir) == False:
             os.makedirs(self.rfidir)           

        if os.path.exists(self.tabledir) == False:
             os.makedirs(self.tabledir)

        if os.path.exists(self.rfiplotdir) == False:
             os.makedirs(self.rfiplotdir)

        if self.cfg_par['rfi']['chunks']['time_enable'] == True:

            self.rfitimedir = self.rfidir+'time_chunks/'
            self.cfg_par[key]['rfitimedir'] = self.rfitimedir

            if os.path.exists(self.rfitimedir) == False:
                 os.makedirs(self.rfitimedir)

            self.timetabledir = self.tabledir+'time_chunks/'
            self.cfg_par[key]['timetabledir'] = self.timetabledir

            if os.path.exists(self.timetabledir) == False:
                 os.makedirs(self.timetabledir)

            timeplotdir_tmp = self.rfiplotdir+'time_chunks/'
            
            if os.path.exists(timeplotdir_tmp) == False:
                 os.makedirs(timeplotdir_tmp)

            self.timeplotdir1d = timeplotdir_tmp+'1D/'
            self.cfg_par[key]['timeplotdir1D'] = self.timeplotdir1d

            if os.path.exists(self.timeplotdir1d) == False:
                 os.makedirs(self.timeplotdir1d)

            self.timeplotdir2d = timeplotdir_tmp+'2D/'
            self.cfg_par[key]['timeplotdir2D'] = self.timeplotdir2d

            if os.path.exists(self.timeplotdir2d) == False:
                 os.makedirs(self.timeplotdir2d)

            self.altazplotdir = self.rfidir+'plots/altaz/'
            self.cfg_par[key]['altazplotdir'] = self.altazplotdir        

            if os.path.exists(self.altazplotdir) == False:
                 os.makedirs(self.altazplotdir)

    def read_args(self,args):

        if args.working_dir:
            self.cfg_par['general']['workdir'] = args.working_dir
        if args.output_dir:
            self.cfg_par['general']['outdir'] = args.output_dir
        if args.input:
            self.cfg_par['general']['msname'] = args.input
        if args.field:
            self.cfg_par['general']['field'] = args.field
        if args.polarization:
            self.cfg_par['rfi']['polarization'] = args.polarization
        if args.telescope:
            self.cfg_par['rfi']['telescope'] = args.telescope
        if args.baseline_cut:
            self.cfg_par['rfi']['baseline_cut'] = args.baseline_cut
        if args.time_step:
            self.cfg_par['rfi']['chunks']['time_enable'] = True
            self.cfg_par['rfi']['chunks']['time_step'] = args.time_step
        if args.spw_av:
            self.cfg_par['rfi']['chunks']['spw_enable'] = True
            self.cfg_par['rfi']['chunks']['spw_width'] = args.spw_av
        if (args.rfimode == 'rms_clip' or args.rfimode == 'rms') :
                self.cfg_par['rfi']['RFInder_mode'] = args.rfimode
                if args.sigma_clip:
                    self.cfg_par['rfi']['rms_clip'] = args.sigma_clip
                if args.frequency_interval:
                    self.cfg_par['rfi']['noise_measure_edges'] = args.frequency_interval


 


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


    def main (self,argv):

        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        parser = ArgumentParser(description='RFInder: package to visualize the flagged RFI in a dataset '
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

        add('-c', '--config',
            type=lambda a: is_valid_file(parser, a),
            default=False,
            help='RFInder configuration file (YAML format)')

        add('-w', '--working_dir',
            type= str,
            default = False,
            help= 'select working directory (MS file assumed to be here)')

        add('-odir', '--output_dir',
            type=str,
            default=False,
            help='select output directory')

        add('-i', '--input',
            type=str,
            default=False,
            help='''input ['MS'] file''')

        add('-fl', '--field',
            type=int,
            default=False,
            help='select field of MS file to analyze')

        add('-tel', '--telescope',
            type=str,
            default=False,
            help='select telescope: meerkat or apertif(WSRT)')

        add('-mode', '--rfimode',
            type=str,
            default=False,
            help='select mode where to investigate RFI: use_flags or rms_clip, rms')

        add('-pol', '--polarization',
            type=str,
            default=False,
            help='select stokes parameter')

        add('-fint', '--frequency_interval',
            nargs='*',
            default=False,
            help='select frequency interval where to measure noise')

        add('-spwAv', '--spw_av',
            type=int,
            default=False,
            help='select average in frequency')

        add('-tStep', '--time_step',
            type=int,
            default=False,
            help='select time step')

        add('-sig', '--sigma_clip',
            type=int,
            default=False,
            help='select sigma clip for rms_clip mode')

        add('-baseCut', '--baseline_cut',
            type=int,
            default=False,
            help='select baseline cut for differential RFI analysis')

        args = parser.parse_args(argv)

        if args.help:  #rfinder -h 
            parser.print_help()

            print ("""\nRun a command. This can be: \nrfinder \nrfinder -c path_to_config_file.yml\n rfinder -i ['ngc1399.ms']""")

            sys.exit(0)

        elif args.config:    #rfinder -c config_file.yml         
            self.logger.info('\t ... Reading your parameter file ... \n')
            # read database here
            files =  args.config
            cfg = open(files)
            self.cfg_par = yaml.load(cfg)

        else: #rfinder  or rfinder -options
            workdir = os.getcwd()
            exists = os.path.isfile(workdir+'/'+DEFAULT_CONFIG)
            if exists:
                self.logger.info('\t ... Reading default parameter file in your directory ... \n')
                file_default = os.path.join(workdir, DEFAULT_CONFIG)
                print file_default
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg)
            else:
                # Keep presets
                self.logger.info('\t ... Reading default installation parameter file ... \n')
                file_default = os.path.join(RFINDER_DIR, DEFAULT_CONFIG)
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg)            
                self.cfg_par['general']['workdir'] = workdir
                self.cfg_par['general']['outdir'] = workdir
                print workdir
                with open(workdir+'/'+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)

                if args.field ==False and args.input==False:

                    self.logger.info('''MSNAME & Field missing\n please edit rfinder_default.yml in your current directory\n
                    or run: rfinder -i ['msname'] -fl number (assuming the observation is located in your current directory)
                    \n''')
                    
                    sys.exit(0)

                else:
                    self.logger.info('''... you gave MSname and field in your first run, 
                        assuming they are in your current directory ...
                    ''')

            if all(x==False for x in vars(args).values()) == False and (args.help==False and args.config == False):
            
                self.logger.info('\t ... Updating arguments given from terminal ... \n')
                self.read_args(args)
                with open(workdir+'/'+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)
 

        self.cfg_par['general']['template_folder'] = os.path.join(RFINDER_PATH,'rfinder/templates')
        self.set_dirs()

        return self
