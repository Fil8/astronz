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

def main (argv):

    print '\n\t************* --- ASTRONZ --- **************\n'

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
        help= 'tools for cosmological calculations')

    add('-a', '--agn',
        #type=str,
        default=False,
        action='store_true',
        help='tools for AGN science')

    add('-hi', '--radioHI',
        #type=bool,
        default=False,
        action='store_true',
        help='''tools for neutral hydrogen science''')

    args = parser.parse_args(argv)


# Get the function to execute for the command
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


