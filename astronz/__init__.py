# Import modules
import os
import sys
import string
import numpy as np


import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap


# get AstronZ install directory
ASTRONZ_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASTRONZ_DIR = ASTRONZ_PATH+'/astronz/'
sys.path.append(os.path.join(ASTRONZ_PATH, 'astronz'))

import agn
import cosmo
import radiohi

import pkg_resources
try:
    __version__ = pkg_resources.require("AstronZ")[0].version
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


    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    parser = ArgumentParser(description='AstronZ: tools to analyse astronomical data'
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
        action='store_true',
        help= 'tools for cosmological calculations')

    add('-a', '--agn',
        action='store_true',
        help='tools for AGN science')

    add('-hi', '--radioHI',
        action='store_true',
        help='''tools for neutral hydrogen science''')

    args = parser.parse_args(argv)

    if args.help:
        print '\n\t************* --- AstronZ : Help --- **************\n'

        print ('\t\t  ... called for help ...\n')
        parser.print_help()

        print ("""\nRun a command. This can be:\n
astronz\t\t(all tools)
astronz -c\t(cosmological tools)
astronz -hi\t(neutral hydroge tools)
astronz -agn \t(AGN science tools)
            """)
        print '\n\t************* --- AstronZ : DONE --- **************\n'

        sys.exit(0)

    elif args.cosmo:

        print ('\n\t************* --- AstronZ : Cosmo --- **************\n')
        print ('\t\t    ... Cosmological Tools ... \n')        
        c = cosmo.Cosmo()
        c.main()

    elif args.radioHI:
        print ('\n\t************* --- AstronZ : HI --- **************\n')
        print ('\t\t... Neutral Hydrogen Tools ... \n')
        hi = radiohi.radioHI()

    elif args.agn:
        print ('\n\t************* --- AstronZ : AGN --- **************\n')
        print ('\t\t   ... AGN Science tools ... \n')
        a = agn.AGN()
        a.main()
    else:
        print '\n\t************* --- AstronZ --- **************\n'
        in1="\t   ... list of the avaliable classes: ...\n"
        in2='''\n\t - c (cosmological tools)
\t - hi (neutral hydrogen tools)
\t - a (AGN science tools)\n
'''
        inp=str(raw_input(in1+in2))

        if inp == 'c':
            print ('\n\t************* --- AstronZ : Cosmo --- **************\n')
            print ('\t\t    ... Cosmological Tools ... \n')        
            c = cosmo.Cosmo()
            c.main()
        elif inp == 'a':
            print ('\n\t************* --- AstronZ : AGN --- **************\n')
            print ('\t\t   ... AGN Science tools ... \n')
            a = agn.AGN()
            a.main()
        elif inp == 'hi':
            print ('\n\t************* --- AstronZ : HI --- **************\n')
            print ('\t\t... Neutral Hydrogen Tools ... \n')
            hi = radiohi.radioHI()
            hi.main()
        else:
            print ('\n\t ... you have not entered an available class function ... \n')
            print ('\t************* --- AstronZ : ERROR --- **************\n')
            sys.exit(0)