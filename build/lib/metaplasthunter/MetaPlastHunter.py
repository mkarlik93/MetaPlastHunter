#!/usr/bin/env python

###############################################################################
#                                                                             #
#    Entry point for MetaPlastHunter                                          #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Michal Karlicki'
__copyright__ = 'Copyright 2018'
__credits__ = ['Michal Karlicki']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Michal Karlicki'
__email__ = 'michal.karlicki@gmail.com'
__status__ = 'Development'

import logging
import argparse
logger = logging.getLogger("MetaPlastHunter")
logging.basicConfig(level=logging.INFO)
import argparse
from time import gmtime, strftime
from bin.bbmap_wrapper_v_1 import Mapping_runner, SAM2coverage, RapidRunner
from bin.taxonomic_assignment_v_1 import Taxonomic_assignment_Runner
from bin.settings import Settings_loader_yaml
import multiprocessing as mp
import os
import sys

class Run:

    """


    Main class for running MetaPlastHunter RC

    MPH for read classification



    """

    def __init__(self,input,input2,output,settings,threads):

        self.list_sra = ""
        self.station_name = ""
        self.settings = settings
        self.threads = threads
        self.input = input
        self.input2 = input2
        self.output = output

    def assign_taxnomomy_to_SAM(self):

        logger.info( " [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
#        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info( " [%s] Starting assign taxa to SAM file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        SAM2coverage(self.input, self.output, self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.input, self.output,self.settings).process()
    #Tests are still needed

    def rapid_taxonomic_assignment(self):
        #Do przerobienia (?)

        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
#        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        RapidRunner(self.input,self.input2,self.output,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.input, self.output,self.settings).process()

    def taxonomic_assigment(self):
        #Ta tez
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file_classification()

        Coverage_utillities(self.settings).add_empirical_treshold()

        logger.info("     [%s] Mapping and generating SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Mapping_runner(self.input,self.input2,self.output, self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.input,self.output, self.settings).process()



def main():

    description = """

Version %s


MetaPlastHunter -


 The efficient and accurate plastid reads classification pipeline.

Quantitative aproach for eukaryotic metagenomics.

Available workflows:

[--taxonomic_classification/-C] Searching, Classification, Visualization

[--rapid_classification, -Acc] Use it to lunch pipeline with in exact k-mer matching preliminary classification

[--sam_assign, -A]   Sequence alignment file (SAM) classification

Obligatory arguments:

Settings - inpute file

Facultative arguments:

threads

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com

Please cite:

https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/


This sofware was written by %s.
""" % (__version__,__author__)

    epilog = """


"""

    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('--taxonomic_classification','-C',action='store_true')
    parser.add_argument('--rapid_classification', '-Acc',action='store_true')
    parser.add_argument('--settings','-S', metavar='settings', type=str)
    parser.add_argument('--sam_assign','-A',action='store_true')
    parser.add_argument('--in_1', metavar='input',type=str)
    parser.add_argument('--in_2',nargs='?',type=str, default="")
    parser.add_argument('--output','-O',type=str)
    parser.add_argument('--threads','-T',nargs='?', type=int,default=mp.cpu_count())

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())


    process = Run(os.path.abspath(args.in_1),os.path.abspath(args.in_2), args.output, args.settings,args.threads)

    if args.taxonomic_classification:

        process.taxonomic_assigment()

    elif args.rapid_classification:

        process.rapid_taxonomic_assignment()

    elif args.sam_assign:

        process.assign_taxnomomy_to_SAM()

    else:

        logger.error("      Please specify pipeline")
        sys.exit()

    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)



if __name__ == "__main__":

    main()
