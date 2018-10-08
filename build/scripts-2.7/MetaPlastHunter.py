#!//anaconda/bin/python
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
from time import gmtime, strftime
import os
from bin.bbmap_wrapper import BBpipe
from bin.bbmap_wrapper import BBpipe_with_bbduk_preliminary
from bin.bbmap_wrapper_v_1 import Mapping_runner
from bin.krakenize import Pipeline_kraken
from bin.taxonomic_assignment_v_1 import Taxonomic_assignment_Runner
from bin.get_data import Pipeline_fetch
from bin.settings import Settings_loader_yaml
from bin.genome_reconstruction import Genome_reconstruction_pipe
import sys
import multiprocessing as mp


logger = logging.getLogger("MetaPlastHunter")
logging.basicConfig(level=logging.INFO)


class Run:

    """
    Main class for running MetaPlastHunter

    """

    def __init__(self,settings,threads):

        self.list_sra = ""
        self.station_name = ""
        self.settings = settings
        self.threads = threads

    def genomic_reconstruction(self):
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info( "     [%s] Starting genomic reconstruction " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Genome_reconstruction_pipe(self.settings).process()

    def full_wf(self):
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info( "     [%s] Downloading data from SRArchive" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_fetch(self.list_sra,self.station_name,self.settings).run()
        logger.info("     [%s] Preliminary classification" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_kraken(self.list_sra, self.station_name,self.settings,self.threads).run()
        logger.info("     [%s] BBmap initial mapping" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        BBpipe(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Taxonomic_assignment_Runner(self.list_sra, self.station_name,self.settings).process()

    #Osobny skrypt
    def fetch_wf(self):
        #Ta funkcja jest okej
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info( "     [%s] Downloading data from SRArchive" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_fetch(self.list_sra,self.station_name,self.settings).run()

    def classification_wf(self):

        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] Preliminary classification" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Pipeline_kraken(self.list_sra, self.station_name,self.settings,self.threads).run()
        logger.info("     [%s] BBmap initial mapping" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        BBpipe(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.list_sra, self.station_name,self.settings).process()

    def recalculation_wf(self):
        #Do wywalenia
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        BBpipe(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.list_sra, self.station_name,self.settings).process()


    #Some tests are still needed
    def accelerated_taxonomic_assignment(self):
        #Do przerobienia
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        BBpipe_with_bbduk_preliminary(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.list_sra, self.station_name,self.settings).process()

    def taxonomic_assigment(self):
        #Ta tez
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file_classification()
        logger.info("     [%s] Mapping and generating SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Mapping_runner(self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Taxonomic_assignment_Runner(self.settings).process()



def main():

    description = """

Version %s


MetaPlastHunter -


 The efficient and accurate plastid reads classification pipeline.

Quantitative aproach for eukaryotic metagenomics.

Available workflows:

[-reconstruction/--R] Reconstruction of large plastid parts

[-download_data/--fetch] Downloading data and directories preparation

[--taxonomic_classification/--C] Classification and visualization

[--rapid_classification, --Acc] Use it to lunch pipeline with in house kmer based preliminary classification

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

    parser.add_argument('--reconstruction','-R',action='store_true')
    parser.add_argument('--taxonomic_classification','-C',action='store_true')
    parser.add_argument('--download_data','-F',action='store_true')
    parser.add_argument('--rapid_classification', '-Acc',action='store_true')
    parser.add_argument('--settings','-S', metavar='settings', type=str)
    parser.add_argument('--threads',"-T",nargs='?', type=int,default=mp.cpu_count())

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()


    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

    process = Run(args.settings,args.threads)

    if args.taxonomic_classification:

        process.taxonomic_assigment()

    elif args.reconstruction:

        process.genomic_reconstruction()


    elif args.download_data:

        process.fetch_wf()

    elif args.rapid_classification:

        process.accelerated_taxonomic_assignment()
        
    else:

        logger.error("      Please specify pipeline")
        sys.exit()

    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)



if __name__ == "__main__":
    main()
