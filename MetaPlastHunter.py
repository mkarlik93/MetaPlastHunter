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

logger = logging.getLogger("MetaPlastHunter")
logging.basicConfig(level=logging.INFO)


class Run:

    """
    Main class for running MetaPlastHunter

    """

    def __init__(self,list_sra, station_name,settings,threads):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings
        self.threads = threads

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
        Run_analysis_sam_lca(self.list_sra, self.station_name,self.settings).process()

    def fetch_wf(self):

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
        Run_analysis_sam_lca(self.list_sra, self.station_name,self.settings).process()

    def recalculation_wf(self):

        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        BBpipe(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Run_analysis_sam_lca(self.list_sra, self.station_name,self.settings).process()


    def classification_with_kmer_method(self):
        logger.info( "     [%s] Testing settings file " % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Settings_loader_yaml(self.settings).yaml_check_settings_file()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        BBpipe_with_bbduk_preliminary(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Taxonomic assignment based on SAM file" % (strftime("%a, %d %b %Y %H:%M:%S +2", gmtime())))
        Run_analysis_sam_lca(self.list_sra, self.station_name,self.settings).process()




if __name__ == "__main__":

    from time import gmtime, strftime
    import sys
    import os
    import argparse
    from src.bbmap_wrapper import BBpipe
    from src.bbmap_wrapper import BBpipe_with_bbduk_preliminary
    from src.krakenize import Pipeline_kraken
    from src.sam_analyzer import Run_analysis_sam_lca
    from src.get_data import Pipeline_fetch
    from src.settings import Settings_loader_yaml
    import multiprocessing as mp



    description = """

Version %s


MetaPlastHunter - The efficient and accurate plastid reads classification pipeline.


Quantitative aproach for eukaryotic metagenomics.


Available workflows:

[-full_wf/--full] From downloading data to classification and visualization

[-download_data/--fetch] Downloading data and directories preparation

[-classification_wf/--classify] Classification and visualization

[-recalculation_wf/--recalculate] Use it for recalculate taxonomic assignemnt based on LCA algorithm




Obligatory arguments:

sra_ids

station_name

settings


Facultative arguments:

threads



If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com


This sofware was written by %s.
""" % (__version__,__author__)

    epilog = """




"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)


    parser.add_argument('-full_wf','--full',action='store_true')
    parser.add_argument('-classification_wf','--classify',action='store_true')
    parser.add_argument('-recalculation_wf','--recalculate',action='store_true')
    parser.add_argument('-download_data','--fetch',action='store_true')
    parser.add_argument('-process_with_kmer_classif','--kmer_classif',action='store_true')
    parser.add_argument('sra_ids', metavar='sra_ids', type=str)
    parser.add_argument('station_name', metavar='station_name', type=str)
    parser.add_argument('settings', metavar='settings', type=str)
    parser.add_argument('threads',nargs='?', type=int,default=mp.cpu_count())

#   parser.add_argument('number', nargs='?', type=int,default=None)
#   parser.add_argument('format', nargs='?', type=str,default="fasta")



    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

    process = Run(args.sra_ids, args.station_name,args.settings,args.threads)

    if args.full:

        process.full_wf()

    elif args.classify:

        process.classification_wf()

    elif args.recalculate:

        process.recalculation_wf()

    elif args.fetch:

        process.fetch_wf()

    elif arg.kmer_classif:

        process.classification_with_kmer_method()

    else:

        logger.error("      Please specify pipeline")
        sys.exit()

    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)
