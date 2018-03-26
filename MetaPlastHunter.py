#!/usr/bin/env python
###############################################################################
#                                                                             #
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





#LET's write whole main



class WholePipeline:

    def __init__(self,list_sra, station_name,settings,threads):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings
        self.threads = threads

    def run(self):


        logger.info( "     [%s] Dowloading data" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_fetch(self.list_sra,self.station_name,self.settings).run()
        logger.info("     [%s] Kraken classification" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_kraken(self.list_sra, self.station_name,self.settings,self.threads).run()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        BBpipe(list_sra,station_name,settings).process()
        logger.info("     [%s] Output analysis" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Run_analysis(list_sra, station_name,4).process()

class Pipeline_without_downloading:

    def __init__(self,list_sra, station_name,settings,threads):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings
        self.threads = threads

    def run(self):

        logger.info("     [%s] Kraken classification" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Pipeline_kraken(self.list_sra, self.station_name,self.settings,self.threads).run()
        logger.info("     [%s] BBtools postprocessing" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        BBpipe(self.list_sra,self.station_name,self.settings).process()
        logger.info("     [%s] Output analysis" % (strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
        Run_analysis(self.list_sra, self.station_name,4,self.settings).process()


if __name__ == "__main__":

    from time import gmtime, strftime
    import sys
    import os
    import argparse
    from src.bbmap_wrapper import BBpipe
    from src.krakenize import Pipeline_kraken
    from src.sam_analyzer import Run_analysis_sam_lca
    from src.get_data import Pipeline_fetch
    import multiprocessing as mp

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)


    description = """

Version %s


The first version of MetaPlastHunter.

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


    parser.add_argument('-analysis_with_dumping_data','--full',action='store_true')
    parser.add_argument('-analysis','--partial',action='store_true')
    parser.add_argument('sra_ids', metavar='sra_ids', type=str)
    parser.add_argument('station_name', metavar='station_name', type=str)
    parser.add_argument('settings', metavar='settings', type=str)
    parser.add_argument('threads',nargs='?', type=int,default=mp.cpu_count())

    parser.add_argument('number', nargs='?', type=int,default=None)
    parser.add_argument('format', nargs='?', type=str,default="fasta")



    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    if args.full:
        WholePipeline(args.list_sra, args.station_name,args.settings,args.threads).run()
    if args.partial:
        Pipeline_without_downloading(args.sra_ids, args.station_name,args.settings,args.threads).run()
    else:
        logger.error("      Please specify pipeline")
        sys.exit()
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)
