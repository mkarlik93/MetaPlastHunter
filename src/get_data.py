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

from settings import *
from time import gmtime, strftime
import sys
import os
import argparse
from multiprocessing import Process
import glob
import logging

logger = logging.getLogger('src.get_data')
logging.basicConfig(level=logging.INFO)


def fastq_dump_sra_file(station_name,list_sra,path):
    command_create_dir = "mkdir %s/%s" % (station_name,list_sra)
    os.system(command_create_dir)
    command = "%sfastq-dump %s --skip-technical -I --split-3" % (path,list_sra)
    os.chdir("%s/%s/" % (station_name,list_sra))
    os.system(command)
    logger.info("Downloading of %s has been started" % list_sra)


class Pipeline_fetch:

    def __init__(self, list_sra, station_name,settings):

        self.list_sra = list_sra
        self.station_name = station_name
        self.path = Settings_loader(mode="fastq-dump",path=settings).read_path()["fastq-dump"]


    def create_station_dir(self):
        command  = "mkdir %s" % (self.station_name)
        os.system(command)
        logger.info("Created direcory named: %s" % (self.station_name))

    def preprocess_sra_id(self):
        'Splits sra_ids line into list of sra ids'
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        return list_sra_ids


    def multiprocess(self):
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        path = self.path
        for i in list_sra_ids:

                job = Process(target=fastq_dump_sra_file, args=(self.station_name,self.list_sra,self.path))
                job.start()
                job.join()

    #CHECK THIS
    def evaluation(self):
        cwd = os.getcwd()
        station_name = self.station_name
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        os.chdir("%s/%s/" % (cwd,station_name))
        for i in list_sra_ids:
            os.chdir( "%s/%s/%s/" % (cwd,station_name,i))
            if len(glob.glob("*.fastq")) == 2:
                logger.info("For "+str(i)+" everything is ok!")
            else:
                logger.info("check it : "+str(i))
            os.chdir("%s/%s/" % (cwd,station_name))
        os.chdir(cwd)

    def fastqc_report(self):
        cwd = os.getcwd()
        station_name = self.station_name
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        os.chdir("%s/%s/" % (cwd,station_name))
        for i in list_sra_ids:
            os.chdir( "%s/%s/%s/" % (cwd,station_name,i))
            if len(glob.glob("*.fastq")) == 2:
                for file in glob.glob("*.fastqc"):
                    command = "fastqc "+i
                    os.system(command)
                    logger.info("For "+str(i)+" everything is ok!")
            else:
                logger.info("check it : "+str(i))
            os.chdir("%s/%s/" % (cwd,station_name))
        os.chdir(cwd)

    def run(self):
        self.create_station_dir()
        self.multiprocess()
        self.evaluation()
#        self.fastqc_report()


if __name__ == "__main__":



    description = """

Version 1.00

Script was designed for getting data from SRA repository and evaluation.

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com
"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('sra_ids', metavar='-sra_ids', type=str)
    parser.add_argument('station_name', metavar='-station_name', type=str)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    run = Pipeline_fetch(args.sra_ids,args.station_name).run()
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)
