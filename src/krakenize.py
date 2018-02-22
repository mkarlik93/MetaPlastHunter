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


from kraken import *

class Pipeline_kraken:

    def __init__(self, list_sra, station_name,settings):
        self.list_sra = list_sra
        self.station_name = station_name
        self.database_dir = database_dir
        self.threads = threads
        self.settings = settings


    def kraken_not_multi(self):
        database_dir = self.database_dir
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        starting_dir = os.getcwd()
        kraken = KrakenRunner(self.threads,self.settings)
        for i in list_sra_ids:
            print "kraken for "+i
            read_name_1 = i+"_1.fastq"
            read_name_2 = i+"_2.fastq"
            dir = "%s/%s/" % (self.station_name,i)
            os.chdir(dir)
            kraken.run_classification(read_name_1, read_name_2,dir)
            kraken.run_report(dir)
#            command = "/opt/kraken/kraken -t %s --db %s --paired %s %s --out-fmt paired --fastq-output --classified-out classif > kraken_out" % (str(self.threads),database_dir,read_name_1,read_name_2)
#            os.system(command)

#            command_report = "/opt/kraken/kraken-report --db %s kraken_out > kraken_report_%s.txt" % (str(self.threads), self.station_name)
#            command_translate = "/opt/kraken/kraken-translate --db /home/karlicki/custom kraken_out > kraken_labels_%s.txt" % (self.station_name)
#            os.system(command_report)
#            os.system(command_translate)
            os.chdir(starting_dir)
            os.remove("%s/%s/%s" % (self.station_name,i,read_name_1))
            os.remove("%s/%s/%s" % (self.station_name,i,read_name_2))

    def run(self):
        self.kraken_not_multi()


    def multiprocess(self):
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        for i in list_sra_ids:
            proc = Process(target=fastq_dump_sra_file, args=(self.station_name,i))
            proc.start()
            print "Downloading has started"

#TO trzeba zmienic sa nowe zmienne
if __name__ == "__main__":

    from time import gmtime, strftime
    import sys
    import os
    import argparse
    from multiprocessing import Process


    description = """

Version 1.0

This script was designed to use Kraken for metagenomic sequence classification.

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com
"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)




    parser.add_argument('sra_ids', metavar='sra_ids', type=str)
    parser.add_argument('station_name', metavar='station_name', type=str)
    parser.add_argument('database_dir', metavar='database_dir', type=str)
    parser.add_argument('threads', metavar='threads', type=int)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    run = Pipeline_kraken(args.sra_ids,args.station_name).run()
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    print "Starting time: "+start
    print "Ending time: "+end
