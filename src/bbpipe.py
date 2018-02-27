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

#Sra_id /name super istotny

from settings import *
import os
import sys

class BBmap:

    def __init__(self,settings):

        if settings == None:
            self.checkForBBmap()
            self.path = ""
            self.db = Settings_loader(mode="ref_base",path=self.path).read_path()["ref_base"]
        else:

            self.settings = settings
            self.path = Settings_loader(mode="bbmap.sh",path=self.settings).read_path()["bbmap.sh"]
            self.db = Settings_loader(mode="bbmap.sh",path=self.settings).read_database()["bbmap_base"]


    def checkForBBmap(self):
        """Check to see if Kraken is on the system before we try to run it."""

        try:
            subprocess.call(['bbmap.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "     [Error] Make sure BBmap is on your system path or set usage to path in settings.txt"
            sys.exit()


    def run(self,sample):
        path = self.path
        db = self.db
        command="%sbbmap.sh fast=t minidentity=0.95 reads=-1  in1=%s_de_complex_R1.fastq in2=%s_de_complex_R2.fastq ref=%s outm1=%s_chloroplasts_reads_R1.fq outm2=%s_chloroplasts_reads_R2.fq ambiguous=best  scafstats=%s_chloroplasts.hitstats out=%s_mapped.sam" % (path, sample, sample, db, sample, sample, sample, sample)
        print "     Running command: [%s]" % command
        os.system(command)



#This comes first -> kraken_out
class BBduk:

    def __init__(self,settings):

        if settings == False:
            self.checkForBBduk()
            self.path = ""
        else:
            self.settings = settings
            self.path = Settings_loader(mode="bbduk.sh",path=self.settings).read_path()["bbduk.sh"]

    def checkForBBduk(self):
        """Check to see if BBduk is on the system before we try to run it."""

        try:
            subprocess.call(['bbduk.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "  [Error] Make sure BBduk is on your system path or set usage to path in settings.txt"
            sys.exit()

    def run(self, sample):
        path = self.path
        command="%sbbduk.sh in1=%s_classif_R1.fastq in2=%s_classif_R2.fastq out1=%s_de_complex_R1.fastq out2=%s_de_complex_R2.fastq  outm=%s_repeat_regions_R1.fq outm2=%s_repeat_regions_R2.fq entropy=0.8 overwrite=true" % (path, sample, sample, sample, sample, sample, sample)
        print "     Running command: [%s]" % command
        os.system(command)

class BBpipe:
    def __init__(self, list_sra, station_name,settings):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings

    def process(self):
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        starting_dir = os.getcwd()
        bbduk = BBduk(self.settings)
        bbmap = BBmap(self.settings)
        for i in list_sra_ids:
            print "Procesing "+i
            dir = "%s/%s/" % (self.station_name,i)
            os.chdir(dir)
            bbduk.run(i)
            bbmap.run(i)
            os.chdir(starting_dir)
