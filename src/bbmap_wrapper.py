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

from cov import Coverage
from settings import *
import os
import sys
import logging

logger = logging.getLogger('src.bbmap_wrapper')
logging.basicConfig(level=logging.INFO)


class BBmap:

    def __init__(self,settings):

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbmap.sh"]
        self.db = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["bbmap_base"]
        self.db_silva = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["silva"]


#    def checkForBBmap(self):
#
#        """Checks that Kraken is on the system before we try to run it."""
#
#        try:
#            subprocess.call(['bbmap.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
#        except:
#            logger.error("   Make sure BBmap is on your system path or set usage to path in settings.txt")
#            sys.exit()

#For filtering out 16sRNA 'ALL' can be just bacterial and archean
    def ssu_n_lsu_rDNA_flitering_run(self,sample):
        command="%sbbmap.sh fast=t minidentity=0.70 reads=-1 in1=%s_de_complex_R1.fastq in2=%s_de_complex_R2.fastq nodisk ref=%s outu1=%s_filtered_chloroplasts_reads_R1.fq outu2=%s_filtered_chloroplasts_reads_R2.fq ambiguous=best out=%s_filtered_mapped.sam" % (self.path, sample, sample, self.db_silva, sample, sample, sample)
        logger.info("     Running command: [%s]" % command)
        os.system(command)


#remapping with smaller database
#    def remap_run(self,sample):
#        command="%sbbmap.sh nodisk minidentity=0.70 idtag=t in1=%s_final_chloroplasts_reads_R1.fq in2=%s_final_chloroplasts_reads_R2.fq ref=tmp_ref_base.fasta outu1=%s_f_chloroplasts_reads_R1.fq outu2=%s_f_chloroplasts_reads_R2.fq ambiguous=all scafstats=%s_final_chloroplasts.hitstats statsfile=%s_final_mapping_stats.txt out=%s_final_mapped.sam bincov=bincov.txt covbinsize=200" % (self.path, sample, sample, sample, sample,sample,sample,sample)
#        logger.info("     Running command: [%s]" % command)
#        os.system(command)

    def remap_run(self,sample):
        command="%sbbmap.sh nodisk minidentity=0.70 idtag=t in1=%s_filtered_chloroplasts_reads_R1.fq in2=%s_filtered_chloroplasts_reads_R2.fq ref=tmp_ref_base.fasta outu1=%s_f_chloroplasts_reads_R1.fq outu2=%s_f_chloroplasts_reads_R2.fq ambiguous=all scafstats=%s_final_chloroplasts.hitstats statsfile=%s_final_mapping_stats.txt out=%s_final_mapped.sam bincov=bincov.txt covbinsize=200" % (self.path, sample, sample, sample, sample,sample,sample,sample)
        logger.info("     Running command: [%s]" % command)
        os.system(command)

    def run(self,sample):

        command="%sbbmap.sh fast=t minidentity=0.70 nodisk reads=-1 idtag=t in1=%s_filtered_chloroplasts_reads_R1.fq in2=%s_filtered_chloroplasts_reads_R2.fq ref=%s scafstats=%s_chloroplasts.hitstats out=%s_mapped.sam bincov=bincov.txt covbinsize=200" % (self.path, sample, sample, self.db,sample, sample)
        logger.info("      Running command: [%s]" % command)
        os.system(command)

#This comes first -> kraken_out
class BBduk:

    def __init__(self,settings):

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbduk.sh"]

        if self.path == "":
            self.checkForBBduk()


    def checkForBBduk(self):

        """Checks to see if BBduk is on the system's path before we try to run it."""

        try:
            subprocess.call(['bbduk.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logger.error("Make sure BBduk is on your system path or set proper path in settings.txt")
            sys.exit()

    def run(self, sample):
        path = self.path
        command="%sbbduk.sh in1=%s_classif_R1.fastq in2=%s_classif_R2.fastq out1=%s_de_complex_R1.fastq out2=%s_de_complex_R2.fastq  outm=%s_repeat_regions_R1.fq outm2=%s_repeat_regions_R2.fq entropy=0.8 overwrite=true" % (path, sample, sample, sample, sample, sample, sample)
        logger.info("     Running command: [%s]" % command)
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

        for sra_id in list_sra_ids:
            logger.info("Procesing "+sra_id)
            dir = "%s/%s/" % (self.station_name,sra_id)
            os.chdir(dir)
            bbduk.run(sra_id)
            bbmap.ssu_n_lsu_rDNA_flitering_run(sra_id)
            bbmap.run(sra_id)
            Coverage('bincov.txt',sra_id+"_chloroplasts.hitstats",self.settings).ref_for_remapping()
            bbmap.remap_run(sra_id)
            os.chdir(starting_dir)
