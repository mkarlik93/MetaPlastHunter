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
from subprocess import Popen, PIPE
import logging

#TODO

logger = logging.getLogger('src.bbmap_wrapper')
logging.basicConfig(level=logging.INFO)
#TODO

#Tworzenie log-a bbmapowego
#Moze lepiej przerobic na nazwe pliku
#Typu: -1 -2 -o jak w spades


class BBmap:

    """ This class is bbmap wrapper for mapping, filtering of metagenomic reads.
        Log goes to seperate file called sample_id+log.txt


     """

    def __init__(self,settings):
        " Below initial parameters taken from settings file "
        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbmap.sh"]
        self.db = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["bbmap_base"]
        self.db_silva = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["silva"]
        self.db_kmers = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["kmers"]
        self.minkmerhits = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["minkmerhits"]
        self.kmer_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["kmer_len"]
        self.remap_min_identity = Settings_loader_yaml(path=self.settings).yaml_handler()["Params"]["min_identity"]

    def filtering_conserved_regions(self,sample):
        " Filters out 16S rDNA and 18S rDNA "

        command="%sbbmap.sh fast=t reads=-1 in1=%s_de_complex_R1.fastq in2=%s_de_complex_R2.fastq nodisk ref=%s outu1=%s_filtered_chloroplasts_reads_R1.fq outu2=%s_filtered_chloroplasts_reads_R2.fq ambiguous=best outm1=%s_filtered_mapped_conserved_1.fq outm2=%s_filtered_mapped_conserved_2.fq  scafstats=filter_ribosomal.stats" % (self.path, sample, sample, self.db_silva, sample, sample, sample,sample)
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr


    def primary_mapping(self,sample):
        " Maps filtered reads to the reference chloroplast database "

        command="%sbbmap.sh fast=t nodisk reads=-1 idtag=t in1=%s_filtered_chloroplasts_reads_R1.fq in2=%s_filtered_chloroplasts_reads_R2.fq ref=%s scafstats=%s_chloroplasts.hitstats out=%s_final_mapped.sam bincov=bincov.txt minidentity=0.70 covbinsize=101" % (self.path, sample, sample, self.db,sample,sample)
        command = command.split(" ")
        logger.info("     Running primary mapping")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr

    def secondary_mapping(self,sample):
        " Maps reads again to the smaller database "

        command="%sbbmap.sh nodisk minidentity=%s idtag=t in1=%s_filtered_chloroplasts_reads_R1.fq in2=%s_filtered_chloroplasts_reads_R2.fq ref=tmp_ref_base.fasta outm1=%s_chloroplasts_reads_R1.fq outm2=%s_final_chloroplasts_reads_R2.fq ambiguous=all scafstats=%s_final_chloroplasts.hitstats statsfile=%s_final_mapping_stats.txt out=%s_final_mapped.sam bincov=bincov.txt covbinsize=200" % (self.path,self.remap_min_identity ,sample, sample, sample,sample,sample,sample,sample)
        command = command.split(" ")
        logger.info("     Running secondary mapping")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr



#This comes first -> kraken_out
class BBduk:
    """
    This class is a wrapper of bbduk programme from bbtools package.
    It contains:
        a) Function for filtering out low complex metagenomic reads
        b) Module for in-house made preliminary classification
    """
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

    def filtering(self,sample):

        """ Filters out low complexity reads """

        path = self.path
        logger.info("     Running  filtering" )
        command="%sbbduk.sh in1=%s_classif_R1.fastq in2=%s_classif_R2.fastq out1=%s_de_complex_R1.fastq out2=%s_de_complex_R2.fastq  outm=%s_repeat_regions_R1.fq outm2=%s_repeat_regions_R2.fq entropy=0.8 overwrite=true" % (path, sample, sample, sample, sample, sample, sample)
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return "DECOMPLEXATION"+"\n"+stderr

    def bbduk_pre_classification(self,sample):

        """

        Module for very fast preliminary classification. It needs specially prepared by kcompress
        from bbtools package

         """

        command="%sbbduk.sh  in1=%s_1.fastq in2=%s_2.fastq outm1=%s_classif_R1.fastq  outm2=%s_classif_R2.fastq minkmerhits=%s k=%s ref=%s" % (self.path,sample, sample,sample,sample,self.db_kmers)
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        logger.info("     Running  preliminary classification")
        return stderr


def log_writing(list,sample_id):

    with open(sample_id+"_log.txt","w") as f:
        for rec in list:
            f.write(rec+"\n")


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

            bblog = []
            logger.info("Procesing "+sra_id)
            dir = "%s/%s/" % (self.station_name,sra_id)
            os.chdir(dir)
            bblog.append(bbduk.filtering(sra_id))
            bblog.append(bbmap.filtering_conserved_regions(sra_id))
            bblog.append(bbmap.primary_mapping(sra_id))
            Coverage('bincov.txt',sra_id+"_chloroplasts.hitstats",self.settings).ref_for_remapping()
            bblog.append(bbmap.secondary_mapping(sra_id))
            log_writing(bblog,sra_id)
            os.chdir(starting_dir)

class BBpipe_with_bbduk_preliminary:

    def __init__(self, list_sra, station_name,settings):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings

        for sra_id in list_sra_ids:

            bblog = []
            logger.info("Procesing "+sra_id)
            dir = "%s/%s/" % (self.station_name,sra_id)
            os.chdir(dir)
            bblog.append(bbduk.bbduk_pre_classification(sra_id))
            bblog.append(bbduk.filtering(sra_id))
            bblog.append(bbmap.filtering_conserved_regions(sra_id))
            bblog.append(bbmap.primary_mapping(sra_id))
            Coverage('bincov.txt',sra_id+"_chloroplasts.hitstats",self.settings).ref_for_remapping()
            bblog.append(bbmap.secondary_mapping(sra_id))
            log_writing(bblog,sra_id)
#            Coverage("bincov_2.txt",sra_id+"_final_chloroplasts.hitstats",self.settings).report_cov()
            os.chdir(starting_dir)


class Single_mapping:

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

            bblog = []
            logger.info("Procesing "+sra_id)
            dir = "%s/%s/" % (self.station_name,sra_id)
            os.chdir(dir)
            bblog.append(bbduk.filtering(sra_id))
            bblog.append(bbmap.filtering_conserved_regions(sra_id))
            bblog.append(bbmap.primary_mapping(sra_id))
            log_writing(bblog,sra_id)
            os.chdir(starting_dir)
