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

logger = logging.getLogger('src.bbmap_wrapper')
logging.basicConfig(level=logging.INFO)

class BBmap:

    """ This class is bbmap wrapper for mapping, filtering of metagenomic reads.
        Log goes to seperate file called sample_id+log.txt

     """

    def __init__(self,input, input2, output, settings,threads):

        """
        Input Parameters
        ----------

        settings : str
            The path to general settings file which keeps hyperparameters

        Calculated gobal variables
        ----------
        _seqid: dictionary

        path: str
            Sample name

        db: dictionary

        db_silva: dictionary

        minkmerhits: dictionary

        kmer_len: graph data structure

            Orginal data structure (weighted digraph) which keeps NCBI tree with proposed taxonomic positions (LTU and TTU)

        remap_min_identity: graph data stucture

            Prunned _lca_graph that keeps only taxonomic positions that have been above the threshold (min support value)

        reads_1: fastq file (or SAM file)

        reads_2: fastq file

        project_name: name of the output

        """

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbmap.sh"]
        self.db = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["bbmap_base"]
        self.db_silva = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["silva"]
        self.db_kmers = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["kmers"]
        self.minkmerhits = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["minkmerhits"]
        self.kmer_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["kmer_len"]
        self.bincov_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Params"]["bincov_len"]
        self.remap_min_identity = Settings_loader_yaml(path=self.settings).yaml_handler()["Params"]["min_identity"]
        self.project_name = output
        self.reads_1 =  input
        self.reads_2 =  input2
        self.threads = str(threads)

    def filtering_conserved_regions(self):

        """ Filters out parts of ribosomal operons based on mapping strategy """

        command="%sbbmap.sh fast=t reads=-1 in1=%s_de_complex_R1.fastq in2=%s_de_complex_R2.fastq ref=%s outu1=%s_filtered_chloroplasts_reads_R1.fq outu2=%s_filtered_chloroplasts_reads_R2.fq ambiguous=best outm1=%s_filtered_mapped_conserved_1.fq outm2=%s_filtered_mapped_conserved_2.fq  scafstats=filter_ribosomal.stats threads=%s" % (self.path, self.project_name, self.project_name, self.db_silva, self.project_name, self.project_name, self.project_name,self.project_name,self.threads)
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr

    def primary_mapping(self):

        """ Maps filtered reads to the reference chloroplast database """

        command="%sbbmap.sh ambiguous=best minidentity=%s nodisk reads=-1 idtag=t in1=%s_filtered_chloroplasts_reads_R1.fq in2=%s_filtered_chloroplasts_reads_R2.fq ref=%s scafstats=%s_chloroplasts.hitstats out=%s_final_mapped.sam bincov=bincov.txt outm1=%s_chloroplasts_reads_R1.fq outm2=%s_final_chloroplasts_reads_R2.fq covbinsize=%s threads=%s" % (self.path,self.remap_min_identity ,self.project_name, self.project_name, self.db,self.project_name,self.project_name,self.project_name,self.project_name, self.bincov_len,self.threads)
        command = command.split(" ")
        logger.info("     Running primary mapping")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr

class BBduk:

    """
    This class is a wrapper of bbduk programme from bbtools package.
    It contains:
        a) Function for filtering out low complex metagenomic reads
    """

    def __init__(self,input,input2,output,settings,threads):

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbduk.sh"]
        self.project_name = output
        self.reads_1 =  input
        self.reads_2 =  input2
        self.db_kmers = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["kmers"]
        self.minkmerhits = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["minkmerhits"]
        self.kmer_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Preliminary classification"]["kmer_len"]
        self.silva_compressed = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["compressed_Silva"]
        self.threads = str(threads)

        if self.path == "":
            self.checkForBBduk()

    def checkForBBduk(self):

        """Checks to see if BBduk is on the system's path before we try to run it."""

        try:
            subprocess.call(['bbduk.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

        except:
            logger.warning("Make sure BBduk is on your system path or set proper path in settings.txt")
            sys.exit()

    def filtering_conserved_regions_based_on_kmers(self):

        """ Filters out parts of ribosomal operons based on exact kmer matching strategy """

        command = "%sbbduk.sh minkmerhits=%s k=%s in1=%s_de_complex_R1.fastq in2=%s_de_complex_R2.fastq ref=%s outu1=%s_filtered_chloroplasts_reads_R1.fq outu2=%s_filtered_chloroplasts_reads_R2.fq  outm1=%s_filtered_mapped_conserved_1.fq outm2=%s_filtered_mapped_conserved_2.fq threads=%s" % (self.path,self.minkmerhits,self.kmer_len,self.project_name, self.project_name, self.silva_compressed, self.project_name, self.project_name, self.project_name,self.project_name,self.threads)
        splited_cmd = command.split(" ")

        process = Popen(splited_cmd,stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()

        return "Filtration conserved reads using kmer catching method"+"\n"+stderr

    def filtering_without_pre_classif(self):

        """ Filters out low complexity reads """

        path = self.path
        logger.info("     Running  filtering" )
        command="%sbbduk.sh in1=%s in2=%s out1=%s_de_complex_R1.fastq out2=%s_de_complex_R2.fastq  outm=%s_repeat_regions_R1.fq outm2=%s_repeat_regions_R2.fq entropy=0.8 overwrite=true threads=%s" % (path, self.reads_1, self.reads_2, self.project_name, self.project_name, self.project_name, self.project_name,str(self.threads))
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return "DECOMPLEXATION"+"\n"+stderr

    def filtering_with_pre_classif(self):

        """ Filters out low complexity reads """

        path = self.path
        logger.info("     Running  filtering" )
        command="%sbbduk.sh in1=%s_classif_R1.fastq in2=%s_classif_R2.fastq out1=%s_de_complex_R1.fastq out2=%s_de_complex_R2.fastq  outm=%s_repeat_regions_R1.fq outm2=%s_repeat_regions_R2.fq entropy=0.8 overwrite=true threads=%s" % (path, self.project_name, self.project_name, self.project_name, self.project_name, self.project_name, self.project_name,str(self.threads))
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return "DECOMPLEXATION"+"\n"+stderr

    def bbduk_pre_classification(self):

        """

        Module for speeding up preliminary classification. It needs file specially prepared by kcompress
        from bbtools package

        """
        path = self.path
        logger.info("     Running  pre classification" )
        command="%sbbduk.sh  in1=%s in2=%s outm1=%s_classif_R1.fastq  outm2=%s_classif_R2.fastq minkmerhits=%s k=%s ref=%s threads=%s" % (self.path,self.reads_1, self.reads_2,self.project_name,self.project_name,self.minkmerhits,self.kmer_len,self.db_kmers,self.threads)
        command = command.split(" ")
        process = Popen(command,shell=False,stdout=PIPE)
        stdout, stderr = process.communicate()
        return "Pre classification"+"\n"+stdout

def log_writing(list,sample_id):

    with open(sample_id+"_log.txt","w") as f:
        for rec in list:
            f.write(rec+"\n")

class Pileup:

    """
    This class is a wrapper of pileu programme from bbtools package.
    It contains:
        a) Function for generating bincov.txt from SAM file which is nessecary in taxonomic classification

    """

    def __init__(self,input,settings):

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["pileup.sh"]
        self.bincov_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Params"]["bincov_len"]
        self.input = input


        if self.path == "":

            self.checkForPileup()

    def checkForPileup(self):

        """Checks to see if BBduk is on the system's path before we try to run it."""

        try:
            subprocess.call(['pileup.sh', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

        except:
            logger.waring("Make sure pileup is on your system path or set proper path in settings.txt")
            sys.exit()

    def prepare_cov_file(self):

        command= "%spileup.sh in=%s bincov=bincov.txt binsize=%s" % (self.path, self.input, str(self.bincov_len))
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr

class Mapping_runner:

    def __init__(self,input,input2,output, settings,threads,mapping):

        self.settings = settings
        self.input  =  input
        self.input2 = input2
        self.output = output
        self.threads = threads
        self.mapping = mapping

    def process(self):

        starting_dir = os.getcwd()

        bbduk = BBduk(self.input,self.input2,self.output,self.settings,self.threads)
        bbmap = BBmap(self.input,self.input2,self.output,self.settings,self.threads)

        project_name = self.output

        grlog = []

        logger.info("Procesing "+project_name)
        dir = "%s/" % (project_name)

        if os.path.isdir(dir):

            logger.warning("This directory exists! Please change the output name!")
            sys.exit()

        else:
            os.mkdir(dir)
            os.chdir(dir)
        bblog = []
        logger.info("Procesing "+project_name)

        if self.mapping is True:

            bblog.append(bbduk.filtering_without_pre_classif())
            bblog.append(bbmap.filtering_conserved_regions())
            bblog.append(bbmap.primary_mapping())
            log_writing(bblog,project_name)
            os.chdir(starting_dir)

        else:

            bblog.append(bbduk.filtering_without_pre_classif())
            bblog.append(bbduk.filtering_conserved_regions_based_on_kmers())
            bblog.append(bbmap.primary_mapping())
            log_writing(bblog,project_name)
            os.chdir(starting_dir)


class RapidRunner:


    def __init__(self,input,input2,output,settings):

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbduk.sh"]
        self.settings = settings
        self.input  =  input
        self.input2 = input2
        self.output = output

    def process(self):

        starting_dir = os.getcwd()

        bbduk = BBduk(self.input,self.input2,self.output,self.settings)
        bbmap = BBmap(self.input,self.input2,self.output,self.settings)

        project_name = self.output

        grlog = []

        logger.info("Procesing "+project_name)
        dir = "%s/" % (project_name)

        if os.path.isdir(dir):

            logger.warning("This directory name exists! Please change the output name!")
            sys.exit()

        else:
            os.mkdir(dir)
            os.chdir(dir)

        bblog = []
        logger.info("Procesing "+project_name)
        bbduk.bbduk_pre_classification()

        bblog.append(bbduk.filtering_with_pre_classif())
        bblog.append(bbmap.filtering_conserved_regions())

        bblog.append(bbmap.primary_mapping())
        log_writing(bblog,project_name)
        os.chdir(starting_dir)

class SAM2coverage:

    def __init__(self,input,output, settings):

        self.settings = settings
        self.input  = input
        self.output = output

    def process(self):

        starting_dir = os.getcwd()

        pileup = Pileup(self.input,self.settings)
        project_name = self.output
        grlog = []
        logger.info("Procesing "+project_name)
        dir = "%s/" % (project_name)

        if os.path.isdir(dir):

            logger.warning("This directory name exists! Please change the output name!")
            sys.exit()

        else:

            os.mkdir(dir)
            os.chdir(dir)

        bblog = []
        logger.info("Procesing "+project_name)
        bblog.append(pileup.prepare_cov_file())
        log_writing(bblog,project_name)
        os.chdir(starting_dir)
