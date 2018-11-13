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


#To moze zostac przeniesione gdzies indziej - inny typ.

from cov import Coverage
from settings import *
import os
import sys
from subprocess import Popen, PIPE
import logging

#TODO

logger = logging.getLogger('src.genome_reconstruction')
logging.basicConfig(level=logging.INFO)

class Genome_reconstruction:


    """

    This class is bbmap wrapper for mapping, filtering of metagenomic reads.
        Log goes to seperate file called sample_id+log.txt

     """

    def __init__(self,settings):

        """ Below initial parameters taken from settings file

        Params are loaded from settings file

        """

        self.settings = settings
        self.path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["bbmap.sh"]

        self.reference = Settings_loader_yaml(path=self.settings).yaml_handler()["Databases and mapping files"]["amino_sketch_database"]

        self.spades_path = Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["spades"]
        if self.spades_path == None:
            self.spades_path = ""

        self.seqtk_path =  Settings_loader_yaml(path=self.settings).yaml_handler()["Software dependencies"]["seqtk"]
        if self.seqtk_path == None:
            self.seqtk_path = ""

        self.project_name = Settings_loader_yaml(path=self.settings).yaml_handler()["Project name"]["name"]
        self.mode = Settings_loader_yaml(path=self.settings).yaml_handler()["Mode type single or paired"]["mode"]

        self.contigs = Settings_loader_yaml(path=self.settings).yaml_handler()["Input contigs"]["contigs"]
        self.orginal_reads_1 =  Settings_loader_yaml(path=self.settings).yaml_handler()["Input reads"]["read_1"]
        self.orginal_reads_2 =  Settings_loader_yaml(path=self.settings).yaml_handler()["Input reads"]["read_2"]


        self.minhits = Settings_loader_yaml(path=self.settings).yaml_handler()["Sketch comparsion"]["min_hits"]
        self.k_len = Settings_loader_yaml(path=self.settings).yaml_handler()["Sketch comparsion"]["k_len"]
        self.project = self.project_name

        if len(self.orginal_reads_2) == 0:
            self.mode == "single"
        else:
            self.mode == "paired"

    def sketch_make_aa(self):

        " Make amino acid sketch from input file "

        command="%sbbsketch.sh in=%s mode=sequence k=12 amino=t out=in_sketch" % (self.path,self.contigs)
        command = command.split(" ")
        process = Popen(command, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        return stderr


    def compare_sketch(self):

        " Sketch comparing  "

        command= """%scomparesketch.sh in=in_sketch ref=%s records=1 minhits=5 amino=t k=12  | grep "WKID" -1 | grep 'Query' | awk '{split($0,a,"\t"); print a[1]}' | awk '{split($0,a,":"); print a[2]}'  | sed 's/^.//' > target_list """ % (self.path, self.reference)
        logger.info("     Sketch comparsion")


        process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)

        stdout, stderr = process.communicate()

        return stderr

    def seqtk(self):

        " Seqtk extraction"

        command="%sseqtk subseq %s target_list > db_for_remapping" % (self.seqtk_path,self.contigs)
        logger.info("     Contig extraction")
        process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
        stdout, stderr = process.communicate()
        return stderr


    def read_extraction(self):

        " Mapping filtered reads to the reference chloroplast database "

        if self.mode == "paired":
            command="%sbbmap.sh  nodisk in1=%s in2=%s ref=db_for_remapping outm1=%s_extracted_chloroplasts_reads_R1.fq outm2=%s_extracted_chloroplasts_reads_R2.fq " % (self.path,self.orginal_reads_1,self.orginal_reads_2,self.project_name,self.project_name)
        if self.mode == "single":
            command="%sbbmap.sh  nodisk in=%s ref=db_for_remapping outm1=%s_extracted_chloroplasts_reads_R1.fq " % (self.path,self.orginal_reads_1,self.project_name)

        logger.info("     Running  read extraction")
        process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
        stdout, stderr = process.communicate()
        return stderr

    def read_correction_elongation(self):

        " Based on k-mer counts, extend by 250 bp "


        if self.mode == "paired":
            command="%stadpole.sh mode=extend ecc=t tossjunk=t extendright=250 extendleft=250 in1=%s_chloroplasts_reads_R1.fq in2=%s_chloroplasts_reads_R2.fq  out1=%s_tmp_R1.fq out2=%s_tmp_R2.fq " % (self.path,self.project_name,self.project_name,self.project_name,self.project_name)
        if self.mode == "single":
            command="%stadpole.sh mode=correct ecc=t tossjunk=t extendright=250 extendleft=250 in=%s_chloroplasts_reads_R1.fq  outm1=%s_tmp_R1.fq " % (self.path,self.project_name,self.project_name)

        logger.info("     Running assembly by Spades")
        process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
        stdout, stderr = process.communicate()
        return stderr


    def final_assembly(self):

        if self.mode == "paired":
            command="%sspades --only-assembler -1 %s_tmp_R1.fq  -2 %s_tmp_R2.fq  -o %s_spades_assembly " % (self.spades_path,self.project_name,self.project_name,self.project_name)
        if self.mode == "single":
            command="%sspades --only-assembler -s %s_tmp_R1.fq -o %s_spades_assembly " % (self.spades_path,self.project_name,self.project_name)

        logger.info("     Running correction")
        process = Popen(command, stdout=PIPE, stderr=PIPE,shell=True)
        stdout, stderr = process.communicate()
        print " Thank you for usage genome reconstruction pipeline!! Check your results! "
        return stderr

def log_writing(list,sample_id):

    with open(sample_id+"_log.txt","w") as f:
        for rec in list:
            f.write(rec+"\n")

class Genome_reconstruction_pipe:

    def __init__(self, settings):

        self.settings = settings

    def process(self):


        starting_dir = os.getcwd()
        gr = Genome_reconstruction(self.settings)

        project_name = gr.project

        grlog = []

        logger.info("Procesing "+project_name)
        dir = "%s/" % (project_name)

        if os.path.isdir(dir):
            os.chdir(dir)
        else:
            os.mkdir(dir)
            os.chdir(dir)
        grlog = []
        print "Happy Hunting!"
        grlog.append(gr.sketch_make_aa())
        grlog.append(gr.compare_sketch())
        grlog.append(gr.seqtk())
        grlog.append(gr.read_extraction())
        grlog.append(gr.read_correction_elongation())
        grlog.append(gr.final_assembly())
        log_writing(grlog,str(project_name))

        os.remove(project_name+'_tmp_R1.fq')
        os.remove(project_name+'_tmp_R2.fq')
