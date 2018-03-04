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



#TODO
#For now this script doesn't work correctly
#Take script named cov.py from the server and modify use modules from kraken_output_analysis
from Bio import SeqIO
import glob
import numpy as np
import sys


class Coverage:

    def __init__(self,filename, histstats):

        self.loaded = self.load_bincov(filename,histstats)
        self.filename = filename
        self.histstats = histstats

    def load_genomes_len(self):
        genomes_dict = {}
        with open(filename,"r") as f:
            for i in f:
                splited = i.split("\t")
                genomes_dict[splited[0]] = splited[1].strip("\n")
                unique_names.add(splited[0])
        return genomes_dict


    def load_bincov(self,filename,hist):

        with open(filename,"r") as f:
            organisms = {}
            for i in f:
                splited = i.split("\t")
                if splited[0] not in organisms and splited[0][0] != "#":
                    organisms[splited[0]] = [float(splited[1])]
                else:
                    if splited[0][0] != "#":
                        organisms[splited[0]].append(float(splited[1]))
        return organisms


    def report_cov(self):
        threshold = 1
        dict_of_genomes = {}
        dict = self.loaded
        for i in dict:
            name = i
            record = dict[i]
            mean =  np.mean(record)
            if mean > 0:
                all_seq = len(record)
                list_tmp = [i for i,v in enumerate(record) if v > threshold]
                covered_part = len(list_tmp)/float(all_seq)

                if covered_part > 0.5:
                    if i is not 0:
                        dict_of_genomes[name] = str((covered_part*100.0))
                else:
                    pass
        if len(dict_of_genomes) == 0:
            print "There is no fully/partially covered genomes"
            sys.exit()

        else:
            with open("covered_genomes.csv","w") as f:
                f.write("Genome name and id,coverage[%]")
                for key in dict_of_genomes:
                    f.write("%s, %s" % (key, dict_of_genomes[key]))
                    return "Report was generated"
