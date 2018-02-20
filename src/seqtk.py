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

import pandas as pd
import glob
from kraken_out_parser import *
from seqtk_wrapper import *

#Do klasy_czas na ten skrypt

def taxonomic_level_to_dict(df):
    reads_dict = pd.Series(df.Taxon.values,index=df.Read_1).to_dict()
    return reads_dict

#Tu nazwy
def parse_unique(df):
    unique_names = set(df["Taxon"].tolist())
    return unique_names

#Ta trzeba teraz naprawic
def segregate(df):
    reads = taxonomic_level_to_dict(df)
    segregated = {}
    unique = parse_unique(df)
    for i in unique:
        segregated[i] = []
    for key in reads:
        segregated[reads[key]].append(key)
    return segregated

#This function is ok
def files_preparation(df):
    dictionary = segregate(df)
    to_save = dictionary.items()
    for i in to_save:
        with open(str("ids_"+i[0]+".1"),"w") as f:
            for id in i[1]:
                f.write(id+"\n")
        print str(i[0])+" has been written! (R)"

        with open(str("ids_"+i[0]+".2"),"w") as f:
            for id in i[1]:
                new_id = id[:-2]
                f.write(new_id+".2"+"\n")
        print str(i[0])+" has been written! (F)"

#Chcialbym napisac workflowy wtedy taka funcja bylaby nie potrzebna
def fromdf2reads_with_files(filename,tax_level,with_unclassif):
    df = taxtree_of_specific_level_to_dataframe('../../kraken_out_masked',4,False)
    files_preparation(df)


#To jest zupelnie nie zrobione
#Tu obsluga globa + klasa seqtk?
#This first!
#Fastq only! - it might be better for assembly
def seqtk_handler(reads_1,reads_2,workDir,settings):
    #workDir is a folder which contains, id_reads, kraken_out etc.
    os.chdir(workDir)
    list_of_ids = glob.glob["ids_*"]
    read_1 = glob.glob["classif_R1."]
    read_2 = glob.glont["classifl_R2."]

    if len(read_1) == 0 or len(read_2) == 0:
        print " [ERROR] There is no files with reads, change the work directory. Check it!"
        sys.exit()

    if len(list_of_ids) == 0:
        print " [ERROR] There is no files with reads id, change the work directory. Check it!"
        sys.exit()
    for reads_1, reads_2 in itertools.combinations(reads_1, reads_2):
        if reads_1[:-2] == reads_2[:-2]:
            print "Running seqtk subseq for %s" % reads_1
            SeqtkRunner(settings).run_seqtk(read_1,reads_1,workDir,".1")
            print "Running seqtk subseq for %s" % reads_2
            SeqtkRunner(settings).run_seqtk(read_2,reads_2,workDir,".2")



if __name__ == "__main__":

    import os
    from time import gmtime, strftime
    import sys
    import argparse

    description = """

Version 1.01

This script is designed for taxonomic-wise reads extraction.

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com
"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)




    parser.add_argument('read_1', metavar='read_1', type=str)
    parser.add_argument('read_2', metavar='read_2', type=str)
    parser.add_argument('labels', metavar='labels', type=str)
    parser.add_argument('ext', metavar='ext', type=str)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

    file_preparation(args.labels)
    seqtk_handler(args.read_1,args.read_2,args.labels,args.ext)
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    print "Starting time: "+start
    print "Ending time: "+end
