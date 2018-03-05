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

from Bio import Entrez
import glob
import os

""" script for chloroplast genomes database  """

def accession_from_file(file):
    file = open(file,"r")
    accession_numbers = []
    for i in file:
        splited = i.split("\t")
        accession_numbers.append(splited[1])


class NCBI_fetch:

    def __init__(self, list_of_ids, path, format, number=None):
        self.list_of_ids = list_of_ids
        self.path = path
        self.format = format
        self.number = number

    def make_db_direcory():
        command = "mkdir "



    def fetcher(self):
        Entrez.email = "mich1@wp.pl"
        number = self.number
        path = self.path
        format = self.format
        list1 = accession_from_file(self.list_of_ids)
        count = 0
        if number == None:
            final_list = list1
        else:
            final_list = random.sample(list1, number)
        for genome_id in final_list:
            count = count + 1
            handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
            record = Entrez.efetch(db="nucleotide", id=genome_id, rettype=str(format), retmode="text")
            x = SeqIO.read(handle, 'genbank')
            name = x.annotations['organism']
            name = name.replace(" ","_")
            name = name.replace("/","")
            filename = name+"_"+str(count)+"_.gb"
            with open(path+filename, 'w') as f:
                f.write(record.read())


"""???"""
    def run_repeat_masker():
        pass



    def cat_files_into_one():
        pass


if __name__ == "__main__":

    from time import gmtime, strftime
    import sys
    import os
    import argparse
    from Bio import Entrez
    from Bio import SeqIO
    import random

    description = """

Version 1.00

Script designed for downloading sequences using NC ids

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com
"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)




    parser.add_argument('file_with_ids', metavar='-file_with_ids', type=str)
    parser.add_argument('location', metavar='-location', type=str)
    parser.add_argument('number', nargs='?', type=int,default=None)
    parser.add_argument('format', nargs='?', type=str,default="fasta")


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    run = NCBI_fetch(args.file_with_ids,args.location,args.format,args.number)
    run.fetcher()

    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    print "Starting time: "+start
    print "Ending time: "+end
