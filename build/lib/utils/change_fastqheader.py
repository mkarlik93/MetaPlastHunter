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


def change_header(read,number_of_pair):
    header_changed_read = []
    for seq_record in SeqIO.parse(read,"fastq"):
        old_header = seq_record.id
        new_header = old_header[:-2]+"/"+str(number_of_pair)
        seq_record.id = new_header
        seq_record.description = ""
        header_changed_read.append(seq_record)
    SeqIO.write(header_changed_read, str(read).rstrip(".fastq")+"_changed.fastq","fastq")
    print "Headers have been changed"

def paired_header(read_1, read_2):
    change_header(read_1,1)
    change_header(read_2,2)

if __name__ == "__main__":


    from Bio import SeqIO
    import os
    from time import gmtime, strftime
    import sys
    import argparse

    description = """

Version 1.0

Script changes fastq headers.

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com



"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)




    parser.add_argument('read_1', metavar='-read1', type=str)
    parser.add_argument('read_2', metavar='-read2', type=str)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    paired_header(args.read_1,args.read_2)
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    print "Starting time: "+start
    print "Ending time: "+end
