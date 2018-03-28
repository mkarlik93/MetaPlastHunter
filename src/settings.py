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
import os
import sys
import subprocess
import logging




class SettingsError(BaseException):
    pass

class Settings_loader:

    def __init__(self,mode="kraken",path='path'):

        self.path = path

        if mode == "kraken":
            self.mode = 'kraken'
        elif mode == "seqtk":
            self.mode = 'seqtk'
#in case of more software
        elif mode == 'fastq-dump':
            self.mode = 'fastq-dump'
        elif mode == 'bbmap.sh':
            self.mode = 'bbmap.sh'

        elif mode == 'bbduk.sh':
            self.mode = 'bbduk.sh'

        elif mode == 'names.dmp':
            self.mode = 'names.dmp'

        elif mode == 'nodes.dmp':
            self.mode =  'nodes.dmp'

        elif mode == 'seqid2taxid.map':
            self.mode = 'seqid2taxid.map'

        elif mode == 'silva':
            self.mode = 'silva'

        elif mode == 'RepeatMasker':
            self.mode = 'RepeatMasker'

        elif mode == 'percentile_treshold':
            self.mode = 'percentile_treshold'

        elif mode == 'min_bin_coverage':
            self.mode = 'min_bin_coverage'

        elif mode == 'bincov4_report':
            self.mode = 'bincov4_report'

        elif mode == 'lca_treshold':
            self.mode = 'lca_treshold'

        else:
            raise SettingsError("Mode %s not understood" % mode)

        self.check_settings()


    def check_settings(self):
        if glob.glob("../settings.txt") == 0:
            print "  [ERROR] There is no settings file"
            sys.exit()


    def read_path(self):
        if self.mode == 'kraken':
            line = 'kraken'
        elif self.mode == 'seqtk':
            line = 'seqtk'
        elif self.mode == 'fastq-dump':
            line = 'fastq-dump'

        elif self.mode == "bbduk.sh":
            line = 'bbduk.sh'

        elif self.mode == 'bbmap.sh':
            line = 'bbmap.sh'


        with open(self.path) as f:
            logger.info("  Loaded settings.txt")
            dict = {}
            for i in f:
                splited = i.split("=")
                if line == splited[0]:
                    dict[splited[0]] = splited[1].strip("\n")
                    print "  Checking for %s" % (splited[0])
                    try:
                        subprocess.call([splited[1].strip("\n")+splited[0], '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
                        logger.info("  Status OK synek!")
                    except:
                        logger.info("  [Error] Make sure %s path is in settings.txt or was set correctly, synek." % splited[0])
                        sys.exit()
            return dict

    def read_database(self):

        if self.mode == 'kraken':
            line = 'kraken_db'
        elif self.mode == 'bbmap.sh':
            line = 'bbmap_base'
        elif self.mode == 'names.dmp':
            line = 'names.dmp'
        elif self.mode == 'nodes.dmp':
            line = 'nodes.dmp'
        elif self.mode == 'seqid2taxid.map':
            line = 'seqid2taxid.map'
        elif self.mode == 'silva':
            line = 'silva'

        with open(self.path) as f:
            dict = {}
            for i in f:
                splited = i.split("=")
                if line == splited[0]:
                    dict[splited[0]] = splited[1].strip("\n")
                    print "  Checking for %s" % (splited[0])
            return dict

    def read_parameters(self):

        if self.mode == 'percentile_treshold':
            line = 'percentile_treshold'

        elif self.mode == 'min_bin_coverage':
            line = 'min_bin_coverage'

        elif self.mode == 'bincov4_report':
            line = 'bincov4_report'

        elif self.mode == 'lca_treshold':
            line = 'lca_treshold'

        with open(self.path) as f:
            dict = {}
            for i in f:
                splited = i.split("=")
                if line == splited[0]:
                    dict[splited[0]] = splited[1].strip("\n")
                    logger.info("  Checking for %s" % (splited[0]))
            return dict
