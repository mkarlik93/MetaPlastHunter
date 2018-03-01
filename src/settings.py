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
            print "  Loaded settings.txt"
            dict = {}
            for i in f:
                splited = i.split("=")
                if line == splited[0]:
                    dict[splited[0]] = splited[1].strip("\n")
                    print "  Checking for %s" % (splited[0])
                    try:
                        subprocess.call([splited[1].strip("\n")+splited[0], '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
                        print "  Status OK synek!"
                    except:
                        print "  [Error] Make sure %s path is in settings.txt or was set correctly, synek." % splited[0]
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
