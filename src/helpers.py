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



#collects table from file and merge it into one large table
def collector_single_specie_oriented(folders, specie_to_find):
    pass



class SettingsError(BaseException):
    pass

class Settings_loader:
    def __init__(self,mode="kraken"):
        if mode == "kraken":
            self.mode = 'kraken'
        elif mode == "seqtk":
            self.mode = 'seqtk'
#in case of more software
#        elif mode == 'align':
#            self.mode = 'align'
#        elif mode == 'fetch':
#            self.mode = 'fetch'
        else:
            raise SettingsError("Mode %s not understood" % mode)

        self.check_settings()


    def check_settings(self):
        if glob.glob("../settings.txt") == 0:
            print "[ERROR] There is no settings file"
            sys.exit()


    def read_path(self):
        if self.mode == 'kraken':
            line = 'kraken'
        elif self.mode == 'seqtk':
            line = 'seqtk'

        with open("../settings.txt") as f:
            dict = {}
            database = "kraken_db"
            for i in f:
                splited = i.split("=")
                if line == splited[0] or database == splited[0] :
                    dict[splited[0]] = splited[1].strip("\n")
            return dict
