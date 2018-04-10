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

import os
import sys
import subprocess
import logging
import yaml

logger = logging.getLogger('src.settings')
logging.basicConfig(level=logging.INFO)


#class SettingsError(BaseException):
#    pass


#    raise SettingsError("Mode %s not understood" % mode)


#YAML


class Settings_loader_yaml:
    def __init__(self,path):
        self.path = path

    def check_software(self,software):
        """Check to see if software is on the system before we try to run it."""

        for key in software:

            if software[key] == None:

                try:
                    subprocess.call([key, '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
                    logger.info("%s has been installed correctly!" % key)
                except:

                    logger.error("Make sure %s is on your system path or set propere to path in settings.txt" % key)
                    sys.exit()
            else:

                try:
                    subprocess.call([software[key]+key, '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
                    logger.info("%s has been installed correctly!" % key)
                except:

                    logger.error("Make sure %s is on your system path or set propere to path in settings.txt" % key)
                    sys.exit()

    def check_db(self,db):

        for key in db:

            if os.path.isfile(db[key]) or os.path.isdir(db[key]):
                logger.info("%s database has been added correctly" % key)

            else:

                logger.error("Check path of %s" % key )
                sys.exit()

    def check_params(self,params):

        for key in params:
            if params[key] == None:
                logger.error("Param cannot be empty! Please fill those fields")
                sys.exit()
            else:
                logger.info("Param %s is fine" % key)

    def yaml_handler(self):
        with open(self.path, 'r') as stream:
            try:
                settings = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return settings


    def yaml_check_settings_file(self):

        settings_to_check = self.yaml_handler()
        self.check_software(settings_to_check["Software dependencies"])
        self.check_db(settings_to_check["Databases and mapping files"])
        self.check_params(settings_to_check["Params"])
