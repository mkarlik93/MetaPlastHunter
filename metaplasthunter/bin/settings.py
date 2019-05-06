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



    """
    Parameters
    ----------
    name : path

        keeps path to the settings file in YAML format

    """


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

    def check_reads(self,reads):
        for key in reads:

            if os.path.isfile(reads[key]) or os.path.isdir(reads[key]):
                logger.info("%s read file has been added correctly" % key)
            else:
                logger.error("Check path of %s" % key )
                sys.exit()

    def check_contig(self,contig):
        for key in contig:

            if os.path.isfile(contig[key]) or os.path.isdir(contig[key]):
                logger.info("%s contig file has been added correctly" % key)

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

                print exc
                sys.exit()
        return settings


    def yaml_check_settings_file(self):

        settings_to_check = self.yaml_handler()
        self.check_software(settings_to_check["Software dependencies"])
        self.check_db(settings_to_check["Databases and mapping files"])
        self.check_params(settings_to_check["Params"])
#        self.check_reads(settings_to_check["Input reads"])
#        self.check_contig(settings_to_check["Input contigs"])

    #Tu bym zrobil jedna wielka funkcje
    def check_software_classification(self,software):

        """Check to see if software is on the system before we try to run it."""

        classification_list = ["bbduk.sh","bbmap.sh"]

        "for reconstruct pipeline"
        reconstruction_list = []

        for key in software:
            if key in classification_list:
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
            else:
                pass

    def check_params_classification(self,params):

        classification_list = ["min_identity","static_coverage_threshold[%]","lca_threshold","min_bin_coverage","bincov4_report"]
        for key in params:
            if key in classification_list:
                if params[key] == None:
                    logger.error("Param cannot be empty! Please fill those fields")
                    sys.exit()
                else:
                    logger.info("%s was set with value %s   " % (key,params[key]))
            else:
                pass

    def check_db_classification(self,db):

        classification_list = ["bbmap_base","seqid2taxid.map","silva"]
        for key in db:

            if key in classification_list:

                if os.path.isfile(db[key]) or os.path.isdir(db[key]):
                    logger.info("%s database has been added correctly" % key)

                else:

                    logger.error("Check path of %s" % key )
                    sys.exit()
            else:
                pass

    def yaml_check_settings_file_classification(self):
        settings_to_check = self.yaml_handler()
        self.check_software_classification(settings_to_check["Software dependencies"])
        self.check_params_classification(settings_to_check["Params"])
        self.check_db_classification(settings_to_check["Databases and mapping files"])
