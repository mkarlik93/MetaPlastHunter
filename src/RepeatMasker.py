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


import subprocess

""" This is wrapper for RepeatMasker """


class RMasker(BaseException):
    pass


def __init__(self,settings):
    self.settings = settings

    if settings == None:
        self.checkForRMasker()
        self.path = ""
    else:

        self.path = Settings_loader(mode="RepeatMasker",path=self.settings).read_path()["RepeatMasker"]

def run_masking(self,fastafile):
    path = self.path
    command = "%sRepeatMasker %s -small" % (path, fastafile)
    os.system(command)


def checkForRMasker(self):
    """Check to see if RepeatMasker is on the system before we try to run it."""

    try:
        subprocess.call(['RepeatMasker', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    except:
        print "  [Error] Make sure RepeatMasker is on your system path or set usage to path in settings.txt"
        sys.exit()
