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


#import skbio.diversity
from math import log as ln
from kraken_out_parser import *
# na poziomie gatunku!!!!




def species_level_shannon_index(df_dict):
    """ Given a hash { 'species': count } , returns the SDI

    >>> sdi({'a': 10, 'b': 20, 'c': 30,})
    1.0114042647073518"""

    def p(n, N):
        """ Relative abundance """
        if n is  0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)

    N = sum(df_dict.values())

    return -sum(p(n, N) for n in df_dict.values() if n is not 0)

species = species4shannon_index('../../kraken_out_masked')
print species
print species_level_shannon_index(species)
