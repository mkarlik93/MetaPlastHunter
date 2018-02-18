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

#Basemap powered part
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd

#for_example_tara_oceans
#We need data table like this:
#        sample_1 sample_2 sample_3
# species
# species
# species

# samples column have to be connected

#TODO
#What do we want?

#take_data
def fromcsv2df(filename):
    df = pd.DataFrame.from_csv(filename, header=0,sep=",", index_col=0)
    return df

def dataframe_deduplication(df):
    deduplicated = df.drop_duplicates(subset='Station identifier [TARA_station#]', keep='first')
    return deduplicated

def load_tara_data():
    my_dataframe = fromcsv2df("../data/tara.csv")
    return my_dataframe
#By jedna nazwe reprezentowaly jeden wspolrzedne (male roznice)

def tara_singletons():
    deduplicated = dataframe_deduplication(my_dataframe)
    return deduplicated

def show_df_headers(df):
    headers = list(df.columns.values)
    return headers

def tara_oceans_lon_lat(df):
    tara_lat = df['Latitude [degrees North]'].tolist()
    tara_lat = map(lambda x: x.replace(",","."), tara_lat)
    tara_lon = df['Longitude [degrees East]']
    tara_lon = map(lambda x: x.replace(",","."), tara_lon)
    return tara_lon, tara_lat

tara_lon, tara_lat = tara_oceans_lon_lat(deduplicated)
#settings = {"world":}


def draw_map(tara_lon, tara_lat):
    #background map
    fig = plt.figure(figsize=(8, 8))
    m = Basemap(projection='cyl', resolution='i',
             llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180,)

    m.shadedrelief()
    m.drawcoastlines(color='gray')
    #m.drawcountries(color='gray')
    #m.drawstates(color='gray')
    m.scatter(tara_lon, tara_lat)
    plt.show()
