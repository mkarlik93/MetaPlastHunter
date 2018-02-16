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
import pandas as pd
import numpy as np
import seaborn
import matplotlib.pyplot as plt

#definitition of unclassified chloroplast
#let's try with these settings
#threshold 0.25
# < unclassified chloroplast

#taxon dictionary
#taxon_levels = {
#"root":1,
#"cellular organisms":2,
#"Eukaryota":3,
# : ,
#4 : ,
#5 : ,
#}


def filtering(kraken_out):
    with open(kraken_out, "r") as f:
        for line in f:
            splited = line.split("\t")
            if splited[0] == "C":
                yield splited

def score(line):
    score_part = line[4].strip("\n")
    score_part = score_part.split(" ")
    assigmnets = {}
    for i in score_part:
        kmers_assignment = i.split(":")
        if kmers_assignment[0] not in assigmnets:
            assigmnets[kmers_assignment[0]] = int(kmers_assignment[1])
        else:
            assigmnets[kmers_assignment[0]] = int(assigmnets[kmers_assignment[0]]) + int(kmers_assignment[1])
    classified = 0
    max_counts = (0,"")
    for key in assigmnets:
        if key != '0' and key != "A":
            if assigmnets[key] > max_counts[0]:
                max_counts = (assigmnets[key],key)
            classified = classified + int(assigmnets[key])
    #this line can be improved
    try:
        score = int(classified)/float(assigmnets["0"]+float(classified))
    except KeyError:
        score = 1

    return [line[1],line[2],str(score),max_counts[1]]

def scoring(list):
    for i in list:
        yield score(i)

def save_file(generator):
    with open("filered_records.txt","w") as f:
        for i in generator:
            f.write(",".join[i]+"\n")

def parse_tax_names_2_dict(namesdmp):
    names = {}
    with open(namesdmp,"r") as f:
        for line in f:
            splited = line.split("|")
            if splited[3] == '\tscientific name\t':
                names[splited[0].strip("\t")] = splited[1].strip("\t")
    return names

def parse_lines_to_names(names_dict,list,threshold):
    dict_init = names_dict
    for record in list:
        adjusted = record[1]
        maximal = record[3]
        if float(record[2]) >= threshold:
            try:
                yield [record[0],dict_init[adjusted],record[2],dict_init[adjusted]]
            except KeyError:
                yield [record[0],str(adjusted),record[2],str(adjusted)]
        else:
            try:
                yield [record[0],"Unclassified chloroplast sequence",record[2],dict_init[adjusted]]
            except KeyError:
                yield [record[0],"Unclassified chloroplast sequence",record[2],str(adjusted)]

def parsed_names_to_csv(filename,generator):
    cwd = os.getcwd()
    with open(str(filename),"w") as f:
        header = "Name of pair,assigned_name,score,kmer_max_count\n"
        f.write(header)
        for line in generator:
            f.write(",".join(line)+"\n")
        print "The file has been written in "+str(cwd)+" directory"

#node.dmp parsing part. From output to dataframe and plots
def parse_tax_tree__2_dict(taxdmp):
    taxtree = {}
    with open(taxdmp,"r") as f:
        for line in f:
            splited = line.split("|")
            taxtree[splited[0].strip("\t")] = splited[1].strip("\t")
    return taxtree

def tree_walk_recurse(taxdict,taxid,list):
    if taxid == '1':
        list.append("1")
        return list
    else:
        list.append(taxid)
        try:
            return tree_walk_recurse(taxdict,taxdict[taxid],list)
        except KeyError:
            if len(list) == 0:
                return "No tax walk for this taxid"

def walk2names(walk,names_dict):
        try:
            walk = walk[::-1]
            for taxid in walk:
                try:
                    yield names_dict[taxid]
                except KeyError:
                    yield "No tax walk for this taxid"
        except TypeError:
            yield "No tax walk for this taxid"

def multiple_walks(walks,names_dict):
    for walk in walks:
        yield walk2names(walk,names_dict)

def line_with_tree(parse_lines_to_names,taxdict,names_dict,threshold):
    for i in parse_lines_to_names:
        if float(i[2]) >= float(threshold):
            tax = [list(walk2names(tree_walk_recurse(taxdict,i[1],[]),names_dict))]
            yield [tax,i[0]]
        else:
            yield [i[0],"Unclassified chloroplast sequence"]

def parse_2_names(filename,out):
    taxdict = parse_tax_tree__2_dict("nodes.dmp")
    names_dict = parse_tax_names_2_dict("names.dmp")
    parsed_names_to_csv(out,parse_lines_to_names(names_dict,scoring(filtering(filename)),0.2))


def parsed_tree_to_csv(filename,generator):
    cwd = os.getcwd()
    with open(str(filename),"w") as f:
        header = "tree,readname\n"
        f.write(header)
        for line in generator:
            f.write(",".join(line)+"\n")
        print "The file has been written in "+str(cwd)+" directory"


def parse_2_taxtree(filename):
    taxdict = parse_tax_tree__2_dict("../../nodes.dmp")
    names_dict = parse_tax_names_2_dict("../../names.dmp")
    line_with_3 = line_with_tree(scoring(filtering(filename)),taxdict,names_dict,1)
    return list(line_with_3)

#Data visualization + counting

def extract_level_from_root_all(taxtree_list, level_from_root, with_unclassif):
    if with_unclassif == True:
        level_from_root = level_from_root -1
        for line in taxtree_list:
            try:
                if line[1] == "Unclassified chloroplast sequence":
                    yield line
                else:
                    if line[0][0][level_from_root] == "Viridiplantae":
                        yield [line[0][0][level_from_root+1],line[1]]
                    else:
                        yield [line[0][0][level_from_root],line[1]]
            except IndexError:
                yield ["Not all of your reads reach this taxonomic level, take another one",line[1]]
    else:
        level_from_root = level_from_root -1
        for line in taxtree_list:
            try:
                if line[1] == "Unclassified chloroplast sequence":
                    pass
                else:
#                    print line
                    if line[0][0][level_from_root] == "Viridiplantae":
                        yield [line[0][0][level_from_root+1],line[1]]
                    else:
                        yield [line[0][0][level_from_root],line[1]]
            except IndexError:
                yield ["Not all of your reads reach this taxonomic level, take another one",line[1]]

#because taxonomic tree of eukaryota is complicated

def extract_species_level(taxtree_list):
    for line in taxtree_list:
        if line[1] == "Unclassified chloroplast sequence":
            pass
        else:
            if len(line[0][0][-1].split(" ")) == 2:
                yield [line[0][0][-1],line[1]]
            else:
                pass

def taxtree_of_specific_level_to_dataframe(filename,distance_from_root,with_unclassif):
    a = extract_level_from_root_all(parse_2_taxtree(filename),distance_from_root,with_unclassif)
    df = pd.DataFrame(a,columns=["Taxon","Read_1"])
    return df

#tables
def count_plot(df):
    seaborn.countplot(y="Taxon",data=df)
    plt.show()

def counts2table(df):
    counts = df["Taxon"].value_counts()
    return counts
    #Prettytable

def percentage_df(df):
    counts = counts2table(df)
    counts = dict(counts)
    df = pd.DataFrame({'Taxon' : counts.keys() , 'counts' : counts.values() })
    return df

#Plots
def percentage_plot_y_oriented(df):
    all = sum(df['counts'])
    ax = seaborn.barplot(y="Taxon",x="counts",data=df, estimator= lambda y: y / float(all) * 100)
    ax.set(xlabel="Percent",ylabel="Taxon")
    plt.show()

def percentage_plot_x_oriented(df):
    all = sum(df['counts'])
    ax = seaborn.barplot(x="Taxon",y="counts",data=df, estimator= lambda x: x / float(all) * 100)
    ax.set(xlabel="Percent",ylabel="Taxon")
    plt.show()

#prepares dict of counts of species
def species4shannon_index(filename,with_unclassif):
    a = extract_species_level(parse_2_taxtree(filename),with_unclassif)
    df = pd.DataFrame(a,columns=["Taxon","Read_1"])
    counts = counts2table(df)
    counts = dict(counts)
    return counts


#a = taxtree_to_dataframe('../../kraken_out_masked',4,False)
#percentage_plot(a)
#new_data_frame =  pd.DataFrame(series2table(a),columns=["Taxon","Readcounts"])
#percent = percentage_df(a)
#percentage_plot(percent)


#print count_plot(a)

#General function
def analyze():
    pass




#parse_2_names("kraken_out_masked","test_filtered_names.csv")



#a =  filtering("Sample_kraken_out.txt")
#scored = list(scoring(a))
#record = list(parse_lines_to_names('names.dmp',scored,0.2))

#print  list(line_with_tree(record,taxdict,names_dict))

#TODO
#Sprawdzic co sie dzieje z odczytami ktore nie maja taksonomii z jakiegos powodu
#Slownik z poziomem taksonomicznym
#Zrobic by nie bylo level, tylko nazwa poziomu taksonomiczego
#Pomyslec nad stopiem integracji z pythonem, biopythonem
