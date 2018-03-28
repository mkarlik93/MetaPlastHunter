#!/usr/bin/env python2.6

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
import matplotlib
matplotlib.use("Agg")

import pysam
from ete3 import Tree
from ete3 import  NCBITaxa
import networkx as nx
from matplotlib import pyplot as plt
from math import log as ln
import pandas as pd
import logging
import os
import sys

from settings import  Settings_loader
from cov import Coverage


logger = logging.getLogger('src.sam_analyzer')
logging.basicConfig(level=logging.INFO)


# TODO:

#Sprawdzic liczby readow
# plus wykresy

"""Voting algorithm was described below
#Algorytm glosujacy
#Potrzeba jeszcze danych o pokryciu
#Jesli para_1 oraz para_2 to ten sam takson i nie ma secondary, zostaw pary


#Jesli para_1 oraz para_2 to primary i ten sam takson i sa secondary -> LCA
#Jesli para_1 oraz para_2 to inny takson -> LCA

"""

class Sam_analyzer:

    def __init__ (self,seqid_map):

        #Cov_list ma byc domyslnie w miejscu operacji
        #"seqid2taxid.map"
        """Params

        param _tree: NCBI tree instance

        param _seqid2taxid: dictonary instance

        param _cov_dict: dictonary instance, contains coverage info for every reference genome
        produced by bbmap during final alignment.

     """
        self._tree = self.get_NCBI_tree(seqid_map)
        self._seqid2taxid = self.seqid2taxid(seqid_map)
        self._cov_dict = self.cov_average("cov_list.txt")


    def get_NCBI_tree(self,seqidmap):

        ncbi = NCBITaxa()
        taxids = []
        with open(seqidmap,"r") as f:
            for line in f:
                splited = line.split("\t")
                taxids.append(splited[1].strip("\n"))
        tree = ncbi.get_topology(taxids,intermediate_nodes=True)
        return tree

    def get_lineage(self,taxid):

        ncbi = NCBITaxa()
        lineage =  ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        lin = [names[taxid] for taxid in lineage]
        return lin

    def extract_id(self,string):
        "In case of custom - chloroplast genome"

        if "chloroplast_seq" not in string:
            splited = string.split(".")
            split_1 = splited[1].split(" ")
            id = splited[0]+"."+split_1[0]

        else:

            id = string.split(" ")[0]

        return id

    def cov_average(self,cov_file):

        cov_dict = {}
        with open(cov_file,"r") as f:
            for line in f:
                splited = line.split(",")
                try:
                    cov_dict[self.extract_id(splited[0])] = splited[2].strip("\n")
                except IndexError:
                    cov_dict[self.extract_id(splited[0])] = splited[1].strip("\n")

        return cov_dict

    def get_NCBI_tree(self,seqidmap):

        ncbi = NCBITaxa()
        taxids = []
        with open(seqidmap,"r") as f:
            for line in f:
                splited = line.split("\t")
                taxids.append(splited[1].strip("\n"))
        tree = ncbi.get_topology(taxids,intermediate_nodes=True)
        return tree


    def seqid2taxid(self,seqidmap):

        seqid2taxid = {}
        with open(seqidmap,"r") as f:
            for line in f:
                splited = line.split("\t")
                seqid2taxid[splited[0]] = splited[1].strip("\n")

        return seqid2taxid

    def transform_read_name(self,readname):
        split_1 = readname.split(".")
        split_2 = split_1[2].split(" ")
        return  "%s.%s.%s" % (split_1[0],split_1[1],split_2[0])

    def sam_parse(self,samfile):

        samfile = pysam.AlignmentFile(samfile)
        dict_reads = {}

        for read in samfile.fetch():
            try:
                if read.is_unmapped:
                    pass
                else:


                    #if read.qname == "ERR538178.142953062.1 H3:C2FGHACXX:1:2310:4346:40753":
                    new_read_name =  self.transform_read_name(read.qname)
                    record = new_read_name[:-2]
                    if new_read_name[-2:] == ".2":
                        pass

                    else:
                        if record not in dict_reads:
                            dict_reads[record] = ([],[])
            #                dict_reads[record].append([read])
            #            else:
            #                dict_reads[record].append([read])

                    if read.is_read1:

                        if read.is_secondary:
                            if new_read_name[-2:] == ".2":
                                pass
                            else:
            #                    print new_read_name+".1",read.reference_name,"secondary"
                                dict_reads[record][0].append(read)
                        else:
            #                print new_read_name+".1",read.reference_name,"primary"
                            dict_reads[record][0].append(read)

                    if read.is_read2:

                        if read.is_secondary:

                            if new_read_name[-2:] == ".2":
                                pass
                            else:
            #                    print new_read_name+".2",read.reference_name,"secondary"
                                dict_reads[record][1].append(read)

                        else:
            #                print new_read_name+".2",read.reference_name,"primary"
                            dict_reads[record][1].append(read)

            #        print read.qname ,read.reference_name

            except ValueError:
                pass
        return dict_reads
    #    samfile.close()

    #    return dict_reads

    def reads_processing(self,samfile,out_file_name):


        logger.info("Processing of sam file has started")
        dict_of_reads = self.sam_parse(samfile)

        with open(out_file_name,"w") as f:

            for record in dict_of_reads:

                r1 = dict_of_reads[record][0]
                r2 = dict_of_reads[record][1]

                try:
                    if len(r1) == 1 and len(r2) == 1:

                        if r1[0].reference_name == r2[0].reference_name:

                            # Output ->

                            f.write("%s\t%s\n" % (record+".1", ";".join(self.get_lineage(self._seqid2taxid[self.extract_id(r1[0].reference_name)]))))

                            f.write("%s\t%s\n" % (record+".2", ";".join(self.get_lineage(self._seqid2taxid[self.extract_id(r1[0].reference_name)]))))

    #                        print record+".1", get_lineage(_seqid2taxid[extract_id(r1[0].reference_name)])
    #                        print record+".2", get_lineage(_seqid2taxid[extract_id(r1[0].reference_name)])



                        else:

                           lca = self._tree.get_common_ancestor(self._seqid2taxid[self.extract_id(r1[0].reference_name)], self._seqid2taxid[self.extract_id(r2[0].reference_name)])

                           f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(lca.taxid))))
                           f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(lca.taxid))))

    #                      print record+".1", get_lineage(lca.taxid)
    #                      print record+".2", get_lineage(lca.taxid)

                    elif len(r1) > 1 and len(r2) == 0:

                        lca = self._tree.get_common_ancestor([self._seqid2taxid[self.extract_id(read.reference_name)] for read in r1])

                        f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(lca.taxid))))
                        #print  record+".1", get_lineage(lca.taxid)

                    elif len(r1) == 0 and len(r2) == 1:

                        f.write("%s\t%s\n" % (record+".2", ";".join(self.get_lineage(self._seqid2taxid[self.extract_id(r2[0].reference_name)]))))
    #                    print record+".2", get_lineage(_seqid2taxid[extract_id(r2[0].reference_name)])



                    elif len(r1) == 0 and len(r2) > 1:
                        lca = self._tree.get_common_ancestor([self._seqid2taxid[self.extract_id(read.reference_name)] for read in r2])

    #                    print  record+".2", get_lineage(lca.taxid)
                        f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(lca.taxid))))

                    elif len(r1) > 1 and len(r2) >= 1:

                        #TU bedzie LCA z wagami -> suma pokrycia genomu pozostalych referencji musi byc wieksza niz najwiekszego z postostalych

                        #Linking with coverage
                        #Case when it's only one mate
            #            if len(r2) == 1:


                        _r1_taxid = [self._seqid2taxid[self.extract_id(read.reference_name)] for read in r1]
                        _r1_cov = [float(self._cov_dict[self.extract_id(read.reference_name)]) for read in r1]
                        _max_r1 = max(_r1_cov)
                        index_max_r1 = _r1_cov.index(_max_r1)

                        without_max = _r1_cov[:index_max_r1] + _r1_cov[index_max_r1+1 :]

                        if len(r2) == 1:

                            _r2 = self._seqid2taxid[self.extract_id(r2[0].reference_name)]
                            if (float(_max_r1) -  sum(cov for cov in without_max)) < 0:

                                lca = self._tree.get_common_ancestor(_r1_taxid + [_r2])

                                f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(lca.taxid))))
                                f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(lca.taxid))))

                                #print record+".1", get_lineage(lca.taxid)
                                #print record+".2", get_lineage(lca.taxid)


                            else:
                                if _r1_taxid[index_max_r1] == _r2:

                                    f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(_r1_taxid[index_max_r1]))))
                                    f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(_r1_taxid[index_max_r1]))))

    #                                print record+".1", get_lineage(_r1_taxid[index_max_r1])
    #                                print record+".2", get_lineage(_r1_taxid[index_max_r1])

                                else:

                                    lca = self._tree.get_common_ancestor(_r1_taxid[index_max_r1],_r2)


                                    f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(lca.taxid))))
                                    f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(lca.taxid))))


    #                                print  record+".1", get_lineage(lca.taxid)
    #                                print  record+".2", get_lineage(lca.taxid)

                        if len(r2) > 1:

                            _r2_taxid = [self._seqid2taxid[self.extract_id(read.reference_name)] for read in r1]
                            _r2_cov = [float(self._cov_dict[self.extract_id(read.reference_name)]) for read in r1]
                            _max_r2 = max(_r1_cov)
                            index_max_r2 = _r1_cov.index(_max_r1)

                            without_max_r2 = _r2_cov[:index_max_r1] + _r1_cov[index_max_r2+1 :]

                            if (float(_max_r1) -  sum(cov for cov in without_max)) > 0 and  (float(_max_r2) -  sum(cov for cov in without_max_r2)) > 0:

                                if _r1_taxid[index_max_r1] == _r2_taxid[index_max_r2]:

                                    f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(_r1_taxid[index_max_r1]))))
                                    f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(_r1_taxid[index_max_r1]))))

                                    #print record+".1", get_lineage(_r1_taxid[index_max_r1])
                                    #print record+".2", get_lineage(_r1_taxid[index_max_r1])

                                else:

                                    lca = self._tree.get_common_ancestor(_r1_taxid[index_max_r1],_r2_taxid[index_max_r2])

                            else:

                                 lca = self._tree.get_common_ancestor(_r1_taxid + _r2_taxid)
                                 f.write("%s\t%s\n" % (record+".1",";".join(self.get_lineage(lca.taxid))))
                                 f.write("%s\t%s\n" % (record+".2",";".join(self.get_lineage(lca.taxid))))

    #                             print record+".1", get_lineage(lca.taxid)
    #                             print record+".2", get_lineage(lca.taxid)


                except ValueError:
                    pass

" LCA POSTPROCESS "


class LCA_postprocess:

    def __init__ (self,i,treshold):

        """Params

        param treshold: To set min number of readsa are needed for support given taxa
        during graph prunning.

        param _lca_graph:  Lca graph instance, takes file produced by Sam_analyzer class and run lca_graph function

        param _graph: Transformed graph instance

        """

        self.treshold = treshold
        self._lca_graph = self.lca_graph(i+"_taxonomic_assignment.txt")
        self._graph = self.lca_graph_analysis(treshold)
        self._sample_name = i


    def get_count(self,filename):
        num_lines = sum(1 for line in open(filename))
        return num_lines

    def open_lca_file(self,filename):

        try:
            count = self.get_count(filename)
        except IOError:
            logger.error("There is no taxonomic_assigment in working directory")
            sys.exit()


        taxa_count_dict = {}


        with open(filename,"r") as f:
            longest = 0
            for line in f:
                count = count + 1
                # first Read name, second lineage

                splited = line.split("\t")

                if len(splited[1].split("/")) > longest:
                    longest = len(splited[1].split(";"))

                if splited[1].strip("\n") not in taxa_count_dict:
                    taxa_count_dict[splited[1].strip("\n")] = 1

                else:
                    taxa_count_dict[splited[1].strip("\n")] = taxa_count_dict[splited[1].strip("\n")] + 1

            return taxa_count_dict,longest

    def lca_graph(self,filename):


        G = nx.Graph()

        taxa_count_dict, longest =  self.open_lca_file(filename)

        good_taxa = {}
        bad_taxa = {}

        G.add_node('root',count = self.get_count(filename))

        for key in taxa_count_dict:


            lineage = key.split(";")
            #edges and nodes
            for i in xrange(len(lineage)-1):

                if lineage[i+1] not in G.nodes():
                    G.add_node(lineage[i+1],count=int(taxa_count_dict[key]))
    #                G[lineage[i+1]]['count'] = taxa_count_dict[key]

                else:

    #                print G[lineage[i+1]]
                    old = G.nodes[lineage[i+1]]['count']
                    G.nodes[lineage[i+1]]['count'] = old + taxa_count_dict[key]
                G.add_edge(lineage[i],lineage[i+1])

        return G

    def lca_graph_analysis(self,treshold):

#        treshold = self.treshold
        count = self._lca_graph.nodes['root']['count']

        _lca_graph_for_iteration = self._lca_graph
        species_treshold = int(count*treshold)
        #node pruning
        to_del = []
        for node in self._lca_graph.nodes():
            if self._lca_graph.nodes[str(node)]['count'] < species_treshold:
                to_del.append(node)
        self._lca_graph.remove_nodes_from(to_del)

        return self._lca_graph

    #Trzeba zapamietywac kazdy wezel w grafie
    #Weight -> moze to byc np. liczba nodow
    #Node tez moze miec atribute

    #Takes graph and return counts of taxa on species level

    def species_level(self):

        " Returns dictonary of counted reads on species level "

        species = {}
        for node in self._graph.nodes():
            if self._graph.degree[node] == 1 and len(node.split(" ")) == 2:
                species[node] = self._graph.nodes[node]['count']
        return species


    def count_plot_species_level(self):

        " Plots barplot of counted reads on species level "
        #Plot_counts

        _species_level = self.species_level()
        plt.bar(range(len(_species_level)), list(_species_level.values()), align='center')
        plt.xticks(range(len(_species_level)), list(_species_level.keys()),rotation=20)
        plt.ylabel("counts")
        plt.xlabel("Species names")
        plt.title("Read counts")
        plt.savefig(self._sample_name+"_species_level.png")

    def pandas_data_frame_species_level(self):

        _species_level = self.species_level()
        df = pd.DataFrame(_species_level.items(),columns=['Taxon', 'Counts'])
        df.to_csv(self._sample_name+"_species_level_table.csv")
#        return df

    #Trzeba dorobic!
    def pandas_data_frame_fourth_level(self):

        #TODO - dodac tu wszystkie z czwartego poziomou
        taxa_to_catch = ["Stramenopiles","Haptophycae","Viridiplantae","Chromista"]
        _species_level = self.species_level()
        df = pd.DataFrame(_species_level.items(),columns=['Taxon', 'Counts'])
        return df

    def tree_reconstruction_of_nodes(self):

        " Returns lineage of node in graph of degree one, except the root"

        avaible_paths = []
        for node in self._graph.nodes():
            if self._graph.degree[node] == 1 and node != 'root':
                avaible_paths.append(nx.shortest_path(self._graph,source="root",target=node))
        return avaible_paths



    def species_level_shannon_index(self):

        " Writes shannon index [shannon.txt] for given sample"

        df_dict = self.species_level()

        def p(n, N):
            """ Relative abundance """
            if n is  0:
                return 0
            else:
                return (float(n)/N) * ln(float(n)/N)

        N = sum(df_dict.values())

        shannon = -sum(p(n, N) for n in df_dict.values() if n is not 0)
        #print "The species level shanon index is %s" % (shannon)
        with open("shannon.txt","w") as f:
            f.write("Shannon index\t%s" % (shannon))
#        return shannon

class Run_analysis_sam_lca:

    def __init__ (self,  list_sra, station_name,settings):

        self.list_sra = list_sra
        self.station_name = station_name
        self.settings = settings

        try:

            self.seqidmap = Settings_loader(mode="seqid2taxid.map",path=self.settings).read_database()["seqid2taxid.map"]
            self.lca_treshold = float(Settings_loader(mode="lca_treshold",path=self.settings).read_parameters()["lca_treshold"])

        except KeyError:

            logger.error("Something went wrong during loading files, check paths")
            sys.exit()

    def process(self):

        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        starting_dir = os.getcwd()
        for i in list_sra_ids:
            logger.info("Procesing "+i)
            dir = "%s/%s/" % (self.station_name,i)
            os.chdir(dir)

            Coverage('bincov.txt',i+"_chloroplasts.hitstats",self.settings).report_cov()
            Sam_analyzer(self.seqidmap).reads_processing(i+"_final_mapped.sam",i+"_taxonomic_assignment.txt")

            lca_postprocess =  LCA_postprocess(i,self.lca_treshold)
            lca_postprocess.pandas_data_frame_species_level()
            lca_postprocess.species_level_shannon_index()
            lca_postprocess.count_plot_species_level()

            os.chdir(starting_dir)


#ANALIZA GRAFU
#print(nx.node_connectivity(_graph,"root","Emiliania huxleyi"))
#betweeness = nx.betweenness_centrality(_graph)

#print sorted(betweeness.values())

# miara dywersyfikacji na podstawie stopnia wezla w grafie
#print species_level()
#for node in graph.nodes():
#    if graph.degree[node] >= 3:
#        print node, graph.nodes[node]['count'],graph.degree[node]

#degree_dict = dict(graph.degree(graph.nodes()))
#print nx.set_node_attributes(graph, degree_dict, 'degree')

# Wychodzi , ze posiadajac stosunkowo malo odczytow, ale rownomiernie rozlozonych mozna latwo przejsc "bramke" pokryciowa
