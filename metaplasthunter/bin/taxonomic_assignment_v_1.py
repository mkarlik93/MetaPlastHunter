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



from ete3 import Tree
from ete3 import  NCBITaxa
import networkx as nx
import pysam
import numpy
import glob
import itertools
import sys
import logging
import pandas as pd
from settings import *
from cov import Coverage
from subprocess import Popen, PIPE


logger = logging.getLogger("src.taxonomic_assignment")
logging.basicConfig(level=logging.INFO)

def sam_get_min_identity_of_covered_genome(target_genome):

    """ Returns identity level of reads of covered genomes at species level min (95% ?????) """

    samfile = glob.glob("*_final_mapped.sam")[0]
    samfile = pysam.AlignmentFile(samfile)
    choosed_reads = []

    for read in samfile.fetch():

        try:
            if read.is_unmapped:
                pass

            else:

                if read.is_secondary:
                    pass

                else:
                    if read.reference_name == target_genome:
                        read_test = read
                        choosed_reads.append(float(read.tags[-1][1]))
        except:
            pass

    return numpy.mean(choosed_reads)

def calculate_average_nucleotide_identity(list_of_genomes):
    pass

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def sam_splitter(input,is_sam,list_of_genomes):

    if is_sam  == True:

        samfile =  pysam.AlignmentFile(input)
    else:
        _samfile = glob.glob("*_final_mapped.sam")[0]
        samfile = pysam.AlignmentFile(_samfile)

    sam_header = samfile.header

    for target_genome in list_of_genomes:

        print target_genome
        good_reads = []
        for read in samfile.fetch():
            try:

                if read.is_unmapped:
                    pass
                else:

                    if read.is_secondary:
                        pass

                    else:
                        if read.reference_name == target_genome:

                            choosed_reads.append(read)
            except:

                pass

        with pysam.AlignmentFile(target_genome+"_single.bam", "wb", header=sam_header) as outf:
            for read in good_reads:
                outf.write(read)

class Taxonomic_assignment(object):


    "Tu musze zmienic "


    ncbi = NCBITaxa()

    def __init__ (self,sample,treshold,seqidmap,sam_type,input):

        """
        Input Parameters
        ----------
        sam_type : boolean
            bincov.txt, draft coverage file, produced during mapping

        input: str
            absolute path of input file

        treshold : int
            XXX

        seqidmap : file
            General settings file which keeps hyperparameters

        Calculated gobal variables
        ----------
        _seqid: dictionary

        _sample_name: str
            project folder name (output)

        _seqid_inverted: dictionary

        _dict_of_reads: dictionary

        lca_assign: dictionary

        _lca_graph: graph data structure
            Orginal data structure (weighted digraph) which keeps NCBI tree with proposed taxonomic positions (LTU and TTU)

        analyzed_graf: graph data stucture
            Prunned _lca_graph that keeps only taxonomic positions that have been above the treshold (min support value)


        """
        self.sam_type = sam_type
        self.input = input
        self._seqid = self.seqid2taxid(seqidmap)
        self._sample_name = sample
        self._seqid_inverted = self._seqid_inverted(self._seqid)
        self._dict_of_reads = self.sam_parse()
        self.treshold = treshold
        self.lca_assign,self.almost_full = self.lca_assignment_just_taxids()
        self.run = self.process_taxonomic_assignment_to_file_easy()
        self._lca_graph = self.lca_graph()
        self.analyzed_graf = self.lca_graph_analysis()

    def merge_two_dicts(x, y):
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z

    def name2taxid(self,name):
        species_name = self._seqid[name.split(" ")[0]]
        return species_name

    def taxid_name_translator(self,taxid):
        return  str(Taxonomic_assignment.ncbi.get_taxid_translator([taxid])[taxid])

    def species_list(self):

        " Returns list of species proper for NCBI taxonomy for species "

        species_taxid_dict = self._seqid
        species_list = []

        for key in species_taxid_dict:
            try:
                species_list.append(self.taxid_name_translator(int(species_taxid_dict[key])))
            except KeyError:
                print key
        return species_list

    #NEW - Haven't tested yet
    def species_identifier(self,input_name):
        " Returns true if given name is in list of species "
        species_name_list = self.species_list()
        if input_name in species_name_list:

            return True
        else:
            return False


    def dict_init(self,dictionary):

        "Creates data structure like this in LTU: {1261581: {'2006943': 0, '1261582': 0, '2006970': 0, '2006944': 0}, developed for The closest db rel"

        new_dict = {}

        for key in dictionary:
            if dictionary[key] not in new_dict:
                new_dict[dictionary[key]] = {key:0}
            else:
                new_dict[dictionary[key]][key] = 0
        return new_dict


    def dict_translator_taxid_name(self,dict):
        #nie infromatywne
        new_dict = {}
        for key in dict:
            new_dict[self.taxid_name_translator(key)] = {}
            for key_sec in dict[key]:
                new_dict[self.taxid_name_translator(key)][self.taxid_name_translator(int(key_sec))] = dict[key][key_sec]
        return new_dict

    def get_NCBI_tree(self,covered_):
        ncbi = self.ncbi
        tree = ncbi.get_topology(covered_,intermediate_nodes=True)
        return tree

    def seqid2taxid(self,seqidmap):

        seqid2taxid = {}
        with open(seqidmap,"r") as f:
            for line in f:
                splited = line.split("\t")
                tmp = splited[1].strip("\n")
                seqid2taxid[splited[0]] = tmp.strip(" ")
        return seqid2taxid

    def _seqid_inverted(self,_seqid):
        inv_map = {v: k for k, v in _seqid.iteritems()}
        return inv_map

    def get_lineage(self,taxid):
        ncbi = self.ncbi

        lineage =  ncbi.get_lineage(taxid)

        names = ncbi.get_taxid_translator(lineage)
        lin = [names[taxid] for taxid in lineage]
        return lin

    def get_lineage_without_tree(self,taxid,ncbi_tree):
        ncbi = ncbi_tree
        lineage =  ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        lin = [names[taxid] for taxid in lineage]
        return lin

    def coverage_file(self,coverage_file):
        covered_ = []
        almost_full = []
        with open(coverage_file,"r") as f:
            for i in f:
                if i.split(" ")[0] != "Genome":
                    coverage =  float((i.split(" ")[-1]).split(",")[1])
                    species_name = self._seqid[i.split(" ")[0]]

                    if coverage <= 90.0:
                        covered_.append(species_name)
                    if coverage > 90.0 :
                        almost_full.append(species_name)
        return covered_,almost_full

    def lca_assignment(self):

        cov,almost_full =  self.coverage_file("cov_list.txt")
        tree = self.get_NCBI_tree(cov)
        dict_of_assigned =  {}
        ncbi = Taxonomic_assignment.ncbi
        for i in itertools.permutations(cov,2):

            lca = tree.get_common_ancestor(i[0],i[1])

            if self.get_lineage_without_tree(i[0],ncbi)[-1] not in dict_of_assigned:
                dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]] = lca.taxid
            else:
                if len(self.get_lineage_without_tree(dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]],ncbi))<  len(self.get_lineage_without_tree(lca.taxid,ncbi)):
                    dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]] = lca.taxid
        return dict_of_assigned,almost_full

    def lca_assignment_just_taxids(self):

        cov,almost_full =  self.coverage_file("cov_list.txt")

        if len(cov) >= 1:
            tree = self.get_NCBI_tree(cov)
            dict_of_assigned =  {}
            ncbi = Taxonomic_assignment.ncbi

            if len(cov) == 1:

                print "It seems that only one cpGenome were detected but with low coverage. Those reads were assigned as Eukaryote"
                dict_of_assigned[cov[0]] = 2759
                return dict_of_assigned,almost_full

            else:

                for i in itertools.permutations(cov,2):

                    lca = tree.get_common_ancestor(i[0],i[1])
                    if i[0] not in dict_of_assigned:
                        dict_of_assigned[i[0]] = lca.taxid
                    else:
                        if len(self.get_lineage_without_tree(dict_of_assigned[i[0]],ncbi)) <  len(self.get_lineage_without_tree(lca.taxid,ncbi)):
                            dict_of_assigned[i[0]] = lca.taxid

                return dict_of_assigned,almost_full
        else:
            dict_of_assigned = {}
            return dict_of_assigned,almost_full

    def _LTU_calc_with_taxa_name(self):

        unique = {}
        lca_assign,almost_full = lca_assignment(self._seqid)
        for key in lca_assign:
            if lca_assign[key] not in unique:
                unique[lca_assign[key]] = [key]
            else:
                unique[lca_assign[key]].append(key)

        return unique,almost_full

    #Ta funkcja do wymianki !!! - BUG!!
    def transform_read_name(self,readname):

        if len(readname.split(" ")) > 1:
            #SRA Reads
            readname = readname.replace(" ","_")
            return readname
#            split_1 = readname.split(".")
#
#            if len(split_1) > 1:
#
#                split_2 = split_1[2].split(" ")
#                return  "%s.%s.%s" % (split_1[0],split_1[1],split_2[0])
#            else:
#                split_2 = readname.split(" ")
#                return  "%s.%s" % (split_2[0],split_2[1])
        else:
            return readname

    def sam_parse(self):

        if self.sam_type  == True:

            samfile =  pysam.AlignmentFile(self.input)
        else:
            _samfile = glob.glob("*_final_mapped.sam")[0]
            samfile = pysam.AlignmentFile(_samfile)
        dict_reads = {}

        for read in samfile.fetch():
            try:
                if read.is_unmapped:
                    pass
                else:

                    #here we have to check -- !
                    #if read.qname == "ERR538178.142953062.1 H3:C2FGHACXX:1:2310:4346:40753":

                    new_read_name =  self.transform_read_name(read.qname)


                    if new_read_name[:-2] == ".1":
                        record = new_read_name[:-2]
                    else:
                        record = new_read_name

                    if new_read_name[-2:] == ".2":

                        pass

                    if record not in dict_reads:
                        dict_reads[record] = ([],[])

                    if read.is_read1:

                        if read.is_secondary:
                            if new_read_name[-2:] == ".2":
                                pass
                            else:
                                dict_reads[record][0].append(read)

                        else:
                            dict_reads[record][0].append(read)
                    if read.is_read2:

                        if read.is_secondary:

                            if new_read_name[-2:] == ".2":
                                pass
                            else:
                                dict_reads[record][1].append(read)
                        else:
                            dict_reads[record][1].append(read)

            except ValueError:
                pass

        return dict_reads

    def reads_processing_easy(self):

        ncbi = self.ncbi
        for record in self.almost_full:

            self.lca_assign[str(record)] = int(record)


        for record in self._dict_of_reads:

            r1 = self._dict_of_reads[record][0]

            r2 = self._dict_of_reads[record][1]

            try:


                try:
                    if len(r1) != 0 and len(r2) != 0:
                        name_r1 = str(self.name2taxid(r1[0].reference_name))
                        name_r2 = str(self.name2taxid(r2[0].reference_name))

                        yield "%s\t%s\n" % (record+".1", ";".join(self.get_lineage(self.lca_assign[name_r1])))
                        yield "%s\t%s\n" % (record+".2", ";".join(self.get_lineage(self.lca_assign[name_r2])))

                    elif len(r1) != 0 and len(r2) == 0:
                        name_r1 = str(self.name2taxid(r1[0].reference_name))
                        yield "%s\t%s\n" % (record+".1", ";".join(self.get_lineage(self.lca_assign[name_r1])))

                    elif len(r1) == 0 and len(r2) != 0:
                        name_r2 = str(self.name2taxid(r2[0].reference_name))
                        yield "%s\t%s\n" % (record+".2", ";".join(self.get_lineage(self.lca_assign[name_r2])))

                except ValueError:
                    print "Value error "+ name_r1
            except KeyError:
                print "Key error "+ name_r1

    def processing_cdr_data_stucture(self):


        for record in self.almost_full:

            self.lca_assign[str(record)] = int(record)

        the_cdr_data_structure =  self.dict_init(self.lca_assign)

        for record in self._dict_of_reads:

            r1 = self._dict_of_reads[record][0]
            r2 = self._dict_of_reads[record][1]

            try:
                try:
                    if len(r1) != 0 and len(r2) != 0 :
                        name_r1 = self.name2taxid(r1[0].reference_name)
                        name_r2 = self.name2taxid(r2[0].reference_name)
                        the_cdr_data_structure[self.lca_assign[name_r1]][name_r1] += 1
                        the_cdr_data_structure[self.lca_assign[name_r2]][name_r1] += 1

                    elif len(r1) != 0 and len(r2) == 0:
                        name_r1 = self.name2taxid(r1[0].reference_name)
                        the_cdr_data_structure[self.lca_assign[name_r1]][name_r1] += 1

                    elif len(r1) == 0 and len(r2) != 0:
                        name_r2 = self.name2taxid(r2[0].reference_name)
                        the_cdr_data_structure[self.lca_assign[name_r2]][name_r2] += 1

                except ValueError:
                    pass
                    #print "Val err "+ r1
            except KeyError:
                pass
        return the_cdr_data_structure

    def process_taxonomic_assignment_to_file_easy(self):
        to_write = list(self.reads_processing_easy())
        #It's better to call function 'write' once than many sam_analyzer times
        open(self._sample_name+"_assignment.txt", "wb").write(''.join(to_write))


    def get_count(self):

        "Takes taxonomic assignment file as an input"

        assignment_file = glob.glob("*_assignment.txt")
        if len(assignment_file) == 0:
            logger.error("There is no taxonomic assignment file, but should be here!")
            sys.exit()
        num_lines = sum(1 for line in open(assignment_file[0]))
        if num_lines == 0:
            logger.warn("There was no classified chloroplast reads, re-run with diffrent settings")
            sys.exit()
        return num_lines

    def open_lca_file(self):
        "Prepares file for lca graph"
        count = self.get_count()
        assignment_file = glob.glob("*_assignment.txt")
        taxa_count_dict = {}
        with open(assignment_file[0],"r") as f:
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

    def lca_graph(self):

        G = nx.DiGraph()

        dict_ltu = self.dict_translator_taxid_name(self.processing_cdr_data_stucture())

        taxa_count_dict, longest =  self.open_lca_file()

        good_taxa = {}
        bad_taxa = {}

        G.add_node('root',count = self.get_count())

        for key in taxa_count_dict:


            lineage = key.split(";")
            #edges and nodes
            for i in xrange(len(lineage)-1):

                if lineage[i+1] not in G.nodes():
                    if lineage[i+1] in dict_ltu:
                        G.add_node(lineage[i+1],count=int(taxa_count_dict[key]),ltu_genomes=dict_ltu[lineage[i+1]])
                    else:
                        G.add_node(lineage[i+1],count=int(taxa_count_dict[key]),ltu_genomes="")

                else:

                    old = G.nodes[lineage[i+1]]['count']
                    G.nodes[lineage[i+1]]['count'] = old + taxa_count_dict[key]
                G.add_edge(lineage[i],lineage[i+1])

        return G

    "Relative min_supporting - threshold * SUM"

    def lca_graph_analysis(self):

        dict_ltu = self.dict_translator_taxid_name(self.processing_cdr_data_stucture())
        count = self._lca_graph.nodes['root']['count']
        species_treshold = int(count*self.treshold)
        #here can appear a logger info message
        print "Species treshold was set on "+str(species_treshold)
        #node pruning
        to_del = []
        #Tu moze jest zla konstrukcja
        #Moze byc - rozpocznij od ostatniego
        for node in self._lca_graph.nodes():
            if self._lca_graph.nodes[str(node)]['count'] < species_treshold:
                to_del.append(node)

        G =  self.lca_graph()
        G_tmp = self.lca_graph()

        sorted_ = sorted((self.lca_graph()).degree, key=lambda x: x[1])

        for node in sorted_:

                if node[0] in to_del:

                    pre = list(G_tmp.predecessors(node[0]))[0]
                    if pre in dict_ltu:
                         old = G.nodes[pre]['ltu_genomes']
                         merged = merge_two_dicts(old, G.nodes[node[0]]["ltu_genomes"])
                         G.nodes[pre]["ltu_genomes"] = merged
                    else:
                        G.nodes[pre]["ltu_genomes"] =  G_tmp.nodes[node[0]]["ltu_genomes"]

        G.remove_nodes_from(to_del)

        return G

    def species_level(self):

        " Returns dictonary of counted reads on the species level  -- tu blad"

        #species predefined dict
        #To polaczyc z funkcja sprawdzajaca czy dany takson to gatunek z listy
        species = {}
        for node in self.analyzed_graf.nodes():
            if self.analyzed_graf.degree[node] == 1 and len(node.split(" ")) == 2:
                species[node] = self.analyzed_graf.nodes[node]['count']
        return species

    def ltu_every_level(self):

        " Returns dictionary with ltu and genomes related with - The closest db relative "

        species = {}
        for node in self.analyzed_graf.nodes():
            if self.analyzed_graf.degree[node] == 1 and node != "root":
                species[node] = (self.analyzed_graf.nodes[node]['count'],self.analyzed_graf.nodes[node]['ltu_genomes'])
        return species

    def ltu_catch(self):

        ltu_after_min_support = []
        for node in self.analyzed_graf.nodes():
            if self.analyzed_graf.degree[node] == 1 and node != "root":
                ltu_after_min_support.append((nx.shortest_path(self.analyzed_graf, "root", target=node),self.analyzed_graf.nodes[node]['count']))
        return ltu_after_min_support

    def pandas_data_frame_species_level(self):
        _species_level = self.species_level()
        df = pd.DataFrame(_species_level.items(),columns=['Taxon', 'Counts'])
        df.to_csv(self._sample_name+"_species_level_table.csv")

    def krona_file_preparing(self,project_name):

        catched_taxa = self.ltu_catch()

        with open(project_name+"_krona.txt","w") as f:
            for taxa in catched_taxa:
                tax_path = "\t".join(taxa[0][2:])
                f.write(str(taxa[1])+"\t"+tax_path+"\n")

        cmd = "ktImportText %s -o %s" % (project_name+"_krona.txt",project_name+"_krona.html")
        subprocess.check_call(cmd,shell=True)

    def sigle_bam_output(self):

        "Splits SAM file into smaller which containes only well-covered genomes - tutaj trzeba uzyc nazwy -> nie taxidu!!!!"

        covered_genomes  = self.almost_full
        sam_splitter(self.input, self.sam_type, covered_genomes)


class Taxonomic_assignment_Runner:

    def __init__ (self,input,output,settings):


        self.input = input
        self.output = output
        self.settings = settings
        self.sam_type = True


        if (self.input).split(".")[1] == "sam":

            self.sam_type = True

        else:

            self.sam_type = False
        try:

            self.seqidmap = Settings_loader_yaml(self.settings).yaml_handler()["Databases and mapping files"]["seqid2taxid.map"]
            self.lca_treshold = Settings_loader_yaml(self.settings).yaml_handler()["Params"]["lca_treshold"]
            self.project_name = output

        except KeyError:

            logger.error("Something went wrong during loading files, check paths")
            sys.exit()

    def process(self):

        starting_dir = os.getcwd()
        project_name = self.project_name
        grlog = []

        logger.info("Procesing "+project_name)

        dir = "%s/" % (project_name)

        if os.path.isdir(dir):
            os.chdir(dir)

        else:
            os.mkdir(dir)
            os.chdir(dir)

        c = Coverage('bincov.txt',project_name+"_chloroplasts.hitstats",self.settings)
        c.getpercentage_cov()
        c.report_cov()
        lca_postprocess =  Taxonomic_assignment(self.project_name,self.lca_treshold,self.seqidmap,self.sam_type,self.input)
#            lca_postprocess.process_taxonomic_assignment_to_file_easy()
        lca_postprocess.pandas_data_frame_species_level()
        lca_postprocess.krona_file_preparing(self.project_name)
        lca_postprocess.sigle_bam_output()
        os.chdir(starting_dir)
