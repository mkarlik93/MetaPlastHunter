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

from ete3 import  NCBITaxa
import networkx as nx
import pysam
import numpy
import glob
import itertools

import pandas as pd
from settings import *
from cov import Coverage
from cov import Coverage_utilities

logger = logging.getLogger("src.taxonomic_assignment")
logging.basicConfig(level=logging.INFO)

def generates_ani(target_species_genome):

    """ Returns identity level of reads of covered genomes at species level min at species level """

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


class Taxonomic_assignment(object):


    ncbi = NCBITaxa()

    def __init__ (self,sample,threshold,seqidmap,sam_type,input,settings,avg_coverage_dict):

        """
        Input Parameters
        ----------
        sam_type : boolean
            bincov.txt, draft coverage file, produced during mapping

        input: str
            absolute path of input file

        threshold : int
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
            Prunned _lca_graph that keeps only taxonomic positions that have been above the threshold (min support value)

        """
        self.settings = settings
        self.sam_type = sam_type
        self.input = input
        self._seqid = self.seqid2taxid(seqidmap)
        self._sample_name = sample
        self._seqid_inverted = self._seqid_inverted(self._seqid)
        self.threshold = threshold
        self.lca_assign,self.almost_full = self.lca_assignment_just_taxids()
        self.taxa_count_dict, self.dict_ltu, self.count = self.sam_parse_experimental()
        self._lca_graph = self.lca_graph()
        self.analyzed_graf = self.lca_graph_analysis()
        self.absolute_coverage_loaded = self.absolute_coverage_loading_to_NCBI()
        self.empirical_coverage_ncbi = self.empirical_coverage_to_NCBI()
        self.avg_coverage_dict = avg_coverage_dict

    def merge_two_dicts(self, x, y):

        """ Merges two dictionaries """

        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z

    def name2taxid(self,name):

        """ Transforms query name to taxid """

        species_name = self._seqid[name.split(" ")[0]]
        return species_name

    def taxid_name_translator(self,taxid):
        return  str(Taxonomic_assignment.ncbi.get_taxid_translator([taxid])[taxid])

    def species_list(self):

        """ Returns list of species proper for NCBI taxonomy for species """

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
        """ Returns true if given name is in list of species """
        species_name_list = self.species_list()
        if input_name in species_name_list:

            return True
        else:
            return False

    def dict_init(self,dictionary):

        """Creates data structure like this in

        LTU: {1261581: {'2006943': 0, '1261582': 0, '2006970': 0, '2006944': 0}

         developed for The closest db relative"""

        new_dict = {}

        for key in dictionary:
            if dictionary[key] not in new_dict:

                new_dict[dictionary[key]] = {key:0}

            else:
                new_dict[dictionary[key]][key] = 0

        return new_dict


    def dict_translator_taxid_name(self,dict):

        """ Transforms dict of taxid names  """
        #nie infromatywne

        new_dict = {}
        for key in dict:
            new_dict[self.taxid_name_translator(key)] = {}

            for key_sec in dict[key]:
                new_dict[self.taxid_name_translator(key)][self.taxid_name_translator(int(key_sec))] = dict[key][key_sec]
        return new_dict

    def get_NCBI_tree(self,covered_):

        """Initialization of NCBI tree - from ete3 packagess"""

        ncbi = self.ncbi
        tree = ncbi.get_topology(covered_,intermediate_nodes=True)
        return tree

    def seqid2taxid(self,seqidmap):

        """Creates dictionary instance with seqid - taxid structure"""

        seqid2taxid = {}
        with open(seqidmap,"r") as f:
            for line in f:
                splited = line.split("\t")
                tmp = splited[1].strip("\n")
                seqid2taxid[splited[0]] = tmp.strip(" ")
        return seqid2taxid

    def _seqid_inverted(self,_seqid):

        """ Returns inverted dictionary """

        inv_map = {v: k for k, v in _seqid.iteritems()}
        return inv_map

    def get_lineage(self,taxid):

        """Generates taxonomical lineage based on taxid"""

        ncbi = self.ncbi
        lineage =  ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        lin = [names[taxid] for taxid in lineage]

        return lin

    def get_lineage_without_tree(self,taxid,ncbi_tree):

        """ Returns lineage (based Tncbi) of given taxid """

        ncbi = ncbi_tree
        lineage =  ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        lin = [names[taxid] for taxid in lineage]
        return lin

    def absolute_coverage_loading_to_NCBI(self):

        """Loads absolute coverage to dict of NCBI taxa names """

        with open("cov_list.txt","r") as f:

            return dict([(self.taxid_name_translator(int(self._seqid[line.split(",")[0]])),line.split(",")[1].strip("\n")) for line in f])

    def empirical_coverage_to_NCBI(self):

        """Loads empirical coverage to dict of NCBI taxa names """

        empirical_threshold = Coverage_utilities(self.settings).load_emprirical_coverage()
        return dict([(self.taxid_name_translator(int(self._seqid[key])),float(empirical_threshold[key].strip("\n"))) for key in empirical_threshold])

    def mean_coverage_to_NCBI(self):

        """Loads mean coverage to dict with NCBI taxa names """

        return dict([(self.taxid_name_translator(int(self._seqid[key])),float(self.avg_coverage_dict[key])) for key in self.avg_coverage_dict])


    def coverage_file(self,coverage_file):

        """ Handles coverage file and seperates
        genomes into two classes LTU and TTU """

        #For every genome MPH uses pre-calculated empirical threshold
        empirical_threshold = Coverage_utilities(self.settings).load_emprirical_coverage()
        covered_ = []
        almost_full = []

        with open(coverage_file,"r") as file:

            for coverage_record in file:

                """ Split coverage file which looks like this: record,coverage  """
                splited_list = coverage_record.split(",")
                coverage =  float(splited_list[1])
                species_name = self._seqid[splited_list[0]]

                if coverage <= float(empirical_threshold[splited_list[0]]):
                    #covered_.append(splited_list[0])
                    covered_.append(species_name)

                if coverage > float(empirical_threshold[splited_list[0]]):
                    #almost_full.append(splited_list[0])
                    almost_full.append(species_name)

        return covered_,almost_full

    def lca_assignment(self):

        """  Finds the longest Lowest Common Ancestor

        between to covered reference genome.  """

        cov,almost_full =  self.coverage_file("cov_list.txt")
        tree = self.get_NCBI_tree(cov)
        dict_of_assigned =  {}
        ncbi = Taxonomic_assignment.ncbi
        for i in itertools.permutations(cov,2):

            #MPH conducts LCA assignmnet of reads usig built-in function from ete3
            lca = tree.get_common_ancestor(i[0],i[1])

            if self.get_lineage_without_tree(i[0],ncbi)[-1] not in dict_of_assigned:
                dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]] = lca.taxid
            else:
                if len(self.get_lineage_without_tree(dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]],ncbi))<  len(self.get_lineage_without_tree(lca.taxid,ncbi)):
                    dict_of_assigned[self.get_lineage_without_tree(i[0],ncbi)[-1]] = lca.taxid

        return dict_of_assigned,almost_full

    def lca_assignment_just_taxids(self):

        """LCA assignmnet workflow"""

        cov,almost_full =  self.coverage_file("cov_list.txt")

        if len(cov) >= 1:
            tree = self.get_NCBI_tree(cov)
            dict_of_assigned =  {}
            ncbi = Taxonomic_assignment.ncbi

            if len(cov) == 1:

                logger.info("It seems that only one cpGenome were detected but with low coverage. Those reads were assigned as Eukaryote")
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

    def _ltu_calc_with_taxa_name(self):

        unique = {}
        lca_assign, almost_full = lca_assignment(self._seqid)

        for key in lca_assign:

            if lca_assign[key] not in unique:
                unique[lca_assign[key]] = [key]
            else:
                unique[lca_assign[key]].append(key)

        return unique, almost_full

    def transform_read_name(self,readname):

        if len(readname.split(" ")) > 1:

            readname = readname.replace(" ","_")
            return readname
        else:
            return readname

    def sam_parse_experimental(self):

        tmp_dict_of_reads = dict([(reference_genome,0) for reference_genome in self._seqid.keys()])


        for record in self.almost_full:

            self.lca_assign[str(record)] = int(record)

        the_cdr_data_structure =  self.dict_init(self.lca_assign)

        count = 0

        if self.sam_type  == True:
            samfile =  pysam.AlignmentFile(self.input)

        else:
            samfile = pysam.AlignmentFile(glob.glob("*_final_mapped.sam")[0])

        for read in samfile.fetch():
            try:
                if read.is_unmapped:
                    pass
                else:
                    if read.reference_name in tmp_dict_of_reads:
                        count = count + 1
                        tmp_dict_of_reads[read.reference_name] = tmp_dict_of_reads[read.reference_name] + 1
                    else:
                        pass

            except ValueError:
                pass
        #dict with read counts

        assigned_reads = dict([(self._seqid[reference_genome],tmp_dict_of_reads[reference_genome]) for reference_genome in tmp_dict_of_reads])
        assigned_reads = dict((k, v) for k, v in assigned_reads.items() if v > 0)

        for taxa in the_cdr_data_structure:

            for reference in the_cdr_data_structure[taxa].keys():
                the_cdr_data_structure[taxa][reference] = assigned_reads[reference]

        taxa_count_dict = dict((str(";".join(self.get_lineage(taxa))),sum(v.values())) for taxa, v in the_cdr_data_structure.items())

        dict_ltu = self.dict_translator_taxid_name(the_cdr_data_structure)

        #Checking reads count

        if count == 0:
            sys.exit("There was no preclassified plastid reads!")

        return taxa_count_dict, dict_ltu, count

    def lca_graph(self):

        "Creates graph representation of taxonomical assigned reads"

        G = nx.DiGraph()
        G.add_node('root',count = self.count)

        for key in self.taxa_count_dict:

            lineage = key.split(";")
            #edges and nodes
            for i in range(len(lineage)-1):

                if lineage[i+1] not in G.nodes():
                    if lineage[i+1] in self.dict_ltu:
                        G.add_node(lineage[i+1],count=int(self.taxa_count_dict[key]),ltu_genomes=self.dict_ltu[lineage[i+1]])

                    else:

                        G.add_node(lineage[i+1],count=int(self.taxa_count_dict[key]),ltu_genomes="")

                else:

                    old = G.nodes[lineage[i+1]]['count']
                    G.nodes[lineage[i+1]]['count'] = old + self.taxa_count_dict[key]
                G.add_edge(lineage[i],lineage[i+1])

        return G

    def lca_graph_analysis(self):

        """ Prunnes graph branches until branch is greater than threshold  """


        #Relative min_supporting - threshold * SUM
        sum = self._lca_graph.nodes['root']['count']
        species_threshold = int(sum*self.threshold)
        logger.info("Species threshold was set on "+str(species_threshold))
        to_del = []

        for node in self._lca_graph.nodes():
            if self._lca_graph.nodes[str(node)]['count'] < species_threshold:
                to_del.append(node)

        G =  self.lca_graph()
        G_tmp = self.lca_graph()

        sorted_ = sorted((self.lca_graph()).degree, key=lambda x: x[1])

        for node in sorted_:

                if node[0] in to_del:

                    pre = list(G_tmp.predecessors(node[0]))[0]

                    if pre in self.dict_ltu:

                         old = G.nodes[pre]['ltu_genomes']

                         merged = self.merge_two_dicts(old, G.nodes[node[0]]["ltu_genomes"])

                         G.nodes[pre]["ltu_genomes"] = merged

                    else:

                         G.nodes[pre]["ltu_genomes"] =  G_tmp.nodes[node[0]]["ltu_genomes"]

        G.remove_nodes_from(to_del)

        return G

    def species_level(self):

        " Returns dictonary of counted reads on the species level  -- tu blad"

        #species predefined dict
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

    def is_LTU(node):

        pass

    def ltu_empirical_threshold_comparsion(self):

        taxonomic_assignment_graph = self.analyzed_graf
        good_taxa  = {}
        empirical_threshold = self.empirical_coverage_to_NCBI()
        avg_coverage = self.mean_coverage_to_NCBI()

        relaxation_parameter = 0.9

        for node in taxonomic_assignment_graph.nodes():

            if 'ltu_genomes' in taxonomic_assignment_graph.nodes[node]:

                if type(taxonomic_assignment_graph.nodes[node]['ltu_genomes']) is dict:

                    try:

                        ltu_sum = sum([float(self.absolute_coverage_loaded[x]) for x in taxonomic_assignment_graph.nodes[node]['ltu_genomes'].keys()])


                        avg_empirical_threshold = relaxation_parameter *  numpy.mean([float(empirical_threshold[x]) for x in taxonomic_assignment_graph.nodes[node]['ltu_genomes'].keys()])

                        if ltu_sum >= avg_empirical_threshold:
                            good_taxa[node] = ltu_sum

                        else:
                            pass
                    except:

                        print "Error  " + taxonomic_assignment_graph.degree[node]

            else:
                print taxonomic_assignment_graph.nodes[node]


        return good_taxa

    def ltu_catch(self):

        ltu_after_min_support = []
        for node in self.analyzed_graf.nodes():
            if self.analyzed_graf.degree[node] == 1 and node != "root":

                ltu_after_min_support.append((nx.shortest_path(self.analyzed_graf, "root", target=node),self.analyzed_graf.nodes[node]['count']))

        return ltu_after_min_support

    def ltu_catch_prunned(self):

        probable_taxa = self.ltu_empirical_threshold_comparsion().keys()
        average_coverage = self.mean_coverage_to_NCBI()
        ltu_after_min_support = []

        for node in self.analyzed_graf.nodes():

            if 'ltu_genomes' in self.analyzed_graf.nodes[node]:

                if type(self.analyzed_graf.nodes[node]['ltu_genomes']) is dict and node in probable_taxa:

                    try:
                        keys_for_average = self.analyzed_graf.nodes[node]['ltu_genomes'].keys()
                        avg_sum = sum(list([average_coverage[key] for key in keys_for_average]))
                        ltu_after_min_support.append((nx.shortest_path(self.analyzed_graf, "root", target=node),self.analyzed_graf.nodes[node]['count'],avg_sum))
                    except TypeError:
                        print "To ten error "+ node

        return ltu_after_min_support

    def pandas_data_frame_species_level(self):
        _species_level = self.species_level()
        df = pd.DataFrame(_species_level.items(),columns=['Taxon', 'Counts'])
        df.to_csv(self._sample_name+"_species_level_table.csv")

    def krona_prunned_count_preparing(self,project_name):

        catched_taxa = self.ltu_catch_prunned()

        with open(project_name+"_count_krona_.txt","w") as f:
            for taxa in catched_taxa:
                tax_path = "\t".join(taxa[0][2:])
                f.write(str(taxa[1])+"\t"+tax_path+"\n")

        cmd = "ktImportText %s -o %s" % (project_name+"_count_krona_.txt",project_name+"_prunned_count_krona.html")

        subprocess.check_call(cmd,shell=True)


    def krona_prunned_avg_coverage(self, project_name):

        catched_taxa = self.ltu_catch_prunned()
        with open(project_name+"_average_coverage_krona.txt","w") as f:
            for taxa in catched_taxa:
                tax_path = "\t".join(taxa[0][2:])
                f.write(str(taxa[2])+"\t"+tax_path+"\n")

        cmd = "ktImportText %s -o %s" % (project_name+"_average_coverage_krona.txt",project_name+"_average_coverage_krona.html")

        subprocess.check_call(cmd,shell=True)


    def krona_file_preparing(self,project_name):

        """ Prepares file formatted for Krona tools """

        catched_taxa = self.ltu_catch()

        with open(project_name+"_krona.txt","w") as f:
            for taxa in catched_taxa:
                tax_path = "\t".join(taxa[0][2:])
                f.write(str(taxa[1])+"\t"+tax_path+"\n")
        #Tu cos by sie przydalo - settings
        cmd = "ktImportText %s -o %s" % (project_name+"_krona.txt",project_name+"_krona.html")
        subprocess.check_call(cmd,shell=True)


class Cleaner:


    def __init__ (self,output):

        self.output = output


    def clean(self):

        os.remove(self.output+"_de_complex_R1.fastq")
        os.remove(self.output+"_de_complex_R2.fastq")
        os.remove(self.output+"_filtered_chloroplasts_reads_R1.fq")
        os.remove(self.output+"_filtered_chloroplasts_reads_R2.fq")

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
            self.lca_threshold = Settings_loader_yaml(self.settings).yaml_handler()["Params"]["lca_threshold"]
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


        #coverage instance
        c = Coverage('bincov.txt',project_name+"_chloroplasts.hitstats",self.settings)
        c.getpercentage_cov()
        c.report_cov()
        avg_coverage = c.average_coverage()
        lca_postprocess = Taxonomic_assignment(self.project_name,self.lca_threshold,self.seqidmap,self.sam_type,self.input,self.settings,avg_coverage)
        lca_postprocess.pandas_data_frame_species_level()
        lca_postprocess.krona_file_preparing(self.project_name)


        lca_postprocess.krona_prunned_count_preparing(self.project_name)

        lca_postprocess.krona_prunned_avg_coverage(self.project_name)


        if self.sam_type == False:

            Cleaner(project_name).clean()

        os.chdir(starting_dir)