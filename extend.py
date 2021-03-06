import igraph
from scipy.stats import binom
import os
import pickle
import pysam
import vcf

from SomaticHaplotype import *

################################################################################
# functions
################################################################################

def create_graph(extended_ps_dict):

  vertex_dict_ps_key = {}
  vertex_count = -1 # graph vertices start from 0

  for ps1 in sorted(extended_ps_dict.keys()):
    if ps1 + "_s1" not in vertex_dict_ps_key:
      vertex_count += 1
      vertex_dict_ps_key[ps1 + "_s1"] = vertex_count
    for ps2 in sorted(extended_ps_dict[ps1].keys()):
      if ps2 + "_s2" not in vertex_dict_ps_key:
        vertex_count += 1
        vertex_dict_ps_key[ps2 + "_s2"] = vertex_count

  vertex_dict_number_key = {v:k for k,v in vertex_dict_ps_key.items()}

  number_of_vertices = len(vertex_dict_number_key.keys())

  phase_set_graph = igraph.Graph()
  phase_set_graph.add_vertices(number_of_vertices)

  edge_tuple_list = []
  edge_weight_list = []

  for ps1 in extended_ps_dict.keys():
    for ps2 in extended_ps_dict[ps1].keys():
      edge_weight = extended_ps_dict[ps1][ps2][-1] # final element of this list
      if edge_weight == 0: # no recommendation, do not form edge
        next
      else:
        edge_tuple_list.append((vertex_dict_ps_key[ps1 + "_s1"], vertex_dict_ps_key[ps2 + "_s2"]))
        edge_weight_list.append(edge_weight)

  phase_set_graph.add_edges(edge_tuple_list)
  phase_set_graph.es["weight"] = edge_weight_list

  return(phase_set_graph)

def extend_phase_sets(ps_dict1, ps_dict2, chrom, start, end):

  # value that will be returned:
  extended_ps_dict1 = {}

  # for each phase set in ps1 check if ps overlaps range
  ps1_phase_sets_in_range = []
  for ps in [ x for x in list(ps_dict1['phase_sets'].keys()) if x.startswith(chrom + ":") ]:
    phase_set_1 = ps_dict1['phase_sets'][ps]

    if phase_set_1.return_FirstVariantPosition() == "NA" or phase_set_1.return_LastVariantPosition() == "NA":
      next
    elif ranges_overlap(phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), chrom, start, end): # ps range overlaps given range
      ps1_phase_sets_in_range.append(ps)

  # for each phase set in ps2 check if ps overlaps range
  ps2_phase_sets_in_range = []
  for ps in [ x for x in list(ps_dict2['phase_sets'].keys()) if x.startswith(chrom + ":") ]:
    phase_set_2 = ps_dict2['phase_sets'][ps]

    if phase_set_2.return_FirstVariantPosition() == "NA" or phase_set_2.return_LastVariantPosition() == "NA":
      next
    elif ranges_overlap(phase_set_2.return_Chromosome(), phase_set_2.return_FirstVariantPosition(), phase_set_2.return_LastVariantPosition(), chrom, start, end): # ps range overlaps given range
      ps2_phase_sets_in_range.append(ps)

  # for each ps1 phase set in range, iterate over each in-range ps in ps2
  for ps1 in ps1_phase_sets_in_range:
    extended_ps_dict1[ps1] = {}
    phase_set_1 = ps_dict1['phase_sets'][ps1]
    for ps2 in ps2_phase_sets_in_range:
      phase_set_2 = ps_dict2['phase_sets'][ps2]

      # check if they overlap
      if ranges_overlap(phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), phase_set_2.return_Chromosome(), phase_set_2.return_FirstVariantPosition(), phase_set_2.return_LastVariantPosition()):

        min_overlap_position, max_overlap_position, length_overlap = ranges_overlap_stats(phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), phase_set_2.return_Chromosome(), phase_set_2.return_FirstVariantPosition(), phase_set_2.return_LastVariantPosition())

        # get list of overlapping variants
        overlapping_variants_list = list(set(phase_set_1.return_Variants().keys()).intersection(phase_set_2.return_Variants().keys()))
        n_variants_overlap = len(overlapping_variants_list)
        n_variants_flip_to_match = 0
        pct_switch = "NA"

        if n_variants_overlap > 0:
          for variant in overlapping_variants_list:
            if phase_set_1.return_Variants()[variant][0].return_Genotype() == phase_set_2.return_Variants()[variant][0].return_Genotype()[::-1]:
              n_variants_flip_to_match += 1
          pct_switch = n_variants_flip_to_match/n_variants_overlap

          # define null hypothesis using switch error rate
          # Switch error rate defined at ~2x10-4
          # Marks, P. et al.Resolving the full spectrum of human genome variation using Linked-Reads.Genome Res.29, 635–645 (2019).
          # PMID: 30894395
          # https://www.ncbi.nlm.nih.gov/pubmed/30894395
          # We set an ultra-conservative (random) short switch error rate at 0.5
          # Even using very conservative short switch error rate of 0.001 was too lenient
          p_short_switch_error = 0.5

          # null is no switch, so probability of switch is short switch error rate
          # p value is right side of probability distribution, or 1 - left side
          # P(Y >= X | No switch) = 1 - P(Y < X | No switch)
          p_value_switch = 1 - binom.cdf(n_variants_flip_to_match - 1, n_variants_overlap, p_short_switch_error)
          # null is switch, so probability of no switch is 1 - short switch error rate
          # p value is left side of probability probability distribution
          # P(Y <= X | Switch)
          p_value_no_switch = binom.cdf(n_variants_flip_to_match, n_variants_overlap, 1 - p_short_switch_error)

          min_p_value = min(p_value_switch, p_value_no_switch)

          if p_value_switch < p_value_no_switch and min_p_value < 0.001 and pct_switch > 0.95:
            recommendation = "Switch"
            graph_weight = 1
          elif p_value_no_switch < p_value_switch and min_p_value < 0.001 and pct_switch < 0.05:
            recommendation = "No Switch"
            graph_weight = 2
          else:
            recommendation = "No Recommendation"
            graph_weight = 0

          extended_ps_dict1[ps1][ps2] = [ps1, phase_set_1.return_Chromosome(), phase_set_1.return_FirstVariantPosition(), phase_set_1.return_LastVariantPosition(), ps2, phase_set_2.return_Chromosome(), phase_set_2.return_FirstVariantPosition(), phase_set_2.return_LastVariantPosition(), min_overlap_position, max_overlap_position, length_overlap, n_variants_overlap, n_variants_flip_to_match, p_value_switch, p_value_no_switch, min_p_value, pct_switch, recommendation, graph_weight] # graph_weight must be final element of list to work with create_graph()
  
  return(extended_ps_dict1)

def super_phase_set_relationships(extended_ps_dict, phase_set_graph):

  graph_clusters = phase_set_graph.clusters()

  super_set_dict_index_key = {}
  super_set_dict_cluster_key = {}

  for cluster_id in range(len(graph_clusters)):
    if cluster_id not in super_set_dict_cluster_key:
      super_set_dict_cluster_key[cluster_id] = []
    for x in graph_clusters[cluster_id]:
      if x not in super_set_dict_index_key:
        super_set_dict_index_key[x] = cluster_id
      if x not in super_set_dict_cluster_key[cluster_id]:
        super_set_dict_cluster_key[cluster_id].append(x)

  vertex_dict_ps_key = {}
  vertex_count = -1 # graph vertices start from 0

  for ps1 in sorted(extended_ps_dict.keys()):
    if ps1 + "_s1" not in vertex_dict_ps_key:
      vertex_count += 1
      vertex_dict_ps_key[ps1 + "_s1"] = vertex_count
    for ps2 in sorted(extended_ps_dict[ps1].keys()):
      if ps2 + "_s2" not in vertex_dict_ps_key:
        vertex_count += 1
        vertex_dict_ps_key[ps2 + "_s2"] = vertex_count

  vertex_dict_number_key = {v:k for k,v in vertex_dict_ps_key.items()}

  phase_set_relationship_dict = {}

  for cluster, index_list in super_set_dict_cluster_key.items():
    base_phase_set_index = None
    for index in index_list:
      if len(index_list) == 1:
        phase_set_relationship_dict[vertex_dict_number_key[index][:-3]] = ["NA"]*4
      else:
        if vertex_dict_number_key[index].endswith("_s1"):
          if base_phase_set_index is None:
            base_phase_set_index = index
            base_phase_set = vertex_dict_number_key[index][:-3]
          phase_set_relationship_dict[vertex_dict_number_key[index][:-3]] = [str(x) for x in [cluster, base_phase_set]]
          phase_set_relationship_dict[vertex_dict_number_key[index][:-3]].extend([str(x) for x in sum_graph_edges(phase_set_graph, index, base_phase_set_index)])

  return(phase_set_relationship_dict)

def ranges_overlap(chr1, start1, end1, chr2, start2, end2):

  if start1 is None: # if no start1, assume it is start2
    start1 = int(start2)
  if end1 is None: # if no end1, assume it is end2
    end1 = int(end2)
    
  if start2 is None: # if no start2, assume it is start1
    start2 = int(start1)
  if end2 is None: # if no end2, assume it is end1
    end2 = int(end1)

  if chr1 != chr2: # ranges on different chromosomes
    return(False)
  elif(end2 < start1 or start2 > end1):
    return(False)
  elif len(set(range(start1, end1 + 1)).intersection(range(start2, end2 + 1))) > 0:
    return(True)
  else:
    return(False)

def ranges_overlap_stats(chr1, start1, end1, chr2, start2, end2):
  overlap = set(range(start1, end1 + 1)).intersection(range(start2, end2 + 1))
  min_overlap_position = min(overlap)
  max_overlap_position = max(overlap)
  length_overlap = len(overlap)
  return([min_overlap_position, max_overlap_position, length_overlap])

def sum_graph_edges(phase_set_graph, vertex_number_1, vertex_number_2):
  all_shortest_paths = phase_set_graph.get_all_shortest_paths(v = vertex_number_1, to = vertex_number_2)

  n_shortest_paths = len(all_shortest_paths)

  if n_shortest_paths == 0:
    return("No connection")
  if n_shortest_paths != 1:
    sys.exit("Number of shortest paths != 1")

  weight_sum = 0
  n_edges = 0
  for x,y in zip(all_shortest_paths[0][:-1], all_shortest_paths[0][1:]):
    weight_sum += phase_set_graph[x,y]
    n_edges += 1
  
  # [0 if even (no switch) 1 if odd (switch), number of reference edges]
  return([weight_sum % 2, int(n_edges/2)])

################################################################################
# main
################################################################################   

def main(args):

  # open up summary file
  summary_file = open(args.sum, 'r')

  # path to input pickle file of first sample
  pickle_path_ps1 = args.ps1
  with open(pickle_path_ps1, 'rb') as pickle_file_ps1:
    bam_phase_set_dictionary_ps1 = pickle.load(pickle_file_ps1)
    vcf_variants_dictionary_ps1 = pickle.load(pickle_file_ps1)

  # path to input pickle file of second sample
  pickle_path_ps2 = args.ps2
  with open(pickle_path_ps2, 'rb') as pickle_file_ps2:
    bam_phase_set_dictionary_ps2 = pickle.load(pickle_file_ps2)
    vcf_variants_dictionary_ps2 = pickle.load(pickle_file_ps2)

  # parse the genomic range argument (which is required)
  chrom = str(args.range.split(":")[0])
  try:
    start = int(args.range.split(":")[1].split("-")[0])
  except:
    start = None
  try:
    end = int(args.range.split(":")[1].split("-")[1])
  except:
    end = None

  # run functions

  extended_phase_set_dictionary = extend_phase_sets(bam_phase_set_dictionary_ps1, bam_phase_set_dictionary_ps2, chrom, start, end)

  extended_phase_set_graph = create_graph(extended_phase_set_dictionary)

  # extended_pairwise_relationships = pairwise_phase_set_relationships(extended_phase_set_dictionary, extended_phase_set_graph)

  super_sets = super_phase_set_relationships(extended_phase_set_dictionary, extended_phase_set_graph)

  # write output for extended_phase_set_dictionary
  header_line = ["ps1", "ps1_Chromosome", "ps1_FirstVariantPosition", "ps1_LastVariantPosition", "ps2", "ps2_Chromosome", "ps2_FirstVariantPosition", "ps2_LastVariantPosition", "min_overlap_position", "max_overlap_position", "length_overlap", "n_variants_overlap", "n_variants_flip_to_match", "p_value_switch", "p_value_no_switch", "min_p_value", "pct_switch", "recommendation", "graph_weight"]
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extend_stats.tsv")
  output_file = open(output_file_path, 'w')
  output_file.write('\t'.join(header_line) + '\n')
  for k1,v1 in extended_phase_set_dictionary.items():
    for k2,v2 in extended_phase_set_dictionary[k1].items():
      if len(header_line) != len(v2):
        sys.exit("Length of header line and number of column not equal.")
      output_file.write('\t'.join([str(x) for x in v2]) + '\n')
  output_file.close()

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extended_phase_sets.tsv")
  output_file = open(output_file_path, 'w')

  header_row_list = summary_file.readline().strip().split()
  header_row_list.extend(["cluster", "cluster_base_phase_set", "switch_status", "n_reference_ps"])
  output_file.write('\t'.join(header_row_list) + '\n')
  for line in summary_file:
    ps_id = line.strip().split()[0]
    if ps_id in super_sets.keys():
      output_file.write('\t'.join(line.strip().split() + super_sets[ps_id]) + '\n')
    else:
      output_file.write('\t'.join(line.strip().split() + ["NA"]*4) + '\n')
  output_file.close()

  #output_file.close()

if __name__ == '__main__':
  main(args)
