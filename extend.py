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

def create_graph(extended_pb_dict):

  vertex_dict_pb_key = {}
  vertex_count = -1 # graph vertices start from 0

  for pb1 in sorted(extended_pb_dict.keys()):
    if pb1 + "_s1" not in vertex_dict_pb_key:
      vertex_count += 1
      vertex_dict_pb_key[pb1 + "_s1"] = vertex_count
    for pb2 in sorted(extended_pb_dict[pb1].keys()):
      if pb2 + "_s2" not in vertex_dict_pb_key:
        vertex_count += 1
        vertex_dict_pb_key[pb2 + "_s2"] = vertex_count

  vertex_dict_number_key = {v:k for k,v in vertex_dict_pb_key.items()}

  number_of_vertices = len(vertex_dict_number_key.keys())

  phase_block_graph = igraph.Graph()
  phase_block_graph.add_vertices(number_of_vertices)

  edge_tuple_list = []
  edge_weight_list = []

  for pb1 in extended_pb_dict.keys():
    for pb2 in extended_pb_dict[pb1].keys():
      edge_weight = extended_pb_dict[pb1][pb2][-1] # final element of this list
      if edge_weight == 0: # no recommendation, do not form edge
        next
      else:
        edge_tuple_list.append((vertex_dict_pb_key[pb1 + "_s1"], vertex_dict_pb_key[pb2 + "_s2"]))
        edge_weight_list.append(edge_weight)

  phase_block_graph.add_edges(edge_tuple_list)
  phase_block_graph.es["weight"] = edge_weight_list

  return(phase_block_graph)

def extend_phase_blocks(pb_dict1, pb_dict2, chrom, start, end):

  # value that will be returned:
  extended_pb_dict1 = {}

  # for each phase block in pb1 check if pb overlapb range
  pb1_phase_blocks_in_range = []
  for pb in [ x for x in list(pb_dict1['phase_blocks'].keys()) if x.startswith(chrom + ":") ]:
    phase_block_1 = pb_dict1['phase_blocks'][pb]

    if phase_block_1.return_FirstVariantPosition() == "NA" or phase_block_1.return_LastVariantPosition() == "NA":
      next
    elif ranges_overlap(phase_block_1.return_Chromosome(), phase_block_1.return_FirstVariantPosition(), phase_block_1.return_LastVariantPosition(), chrom, start, end): # pb range overlapb given range
      pb1_phase_blocks_in_range.append(pb)

  # for each phase block in pb2 check if pb overlapb range
  pb2_phase_blocks_in_range = []
  for pb in [ x for x in list(pb_dict2['phase_blocks'].keys()) if x.startswith(chrom + ":") ]:
    phase_block_2 = pb_dict2['phase_blocks'][pb]

    if phase_block_2.return_FirstVariantPosition() == "NA" or phase_block_2.return_LastVariantPosition() == "NA":
      next
    elif ranges_overlap(phase_block_2.return_Chromosome(), phase_block_2.return_FirstVariantPosition(), phase_block_2.return_LastVariantPosition(), chrom, start, end): # pb range overlapb given range
      pb2_phase_blocks_in_range.append(pb)

  # for each pb1 phase block in range, iterate over each in-range pb in pb2
  for pb1 in pb1_phase_blocks_in_range:
    extended_pb_dict1[pb1] = {}
    phase_block_1 = pb_dict1['phase_blocks'][pb1]
    for pb2 in pb2_phase_blocks_in_range:
      phase_block_2 = pb_dict2['phase_blocks'][pb2]

      # check if they overlap
      if ranges_overlap(phase_block_1.return_Chromosome(), phase_block_1.return_FirstVariantPosition(), phase_block_1.return_LastVariantPosition(), phase_block_2.return_Chromosome(), phase_block_2.return_FirstVariantPosition(), phase_block_2.return_LastVariantPosition()):

        min_overlap_position, max_overlap_position, length_overlap = ranges_overlap_stats(phase_block_1.return_Chromosome(), phase_block_1.return_FirstVariantPosition(), phase_block_1.return_LastVariantPosition(), phase_block_2.return_Chromosome(), phase_block_2.return_FirstVariantPosition(), phase_block_2.return_LastVariantPosition())

        # get list of overlapping variants
        overlapping_variants_list = list(set(phase_block_1.return_Variants().keys()).intersection(phase_block_2.return_Variants().keys()))
        n_variants_overlap = len(overlapping_variants_list)
        n_variants_flip_to_match = 0
        pct_switch = "NA"

        if n_variants_overlap > 0:
          for variant in overlapping_variants_list:
            if phase_block_1.return_Variants()[variant][0].return_Genotype() == phase_block_2.return_Variants()[variant][0].return_Genotype()[::-1]:
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

          extended_pb_dict1[pb1][pb2] = [pb1, phase_block_1.return_Chromosome(), phase_block_1.return_FirstVariantPosition(), phase_block_1.return_LastVariantPosition(), pb2, phase_block_2.return_Chromosome(), phase_block_2.return_FirstVariantPosition(), phase_block_2.return_LastVariantPosition(), min_overlap_position, max_overlap_position, length_overlap, n_variants_overlap, n_variants_flip_to_match, p_value_switch, p_value_no_switch, min_p_value, pct_switch, recommendation, graph_weight] # graph_weight must be final element of list to work with create_graph()
  
  return(extended_pb_dict1)

def super_phase_block_relationships(extended_pb_dict, phase_block_graph):

  graph_clusters = phase_block_graph.clusters()

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

  vertex_dict_pb_key = {}
  vertex_count = -1 # graph vertices start from 0

  for pb1 in sorted(extended_pb_dict.keys()):
    if pb1 + "_s1" not in vertex_dict_pb_key:
      vertex_count += 1
      vertex_dict_pb_key[pb1 + "_s1"] = vertex_count
    for pb2 in sorted(extended_pb_dict[pb1].keys()):
      if pb2 + "_s2" not in vertex_dict_pb_key:
        vertex_count += 1
        vertex_dict_pb_key[pb2 + "_s2"] = vertex_count

  vertex_dict_number_key = {v:k for k,v in vertex_dict_pb_key.items()}

  phase_block_relationship_dict = {}

  for cluster, index_list in super_set_dict_cluster_key.items():
    base_phase_block_index = None
    for index in index_list:
      if len(index_list) == 1:
        phase_block_relationship_dict[vertex_dict_number_key[index][:-3]] = ["NA"]*4
      else:
        if vertex_dict_number_key[index].endswith("_s1"):
          if base_phase_block_index is None:
            base_phase_block_index = index
            base_phase_block = vertex_dict_number_key[index][:-3]
          phase_block_relationship_dict[vertex_dict_number_key[index][:-3]] = [str(x) for x in [cluster, base_phase_block]]
          phase_block_relationship_dict[vertex_dict_number_key[index][:-3]].extend([str(x) for x in sum_graph_edges(phase_block_graph, index, base_phase_block_index)])

  return(phase_block_relationship_dict)

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

def sum_graph_edges(phase_block_graph, vertex_number_1, vertex_number_2):
  all_shortest_paths = phase_block_graph.get_all_shortest_paths(v = vertex_number_1, to = vertex_number_2)

  n_shortest_paths = len(all_shortest_paths)

  if n_shortest_paths == 0:
    return("No connection")
  if n_shortest_paths != 1:
    sys.exit("Number of shortest paths != 1")

  weight_sum = 0
  n_edges = 0
  for x,y in zip(all_shortest_paths[0][:-1], all_shortest_paths[0][1:]):
    weight_sum += phase_block_graph[x,y]
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
  pickle_path_pb1 = args.pb1
  with open(pickle_path_pb1, 'rb') as pickle_file_pb1:
    bam_phase_block_dictionary_pb1 = pickle.load(pickle_file_pb1)
    vcf_variants_dictionary_pb1 = pickle.load(pickle_file_pb1)

  # path to input pickle file of second sample
  pickle_path_pb2 = args.pb2
  with open(pickle_path_pb2, 'rb') as pickle_file_pb2:
    bam_phase_block_dictionary_pb2 = pickle.load(pickle_file_pb2)
    vcf_variants_dictionary_pb2 = pickle.load(pickle_file_pb2)

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

  extended_phase_block_dictionary = extend_phase_blocks(bam_phase_block_dictionary_pb1, bam_phase_block_dictionary_pb2, chrom, start, end)

  extended_phase_block_graph = create_graph(extended_phase_block_dictionary)

  # extended_pairwise_relationships = pairwise_phase_block_relationships(extended_phase_block_dictionary, extended_phase_block_graph)

  super_sets = super_phase_block_relationships(extended_phase_block_dictionary, extended_phase_block_graph)

  # write output for extended_phase_block_dictionary
  header_line = ["pb1", "pb1_Chromosome", "pb1_FirstVariantPosition", "pb1_LastVariantPosition", "pb2", "pb2_Chromosome", "pb2_FirstVariantPosition", "pb2_LastVariantPosition", "min_overlap_position", "max_overlap_position", "length_overlap", "n_variants_overlap", "n_variants_flip_to_match", "p_value_switch", "p_value_no_switch", "min_p_value", "pct_switch", "recommendation", "graph_weight"]
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extend_stats.tsv")
  output_file = open(output_file_path, 'w')
  output_file.write('\t'.join(header_line) + '\n')
  for k1,v1 in extended_phase_block_dictionary.items():
    for k2,v2 in extended_phase_block_dictionary[k1].items():
      if len(header_line) != len(v2):
        sys.exit("Length of header line and number of column not equal.")
      output_file.write('\t'.join([str(x) for x in v2]) + '\n')
  output_file.close()

  # write output and close files
  os.makedirs(args.output_directory, exist_ok = True)
  output_file_path = os.path.join(args.output_directory, args.output_prefix + ".extended_phase_blocks.tsv")
  output_file = open(output_file_path, 'w')

  header_row_list = summary_file.readline().strip().split()
  header_row_list.extend(["cluster", "cluster_base_phase_block", "switch_status", "n_reference_pb"])
  output_file.write('\t'.join(header_row_list) + '\n')
  for line in summary_file:
    pb_id = line.strip().split()[0]
    if pb_id in super_sets.keys():
      output_file.write('\t'.join(line.strip().split() + super_sets[pb_id]) + '\n')
    else:
      output_file.write('\t'.join(line.strip().split() + ["NA"]*4) + '\n')
  output_file.close()

  #output_file.close()

if __name__ == '__main__':
  main(args)
