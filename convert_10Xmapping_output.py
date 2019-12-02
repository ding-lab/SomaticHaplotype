import sys

input_filepath = sys.argv[1]
output_filepath = sys.argv[2]

input_file = open(input_filepath, 'r')

variant_dictionary = {}

for line in input_file:
  if line.startswith("chrChromosome") or line.startswith("chr\t"):
    continue
  if line.startswith("chr"):
    chrom, pos, ref, alt, depth, ref_or_alt = line.strip().split()
    current_variant_key = ":".join([str(x) for x in [chrom, pos, ref, alt]])
    if ref_or_alt == "Ref-Support":
      current_ref_or_alt = "ref"
    else:
      current_ref_or_alt = "alt"
    if current_variant_key not in variant_dictionary:
      variant_dictionary[current_variant_key] = {'ref1':[], 'ref2':[], 'refNone':[], 'alt1':[], 'alt2':[], 'altNone':[]}
  else:
    line_split = line.strip().split()
    haplotype = "None"
    for x in line_split:
      if x.startswith("HP:i:"):
        haplotype = x.split(":")[2]
    for x in line_split:
      if x.startswith("BX:Z:"):
        this_barcode = x.split(":")[2]
        if this_barcode not in variant_dictionary[current_variant_key][current_ref_or_alt+haplotype]:
          variant_dictionary[current_variant_key][current_ref_or_alt+haplotype].append(this_barcode)
input_file.close()

output_file = open(output_filepath, 'w')
output_file.write('\t'.join(['variant_key','chromosome','position','ref', 'alt', 'ref_barcodes_H1', 'ref_barcodes_H2', 'ref_barcodes_None', 'alt_barcodes_H1', 'alt_barcodes_H2', 'alt_barcodes_None']) + '\n')
for k,v in variant_dictionary.items():
  chrom, pos, ref, alt = k.split(":")
  if len(v['ref1']) == 0:
    ref_barcodes_H1 = "NA"
  else:
    ref_barcodes_H1 = ';'.join(v['ref1'])
  if len(v['ref2']) == 0:
    ref_barcodes_H2 = "NA"
  else:
    ref_barcodes_H2 = ';'.join(v['ref2'])
  if len(v['refNone']) == 0:
    ref_barcodes_None = "NA"
  else:
    ref_barcodes_None = ';'.join(v['refNone'])
  if len(v['alt1']) == 0:
    alt_barcodes_H1 = "NA"
  else:
    alt_barcodes_H1 = ';'.join(v['alt1'])
  if len(v['alt2']) == 0:
    alt_barcodes_H2 = "NA"
  else:
    alt_barcodes_H2 = ';'.join(v['alt2'])
  if len(v['altNone']) == 0:
    alt_barcodes_None = "NA"
  else:
    alt_barcodes_None = ';'.join(v['altNone'])
  output_file.write('\t'.join([k, chrom, pos, ref, alt, ref_barcodes_H1, ref_barcodes_H2, ref_barcodes_None, alt_barcodes_H1, alt_barcodes_H2, alt_barcodes_None]) + '\n')
output_file.close()
