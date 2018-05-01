"""
read the openSNP genotype file
"""

from collections import defaultdict
import json

# the grep unearthed some weirdly formatted "snps", to ignore need to make sure
# that lines have the right name
snps_of_interest = []
for line in open('omniex_missense_acmg.txt'):
    snp_name = line.strip().split("\t")[0]
    if snp_name != "Name":
        snps_of_interest.append(snp_name)


# read the unprocessed file and just get absolute number of each genotype/rsid
# in addition generate sum of all genotypes for each rsid
snps = defaultdict(dict)
snps_genotype_counts = {}
with open("opensnp_genotypes.csv") as infile:
    for line in infile:
        line_list = line.strip().split("\t")
        if len(line_list) == 3 and line_list[1] in snps_of_interest:
            _, rsid, genotype = line_list
            if genotype in snps[rsid].keys():
                snps[rsid][genotype] += 1
            else:
                snps[rsid][genotype] = 1
            if rsid in snps_genotype_counts.keys():
                snps_genotype_counts[rsid] += 1
            else:
                snps_genotype_counts[rsid] = 1

# convert absolute numbers into relative ones
snps_relative_frequencies = defaultdict(dict)
for rsid in snps:
    for genotype in snps[rsid]:
        relative_frequency = snps[rsid][genotype] / snps_genotype_counts[rsid]
        snps_relative_frequencies[rsid][genotype] = relative_frequency

# put all into a single dict for export to JSON
output_dictionary = defaultdict(dict)

for rsid in snps:
    output_dictionary[rsid]["absolute_values"] = snps[rsid]
    output_dictionary[rsid]["relative_values"] = snps_relative_frequencies[rsid]
    output_dictionary[rsid]["number_of_observations"] = snps_genotype_counts[rsid]

print(json.dumps(output_dictionary,indent=2))
