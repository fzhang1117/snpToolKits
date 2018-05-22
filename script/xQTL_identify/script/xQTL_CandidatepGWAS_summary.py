## xQTL_CandidatepGWAS_summary.py ##
## Summary Candidate Gene and pGWAS information as a sheet ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-22 ##

import sys, re
from itertools import islice

fl_CandidateGene = sys.argv[1]
fl_pGWASSummary = sys.argv[2]

with open(fl_CandidateGene, 'r') as fh_CandidateGene:
    dic_CandidateGene = {}
    for line in islice(fh_CandidateGene, 1, None):
        line = line.strip('\n').split('\t')
        gene = line[8]
        ## output = [cpd, condition, xQTL, xQTL_chr, xQTL_start, xQTL_end, leadsnp, leadp, CandidateGene, GeneLocation, Distance, Pathway, PathwayAnnotation, v4_gene_model]
        output = line[0:10] + line[12: 14] + [line[15]]
        if dic_CandidateGene.get(gene) is None:
            dic_CandidateGene[gene] = [output]
        else:
            dic_CandidateGene[gene].append(output)

with open(fl_pGWASSummary, 'r') as fh_pGWASSummary:
    result_merge = []
    for line in islice(fh_pGWASSummary, 1, None):
        line = line.strip('\n').split('\t')
        gene = line[7]
        if dic_CandidateGene.get(gene) is not None:
            value = dic_CandidateGene[gene]
            for line2 in dic_CandidateGene[gene]:
                entry = [gene] + line[0: 7] + line2 + line[9: ]
                result_merge.append(entry)

dic_summary, trait_list = {}, []

for line in result_merge:
    gene, cpd, condition = line[0], line[8], line[9]
    key = gene + '-' + cpd + '-' + condition
    xQTL, leadp_cpd, location, pathway, pathwayAnno, ZmV4, genesymbol, altName, zma_alias, zma_description, ara_ortholog, ara_alias, ara_annotation = line[10], line[15], line[17], line[18], line[19], line[20], line[21], line[26], line[27], line[28], line[29], line[30], line[31]
    key_anno = [xQTL, leadp_cpd, location, pathway, pathwayAnno, ZmV4, genesymbol, altName, zma_alias, zma_description, ara_ortholog, ara_alias, ara_annotation]
    trait, snp_pGWAS, leadp_pGWAS = line[1], line[2], float(line[5])
    trait_list.append(trait)
    if dic_summary.get(key) is None:
        dic_summary[key] = {'key_anno': key_anno, 'summary_snp':{}, 'summary_leadp':{}}
        dic_summary[key]['summary_snp'][trait] = 1
        dic_summary[key]['summary_leadp'][trait] = leadp_pGWAS
    else:
        if dic_summary[key]['summary_snp'].get(trait) is None:
            dic_summary[key]['summary_snp'][trait] = 1
            dic_summary[key]['summary_leadp'][trait] = leadp_pGWAS
        else:
            dic_summary[key]['summary_snp'][trait] += 1
            leadp_pGWAS_temp = dic_summary[key]['summary_leadp'][trait]
            dic_summary[key]['summary_leadp'][trait] = min(leadp_pGWAS, leadp_pGWAS)

trait_list = sorted(list(set(trait_list)))
print dic_summary
print trait_list
