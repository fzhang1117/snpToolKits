## xQTL_CandidatepGWAS_summary.py ##
## Summary Candidate Gene and pGWAS information as a sheet ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-22 ##

import sys, re
from itertools import islice

fl_CandidateGene = sys.argv[1]
fl_pGWASSummary = sys.argv[2]
fl_cpdID = sys.argv[3]
fl_output = sys.argv[4]

with open(fl_cpdID, 'r') as fh_cpdID:
    dic_cpdID = {}
    for line in fh_cpdID:
        line = line.strip('\n').split('\t')
        if dic_cpdID.get(line[0]) is None:
            dic_cpdID[line[0]] = line[1]

with open(fl_CandidateGene, 'r') as fh_CandidateGene:
    dic_CandidateGene = {}
    for line in islice(fh_CandidateGene, 1, None):
        line = line.strip('\n').split('\t')
        gene = line[8]
        ## output = [cpd, condition, xQTL, xQTL_chr, xQTL_start, xQTL_end, leadsnp, leadp, CandidateGene, GeneLocation, Distance, Pathway, PathwayAnnotation, v4_gene_model]
        output = line[0:11] + line[12: 14] + [line[15]]
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
    gene, cpdID, condition = line[0], line[8], line[9]
    cpd = dic_cpdID[cpdID]
    key = gene + '-' + cpdID + '-' + condition
    xQTL, leadp_cpd, location, distance, pathway, pathwayAnno, ZmV4, genesymbol, altName, zma_alias, zma_description, ara_ortholog, ara_alias, ara_annotation = line[10], line[15], line[17], line[18], line[19], line[20], line[21], line[22], line[27], line[28], line[29], line[30], line[31], line[32]
    key_anno = [cpd, xQTL, leadp_cpd, location,distance, pathway, pathwayAnno, ZmV4, genesymbol, altName, zma_alias, zma_description, ara_ortholog, ara_alias, ara_annotation]
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
            dic_summary[key]['summary_leadp'][trait] = min(leadp_pGWAS_temp, leadp_pGWAS)

trait_list = sorted(list(set(trait_list)))

title = ['gene', 'cpd', 'condition', 'xQTL', 'p_cpd', 'gene:leadp_location', 'gene:leadp_distance', 'pathwayID', 'pathwayAnno', 'gene_v4', 'genesymbol', 'altName', 'zma_alias', 'zma_description', 'ara_ortholog', 'ara_alias', 'ara_annotation'] + ['snpnum:' + i for i in trait_list] + ['minp:' + i for i in trait_list ]
output = [title]
for key in dic_summary.keys():
    key_split = key.split('-')
    key_anno = dic_summary[key]['key_anno']
    dic_snp, dic_p = {}, {}
    for trait in trait_list:
        if dic_summary[key]['summary_snp'].get(trait) is not None:
            dic_snp[trait] = dic_summary[key]['summary_snp'][trait]
            dic_p[trait] = dic_summary[key]['summary_leadp'][trait]
        else:
            dic_snp[trait] = 0
            dic_p[trait] = '-'
    entry_snp = [dic_snp[i] for i in sorted(dic_snp.keys())]
    entry_p = [dic_p[i] for i in sorted(dic_p.keys())]
    entry = key_split + key_anno + entry_snp + entry_p
    entry = [str(i) for i in entry]
    output.append(entry)

with open(fl_output, 'w') as fh_output:
    output = ['\t'.join(i) for i in output]
    output = '\n'.join(output)
    fh_output.write(output)
