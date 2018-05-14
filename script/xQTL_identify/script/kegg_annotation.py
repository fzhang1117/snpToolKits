## kegg_annotation.py ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-11 ##

"""
Input the query list and output kegg pathway
"""
import sys, os, re
from itertools import islice

fl_gene = sys.argv[1]
fl_out = sys.argv[2]
def my_dicpathwaybuild(path_pathway, path_list):
    with open(path_list, 'r') as fh_pathwaylist:
        list_pathway = fh_pathwaylist.read().strip('\n').split('\n')
        dic_pathway = {}
        for line in list_pathway:
            line = line.split('\t')
            if dic_pathway.get(line[0]) is None:
                dic_pathway[line[0]] = line[1]

    fls_gene = [path_pathway + i for i in os.listdir(path_pathway)]
    dic_gene = {}
    #pathway = [re.split('.txt|=', i)[1] for i in os.listdir(path_gene)]
    for fl in fls_gene:
        pathway = re.split('/genes_|.txt', fl)[1]
        with open(fl, 'r') as fh:
            for line in islice(fh, 1, None):
                line = line.strip('\n').split('\t')
                ## [Entrez_id, Maize_v3, Maize_v4, KO, EC, Gene_description]
                if dic_gene.get(line[1]) is None:
                    dic_gene[line[1]] = [pathway]
                else:
                    dic_gene[line[1]].append(pathway)
    dic_gene.pop('-')
    return(dic_pathway, dic_gene)

res_dicpathwaybuild = my_dicpathwaybuild(path_pathway = "../data/kegg_gene_info_zma/", path_list = "../data/zma_kegg_list.txt")
dic_pathway, dic_gene = res_dicpathwaybuild[0], res_dicpathwaybuild[1]

with open(fl_out, 'w') as fh_out:
    output = []
    with open(fl_gene, 'r') as fh_gene:
        for gene in fh_gene:
            gene = gene.strip('\n')
            if dic_gene.get(gene) is not None:
                pathways = dic_gene[gene]
                ## will be more faster
                annotations = [dic_pathway[i] for i in pathways]
                pathways, annotations = '|'.join(pathways), '|'.join(annotations)
                output.append([gene, pathways, annotations])
            else:
                output.append([gene, '-', '-'])
        output = ['\t'.join(i) for i in output]
        fh_out.write('\n'.join(output))
