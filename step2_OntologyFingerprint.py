# *_* coding: utf-8 *_*
# @Author  : zong hui

# Construct Gene-Ontology Matrix
# Caculate Pvalue by hypergeometric test
# Caculate Adjusted-Pvalue

# import pandas as pd
import scipy.stats as stats
import math
import csv
import time

def build_matrix(genemapping, ontologymapping, abstract_number):
    gene = []
    gene_abstract = {}
    with open(genemapping, 'r') as f1:
        for line in f1:
            l = line.strip().split(",")
            gene.append(l[0])
            gene_abstract[l[0]] = l[1].split(";")

    ontology = []
    ontology_abstract = {}
    with open(ontologymapping, 'r') as f2:
        for line in f2:
            l = line.strip().split(",")
            ontology.append(l[0])
            ontology_abstract[l[0]] = l[1].split(";")
    print("number of pubmeds:{}".format(abstract_number))
    print("number of genes:{}".format(len(gene)))
    print("number of ontology:{}".format(len(ontology)))

    gene_ontology_matrix = {}
    for g in gene:
        gene_ontology_matrix[g] = {}
        for o in ontology:
            g_abstract = set(gene_abstract[g])
            o_abstract = set(ontology_abstract[o])

            gene_ontology_matrix[g][o] = {}
            a = len(g_abstract & o_abstract)
            b = len(g_abstract - o_abstract)
            c = len(o_abstract - g_abstract)
            d = abstract_number - a - b - c
            gene_ontology_matrix[g][o]["a"] = a
            gene_ontology_matrix[g][o]["b"] = b
            gene_ontology_matrix[g][o]["c"] = c
            gene_ontology_matrix[g][o]["d"] = d

    return gene, ontology, gene_ontology_matrix
# build_matrix("./step1_mapping/chol/mapping_gene_abstract.txt", "./step1_mapping/chol/mapping_ontology_abstract.txt", abstract_number = 4642)


def caculatePvalue(genemapping, ontologymapping, abstract_number, pvalue_file):

    gene, ontology, gene_ontology_matrix = build_matrix(genemapping, ontologymapping, abstract_number)

    if pvalue_file: ################pvalufile
        csvfile = open(pvalue_file, 'w', newline='')
        writer = csv.writer(csvfile)
        writer.writerow([""] + gene)

    start_time = time.process_time()
    for o in ontology:
        used_time = round(time.process_time() - start_time, 2)
        if ontology.index(o) % 100 == 0: print("processing: {}/{}, time using: {}s".format(ontology.index(o), len(ontology), used_time))

        o_pvalue = [o]
        for g in gene:

            a = gene_ontology_matrix[g][o]["a"]
            b = gene_ontology_matrix[g][o]["b"]
            c = gene_ontology_matrix[g][o]["c"]
            d = gene_ontology_matrix[g][o]["d"]

            if a == 0:
                pvalue = 1.0
                gene_ontology_matrix[g][o]["pvalue"] = pvalue
            else:
                oddsration, pvalue = stats.fisher_exact([[a, b], [c, d]])
                gene_ontology_matrix[g][o]["pvalue"] = pvalue

            o_pvalue.append(pvalue)

        if pvalue_file:  ################pvalufile
            writer.writerow(o_pvalue) ################pvalufile
    if pvalue_file:         ################pvalufile
        csvfile.close() ################pvalufile

    return gene, ontology, gene_ontology_matrix
# caculatePvalue("./step1_mapping/lihc/mapping_gene_abstract.txt", "./step1_mapping/lihc/mapping_ontology_abstract.txt", 27516, "./step2_ontofing/lihc/pvalue.csv")
# caculatePvalue("./step1_mapping/stad/mapping_gene_abstract.txt", "./step1_mapping/stad/mapping_ontology_abstract.txt", 25514, "./step2_ontofing/stad/pvalue.csv")
# caculatePvalue("./step1_mapping/chol/mapping_gene_abstract.txt", "./step1_mapping/chol/mapping_ontology_abstract.txt", 4641, "./step2_ontofing/chol/pvalue.csv")
# caculatePvalue("./step1_mapping/esca/mapping_gene_abstract.txt", "./step1_mapping/esca/mapping_ontology_abstract.txt", 6001, "./step2_ontofing/esca/pvalue.csv")
# caculatePvalue("./step1_mapping/paad/mapping_gene_abstract.txt", "./step1_mapping/paad/mapping_ontology_abstract.txt", 18134, "./step2_ontofing/paad/pvalue.csv")
# caculatePvalue("./step1_mapping/crc/mapping_gene_abstract.txt", "./step1_mapping/crc/mapping_ontology_abstract.txt", 63789, "./step2_ontofing/crc/pvalue.csv")

def caculateAdjustedPvalue(genemapping, ontologymapping, abstract_number, adjusted_file):

    gene, ontology, gene_ontology_matrix = caculatePvalue(genemapping, ontologymapping, abstract_number, pvalue_file=False)

    csvfile = open(adjusted_file, 'w', newline='')
    writer = csv.writer(csvfile)
    writer.writerow([""] + gene)

    start_time = time.process_time()
    for o in ontology:

        used_time = round(time.process_time() - start_time, 2)
        if ontology.index(o) % 100 == 0: print("processing: {}/{}, time using: {}s".format(ontology.index(o), len(ontology), used_time))

        o_adjustedpvalue = [o]
        for g in gene:
            pvalue = gene_ontology_matrix[g][o]["pvalue"]
            if pvalue == 1.0:
                adjusted_pvalue = 1.0

            if pvalue < 1.0:
                numerator = 0
                denominator = 0
                for k in ontology:
                    if gene_ontology_matrix[g][k]["pvalue"] < 1.0:
                        numerator += 1
                    if gene_ontology_matrix[g][k]["pvalue"] < 1.0 and gene_ontology_matrix[g][k]["pvalue"] >= pvalue:
                        denominator += 1
                adjusted_pvalue = min(1.0, pvalue * float(numerator) / float(denominator))

            gene_ontology_matrix[g][o]["adjusted_pvalue"] = adjusted_pvalue

            o_adjustedpvalue.append(adjusted_pvalue)
        writer.writerow(o_adjustedpvalue)
    csvfile.close()
    return gene, ontology, gene_ontology_matrix

# caculateAdjustedPvalue("./step1_mapping/lihc/mapping_gene_abstract.txt", "./step1_mapping/lihc/mapping_ontology_abstract.txt", 27516, "./step2_ontofing/lihc/adjusted_pvalue.csv")
# caculateAdjustedPvalue("./step1_mapping/stad/mapping_gene_abstract.txt", "./step1_mapping/stad/mapping_ontology_abstract.txt", 25514, "./step2_ontofing/stad/adjusted_pvalue.csv")
# caculateAdjustedPvalue("./step1_mapping/chol/mapping_gene_abstract.txt", "./step1_mapping/chol/mapping_ontology_abstract.txt", 4641, "./step2_ontofing/chol/adjusted_pvalue.csv")
# caculateAdjustedPvalue("./step1_mapping/esca/mapping_gene_abstract.txt", "./step1_mapping/esca/mapping_ontology_abstract.txt", 6001, "./step2_ontofing/esca/adjusted_pvalue.csv")
# caculateAdjustedPvalue("./step1_mapping/paad/mapping_gene_abstract.txt", "./step1_mapping/paad/mapping_ontology_abstract.txt", 18134, "./step2_ontofing/paad/adjusted_pvalue.csv")
# caculateAdjustedPvalue("./step1_mapping/crc/mapping_gene_abstract.txt", "./step1_mapping/crc/mapping_ontology_abstract.txt", 63789, "./step2_ontofing/crc/adjusted_pvalue.csv")



