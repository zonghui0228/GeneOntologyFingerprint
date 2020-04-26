# *_* coding: utf-8 *_*
# @Author  : zong hui

# ========================================================
# get mapping of genes and pubmed abstracts

# ========================================================
# file format:
#
# gene1,abstract1;abstract2;......
# gene2,abstract3;abstract4;......
# ......

# ========================================================
import os
import xlrd

def gene_mapping(infile, outfile):
    # load GIDB pubmed-gene file
    data = xlrd.open_workbook(infile, "utf-8")
    sheet = data.sheets()[0]
    nrows = sheet.nrows
    ncols = sheet.ncols

    genes = []
    abstracts = []
    gene_abstract = {}
    for r in range(1, nrows):
        gene_name = sheet.cell(r, 0).value
        abstracts_count = sheet.cell(r, 1).value

        genes.append(gene_name)
        gene_abstract[gene_name] = []

        for c in range(2, ncols):
            if sheet.cell(r, c).value:
                gene_abstract[gene_name].append(sheet.cell(r, c).value.replace(".txt", ""))
                abstracts.append(sheet.cell(r, c).value.replace(".txt", ""))
        if len(gene_abstract[gene_name]) != int(abstracts_count):
            print("[INFO]: number not match!, {}".format(gene_name))
            break

    print("genes:{}, pubmeds:{}".format(len(set(genes)), len(set(abstracts))))
    with open(outfile, 'w') as outf:
        for k, v in gene_abstract.items():
            outf.write("{},{}\n".format(k, ';'.join(v)))
    return gene_abstract


if __name__ == "__main__":
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\肝癌修回更新.xlsx", "./step1_mapping/lihc/mapping_gene_abstract.txt")
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\胃癌修回更新.xlsx", "./step1_mapping/stad/mapping_gene_abstract.txt")
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\胆管癌修回更新.xlsx", "./step1_mapping/chol/mapping_gene_abstract.txt")
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\食道癌修回更新.xlsx", "./step1_mapping/esca/mapping_gene_abstract.txt")
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\胰腺癌修回更新.xlsx", "./step1_mapping/paad/mapping_gene_abstract.txt")
    # gene_mapping("E:\\data\\OntologyFingerprintData\\GIDB_gene_pubmed\\肠癌修回更新.xlsx", "./step1_mapping/crc/mapping_gene_abstract.txt")
    print("Done!")

