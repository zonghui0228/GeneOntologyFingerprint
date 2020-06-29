# *_* coding: utf-8 *_*
# @Author  : zong hui

# 胆管癌，Biliary Tract Cancer, ====CHOL====
# 肠癌，Colorectal Cancer, ====CRC====
# 食道癌，Esophagus Cancer, ====ESCA====
# 肝癌，Biliary Tract Cancer, ====LIHC====
# 胰腺癌，Pancreatic Cancer, ====PAAD====
# 胃癌，Stomach Cancer, ====STAD====

import math
def gene_gene_similarity(infile, outfile):
    # calculate gene similarity based on modified inner product.
    print("Step1: Load Gene Ontology Fingerprint, from file:>>>>>>>>>>{}".format(infile))
    gene_ontology = {}  # 建立gene fingerprint.
    with open(infile, 'r') as inf:
        genes = inf.readline().strip().split(",")[1:]
        gene_ontology = {g:{} for g in genes}
        for line in inf:
            l = line.strip().split(",")
            o = l[0]
            values = l[1:]
            for i in range(len(values)):
                if float(values[i]) < 1.0:
                    gene_ontology[genes[i]][o] = float(values[i])
                else:
                    pass
    print("there are {} genes.".format(len(gene_ontology)))

    # 计算gene_gene pairs的相似性，方法是修饰过的內积.
    print("Step2: Calculate Genes Similarity Score by Modified Inner Product...")
    gene_gene_score = {}
    genes = list(gene_ontology.keys())
    score_0 = []
    for i in range(len(genes)):
        if i % 100 == 0:
            print("[INFO]: processing: {}/{}".format(i, len(genes)))
        i_go = gene_ontology[genes[i]].keys()
        for j in range(i + 1, len(genes)):
            j_go = gene_ontology[genes[j]].keys()
            gene_pair = genes[i] + '_' + genes[j]
            both_ij = set(i_go) & set(j_go)
            only_i = set(i_go) - set(j_go)
            only_j = set(j_go) - set(i_go)
            if len(both_ij) == 0:  # 如果两个gene的共有的oncology为0，则相似度也为0
                # gene_pairs_score = 0.0
                # gene_gene_score[gene_pair] = gene_pairs_score
                pass
            else:
                gene_pairs_score = 0.0
                gene_gene_score_numerator, gene_gene_score_demominator = 0.0, 0.0  # 初始化相似度计算公式的分子和分母
                for b in both_ij:
                    # print genes[i],b,gene_go[genes[i]][b], genes[j],b,gene_go[genes[j]][b]
                    if gene_ontology[genes[i]][b] != 0.0 and gene_ontology[genes[j]][b] != 0.0:
                        gene_gene_score_numerator += math.log(gene_ontology[genes[i]][b]) * math.log(gene_ontology[genes[j]][b]) # 相似度计算公式的分子
                    else:  # 每个gene与自身的fingerprint oncology的相关性都不会为0，此处为以防万一，若score_0里有东西，则需注意了。
                        if gene_ontology[genes[i]][b] == 0.0:
                            score_0.append(genes[i] + '-' + b)
                        if gene_ontology[genes[j]][b] == 0.0:
                            score_0.append(genes[j] + '-' + b)
                if gene_gene_score_numerator != 0:
                    gene_gene_score_demominator = max(1.0, 0.5 * (len(only_i) + len(only_j)))  # 相似度计算公式的分母
                    gene_pairs_score = gene_gene_score_numerator / gene_gene_score_demominator
                    gene_gene_score[gene_pair] = gene_pairs_score

    gene_gene_score = sorted(gene_gene_score.items(),key=lambda item: item[1],reverse=True)
    score_0 = list(set(score_0))
    print("following gene & ontology pvalue equals to 0:")
    print("".format(len(score_0)), score_0)
    print("there are {} gene pairs similarity scores larger than 0.".format(len(gene_gene_score)))

    # 输出genepairs的相似度结果
    print("Step3: Output Genes Similarity Results, file:>>>>>>>>>>{}".format(outfile))
    with open(outfile, 'w') as outf1:
        outf1.write("gene1,gene2,similarity_score\n")
        for (gene1_gene2, score) in gene_gene_score:
            gene1_2 = gene1_gene2.split('_')
            outf1.write("{},{},{}\n".format(gene1_2[0], gene1_2[1],
                                            str(score)))
    print("Finished!")

# gene_gene_similarity(infile="./step2_ontofing/chol/adjusted_pvalue.csv", outfile="./step3_similarityscore/chol/genepairs_similarity.txt")
# gene_gene_similarity(infile="./step2_ontofing/crc/adjusted_pvalue.csv", outfile="./step3_similarityscore/crc/genepairs_similarity.txt")
# gene_gene_similarity(infile="./step2_ontofing/esca/adjusted_pvalue.csv", outfile="./step3_similarityscore/esca/genepairs_similarity.txt")
# gene_gene_similarity(infile="./step2_ontofing/lihc/adjusted_pvalue.csv", outfile="./step3_similarityscore/lihc/genepairs_similarity.txt")
# gene_gene_similarity(infile="./step2_ontofing/paad/adjusted_pvalue.csv", outfile="./step3_similarityscore/paad/genepairs_similarity.txt")
# gene_gene_similarity(infile="./step2_ontofing/stad/adjusted_pvalue.csv", outfile="./step3_similarityscore/stad/genepairs_similarity.txt")

####################使用从GSEA下载的数据（c5.all.v6.2.symbols），使用同样的方法计算相似性。####################
def gsea_gene_similarity(infile):
    gene_go = {}
    with open(infile, 'r') as f:
        for line in f:
            l = line.strip().split("\t")
            for g in l[2:]:
                if g not in gene_go: gene_go[g] = [l[0]]
                else: gene_go[g].append(l[0])
    print(len(gene_go))
    return gene_go
# gene_go = gsea_gene_similarity("./knol/c5.all.v6.2.symbols.gmt")

def add_annotation_to_similarity_score(infile, outfile):
    gene_go = gsea_gene_similarity("./knol/c5.all.v6.2.symbols.gmt")
    with open(outfile, "w") as outf:
        with open(infile, "r") as inf:
            header = inf.readline().strip().split(",")
            outf.write("{},{},{},{}\n".format(header[0], header[1], header[2], "GSEA"))
            n = 0
            for line in inf:
                n += 1
                if n%1000==0: print(n)
                l = line.strip().split(",")
                if l[0] in gene_go and l[1] in gene_go:
                    shared_go = set(gene_go[l[0]]) & set(gene_go[l[1]])
                    outf.write("{},{},{},{}\n".format(l[0], l[1], l[2], len(shared_go)))
                else:
                    outf.write("{},{},{},{}\n".format(l[0], l[1], l[2], "NA"))

# add_annotation_to_similarity_score("./step3_similarityscore/chol/genepairs_similarity.txt", "./step3_similarityscore/chol/genepairs_similarity_gsea.txt")
# add_annotation_to_similarity_score("./step3_similarityscore/crc/genepairs_similarity.txt", "./step3_similarityscore/crc/genepairs_similarity_gsea.txt")
# add_annotation_to_similarity_score("./step3_similarityscore/esca/genepairs_similarity.txt", "./step3_similarityscore/esca/genepairs_similarity_gsea.txt")
# add_annotation_to_similarity_score("./step3_similarityscore/lihc/genepairs_similarity.txt", "./step3_similarityscore/lihc/genepairs_similarity_gsea.txt")
# add_annotation_to_similarity_score("./step3_similarityscore/paad/genepairs_similarity.txt", "./step3_similarityscore/paad/genepairs_similarity_gsea.txt")
# add_annotation_to_similarity_score("./step3_similarityscore/stad/genepairs_similarity.txt", "./step3_similarityscore/stad/genepairs_similarity_gsea.txt")


# 筛选前p(百分比)的gene-gene similarity score,并保留有GSEA注释的例子!.
def percentile(infile, outfile, p):
    with open(infile, 'r') as inf:
        lines = inf.readlines()
    cadidate_lines = lines[:int(len(lines)*p)] # 筛选前p%的数据
    print("gene-gene similarity scores in {}%：{}".format(p*100, len(cadidate_lines)))
    n = 0
    with open(outfile, 'w') as outf:
        for line in cadidate_lines:
            l = line.strip().split(",")
            if l[3] not in ["0", "NA"]:
                outf.write("{},{},{}\n".format(l[0],l[1],l[2]))
                n += 1
    print("gene-gene similarity scores in {}% and annotated in GSEA: {}".format(p*100, n))

# percentile("./step3_similarityscore/chol/genepairs_similarity_gsea.txt","./step3_similarityscore/chol/genepairs_similarity_gsea_001.txt",p=0.01)
# percentile("./step3_similarityscore/crc/genepairs_similarity_gsea.txt","./step3_similarityscore/crc/genepairs_similarity_gsea_001.txt",p=0.01)
# percentile("./step3_similarityscore/esca/genepairs_similarity_gsea.txt","./step3_similarityscore/esca/genepairs_similarity_gsea_001.txt",p=0.01)
# percentile("./step3_similarityscore/lihc/genepairs_similarity_gsea.txt","./step3_similarityscore/lihc/genepairs_similarity_gsea_001.txt",p=0.01)
# percentile("./step3_similarityscore/paad/genepairs_similarity_gsea.txt","./step3_similarityscore/paad/genepairs_similarity_gsea_001.txt",p=0.01)
# percentile("./step3_similarityscore/stad/genepairs_similarity_gsea.txt","./step3_similarityscore/stad/genepairs_similarity_gsea_001.txt",p=0.01)