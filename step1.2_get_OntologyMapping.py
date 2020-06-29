# *_* coding: utf-8 *_*
# @Author  : zong hui

# ========================================================
# get mapping of ontology(GO terms) and pubmed abstracts(linked to at least one gene)

# ========================================================
# file format:
# ontology1,abstract1;abstract2;......
# ontology2,abstract3;abstract4;......
# ......


import os
import MetaMapResParser



# ==========================================================
# 获得UMLS中的GO信息
def get_UMLS_GOinfo(infile, outfile):
    go_cui = {}
    cuis, gos = set(), set()
    with open(infile, 'r', encoding='UTF-8') as f:
        n = 0
        for line in f:
            n += 1
            if n % 10000 == 0: print(n)
            if "|||GO:" in line:
                cui = line.split('|')[0]
                go = line.split("|GO|")[0].split("|||")[1]

                cuis.add(cui)
                gos.add(go)
                if go not in go_cui:
                    go_cui[go] = set([cui])
                else:
                    go_cui[go].add(cui)
            else:
                pass
    print(len(set(cuis))) # 70303
    print(len(set(gos))) # 47200
    with open(outfile, 'w') as f2:
        for go, cui in go_cui.items():
                f2.write("{}\t{}\n".format(go, ','.join(list(cui))))
# get_UMLS_GOinfo("D:\\UMLS\\2018AB\\2018AB\\META\\MRCONSO.RRF", "./knol/UMLS_GO_CUI.txt")


# ==========================================================
# filter out the abstract,, get the pubmed set linked at least one gene
def pubmed(gene_mapping):
    pubmed_set = set()
    with open(gene_mapping, 'r') as f:
        for line in f:
            pubmeds = line.strip().split(",")[1].split(";")
            for p in pubmeds:
                pubmed_set.add(p)
    print("number of pubmed: {}".format(len(pubmed_set)))
    return pubmed_set

# prepare ontology mapping data
def ontology_mapping(folder_path, outfile, pubmeds):
    gos = set()
    cuis = set()
    cui_go = {}
    # load UMLS GO id&cui
    with open("./knol/UMLS_GO_CUI.txt", 'r') as f1:
        for line in f1:
            l = line.strip("\n").split("\t")
            go = l[0]
            cui = l[1].split(",")

            gos.add(go)
            for c in cui:
                cuis.add(c)
                if c not in cui_go:
                    cui_go[c] = set([go])
                else:
                    cui_go[c].add(go)
    print("Count: allGO: {}, allCUI: {}".format(len(gos), len(cuis))) # 47200, 70303

    ontology_abstract = {}
    n = 0
    files = os.listdir(folder_path)
    for infile in files:
        n += 1
        if n % 1000 == 0: print("[INFO]: processing: {}/{}".format(n, len(files)))

        pubmed_id = infile[:-4]
        if pubmed_id in pubmeds:
            infile_path = os.path.join(folder_path, infile)
            terms = MetaMapResParser.single_parse(infile_path)
            abstract_cuis = [term[1] for term in terms]
            abstract_go_cuis = set(abstract_cuis) & cuis
            if abstract_go_cuis:
                for c in abstract_go_cuis:
                    abstract_gos = cui_go[c]
                    for g in abstract_gos:
                        if g not in ontology_abstract:
                            ontology_abstract[g] = [pubmed_id]
                        else:
                            ontology_abstract[g].append((pubmed_id))
            else:
                pass
    print("number of GO linked to these pubmeds: {}".format(len(ontology_abstract)))

    with open(outfile, 'w') as f3:
        for k, v in ontology_abstract.items():
            f3.write("{},{}\n".format(k, ';'.join(list(set(v)))))
    return ontology_abstract

if __name__ == "__main__":
    # ===========================================LIHC================================================
    # folder_path = "./medline_raw_metamap/LIHC_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/lihc/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/lihc/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    # ===========================================STAD================================================
    # folder_path = "./medline_raw_metamap/STAD_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/stad/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/stad/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    # ===========================================CHOL================================================
    # folder_path = "./medline_raw_metamap/CHOL_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/chol/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/chol/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    # ===========================================ESCA================================================
    # folder_path = "./medline_raw_metamap/ESCA_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/esca/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/esca/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    # ===========================================PAAD================================================
    # folder_path = "./medline_raw_metamap/PAAD_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/paad/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/paad/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    # ===========================================CRC================================================
    # folder_path = "./medline_raw_metamap/CRC_metamap_cancerfilter_all"
    # outfile = "./step1_mapping/crc/mapping_ontology_abstract.txt"
    # pubmeds = pubmed("./step1_mapping/crc/mapping_gene_abstract.txt")
    # ontology_mapping(folder_path, outfile, pubmeds)

    print("Done!")
