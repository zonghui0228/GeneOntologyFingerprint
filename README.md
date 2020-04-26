# A Multi-Resource Knowledge Empowered the Gene Ontology Fingerprint Approach for Disease Biomarker Prediction

This repository mainly about codes and data of Gene Ontology Fingerprint  (GOF) Construction.

## Forward gene ontology fingerprint (GOF) construction

### step1. Data retrieval and corpus construction

In the first step, data retrieval and corpus construction, we used the Entrez Programming Utilities (E-Utilities) tool to search and retrieve articles in MEDLINE format from the PubMed database. GO terms with functional annotations (go.obo, version released on November 26, 2018), including biological process, molecular function and cellular component terms, were downloaded from the GO Consortium. After data collection, a highly configurable NLP tool, MetaMap (2018 Release) , was applied to identify gene and GO term entities in the corpus, and their corresponding concept unique identifiers (CUIs) defined in the Unified Medical Language System (UMLS) Metathesaurus  were obtained.

#### get abstract-gene mapping

The input file of getting abstract-gene mapping was downloaded from the GIDB database. The format of input file is .xlsx, which contains the PubMed ID,  genes identified in the abstract by MetaMap. Totally six types of gastrointestinal (GI) cancer related files were downloaded, including oesophageal (ESCA), colorectal (COADREAD), liver (LIHC), pancreatic (PAAD), stomach (STAD) and bile duct (CHOL).

the output are in the **./step1_mapping/../mapping_gene_abstract.txt** with following format:

```
gene1,abstract1;abstract2;......
gene2,abstract3;abstract4;......
... ...
```

the related code is **step1.1_get_GeneMapping.py**

#### get abstract-ontology mapping

First, we extracted the Gene Ontology terms corresponding concept unique identifiers (CUIs) defined in the Unified Medical Language System (UMLS) Metathesaurus (**./knol/UMLS_GO_CUI.txt**). Second, we processed the abstract by MetaMap, a highly configurable NLP tool. Third, we automatically identified the GO term entities in the abstract.

the results are in the **./step_mapping/lihc/mapping_ontology_abstract.txt** with following format:

```
GO1,abstract1;abstract2;......
GO2,abstract3;abstract4;......
... ...
```

the related code is **step1.2_get_OntologyMapping.py**

### step2. Generation of the GOF

#### calculate P-value

Given a gene and a GO term recognized in the corpus, we first calculated the statistical significance of the connectivity of each gene-GO pair using the unadjusted p-value generated from the hypergeometric test.

The mathematical formula is defined in paper, and the code is **caculatePvalue()** defined in **step2_OntologyFingerprint.py**

#### adjust P-value

Second, to convert the p-value into the adjusted p-value, we took into account several noised non-significant gene-GO relationships with the following three situations: i) abstracts that mentioned many genes; ii) abstracts that included general GO terms (e.g., “cell”, GO:0005623); and iii) hotspot genes or well-studied genes (e.g., TP53, and ERBB2) linked to massive biomedical studies. To eliminate these noised non-significant gene-GO relationships, we adopted the adjusted p-value defined by Tsoi LC.

The mathematical formula is defined in paper, and the code is **caculateAdjustedPvalue()** defined in **step2_OntologyFingerprint.py**

### step3. Calculate gene-gene similarity score (GGSS) 

The gene-gene similarity score (GGSS) was computed based on a modified inner product algorithm. The GGSS reflects the similarity between two genes, indicating their similar ontology fingerprints in a certain disease.

The mathematical formula is defined in paper, and the code is  **step3_GeneSimilarityScore.py**.

Furthermore, we determined the GGSS percentile threshold based on prior knowledge of the gene set enrichment analysis (GSEA) annotation file downloaded from https://www.gsea-msigdb.org/gsea/index.jsp (**./knol/c5.all.v6.2.symbols.gmt**).

For each gene-gene pair, we defined a quaternion GGSS (A, B, SS, L), where A and B represent genes A and B, respectively, SS represents the similarity score between genes A and B, and L represents whether genes A and B were annotated in the same GO process in the GSEA. If so, L was labelled 1; otherwise, it was labelled 0. We ranked gene-gene pairs in descending order of the SS and removed those whose SSs were equal to zero.

## **Evaluation of the discrimination ability of the identified gene expression pattern**

To assess the identified genes that might be used for GI cancer classification, RNA-Seq data from the TCGA and ICGC databases were classified into two groups: the hepatobiliary tumours versus adjacent tissues group and the hepatobiliary tumours versus other GI cancer tissues group.

We generated six different machine learning models with the library ‘sklearn’ in Python: decision tree (DT), k nearest neighbour (kNN), logistic regression (LR), naive Bayesian (NB), random forest (RF) and support vector machine (SVM). The performance of the model was evaluated using the expression of the identified genes as classifiers and the corresponding AUC of each model to determine the discrimination power of the two-classification tasks. We preprocessed the data with MinMaxScale, which scaled the gene expression value to a range between zero and one. Measurements such as accuracy, sensitivity, specificity and F1 score were also used to estimate the performance of the model. To avoid the batch effect of cross-platform gene expression data, we performed a 5-fold cross-validation on the TCGA and ICGC datasets, respectively.

The code and detailed information are in **./machine_learning_classifier/main.ipynb**.



