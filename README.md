# Forward Gene Ontology Fingerprint Construction

This repository mainly about codes and data of forward Gene Ontology Fingerprint  (GOF) construction of the study: **A Multi-Resource Knowledge Empowered Gene Ontology Fingerprint Reveals Distinct Biomarker Clusters for Hepatocellular Carcinoma and Cholangiocarcinoma**.

***

In order to construct GOF, we collected the Pubmed abstracts of GI cancers, and filter out the abstract without any Gene mentioned,  then extracted the Gene and GO terms.

the table is descriptive statistics of abstracts and entities for GOF construction. for example, there are **27516** abstracts of LIHC with at least one Gene, and we totally extracted **5887** unique Gene and **3495** unique GO terms.

| **Cancer**                  | **Abstracts** | **Gene entities** | **GO term entities** |
| :-------------------------- | ------------: | ----------------: | -------------------: |
| Liver cancer (LIHC)         |        27,516 |             5,887 |                3,495 |
| Stomach Cancer (STAD)       |        25,514 |             5,611 |                2,659 |
| Biliary Tract Cancer (CHOL) |         4,641 |             1,837 |                1,035 |
| Esophagus Cancer (ESCA)     |         6,001 |             2,536 |                1,278 |
| Pancreatic Cancer (PAAD)    |        18,134 |             4,542 |                2,639 |
| Colorectal Cancer(CRC)      |        63,789 |             8,233 |                4,598 |
| **Total**                   |   **13,8035** |        **10,589** |            **5,984** |

***

There are three steps to construct GOF, includes:

* step1. Data retrieval and corpus construction
  * get abstract-gene mapping
  * get abstract-ontology mapping
* step2. Generation of the GOF
  * calculate P-value
  * adjust P-value
* step3. Calculate gene-gene similarity score (GGSS) 

***

### step1. Data retrieval and corpus construction

In the first step, data retrieval and corpus construction, we used the Entrez Programming Utilities (E-Utilities) tool to search and retrieve articles in MEDLINE format from the PubMed database. GO terms with functional annotations (**go.obo**, version released on November 26, 2018), including biological process, molecular function and cellular component terms, were downloaded from the GO Consortium. After data collection, a highly configurable NLP tool, MetaMap (2018 Release) , was applied to identify gene and GO term entities in the corpus, and their corresponding concept unique identifiers (CUIs) defined in the Unified Medical Language System (UMLS) Metathesaurus  were obtained.

#### get abstract-gene mapping

The input file of getting abstract-gene mapping was downloaded from the GIDB database. The format of input file is .xlsx, which contains the PubMed ID,  genes identified in the abstract by MetaMap. Totally six types of gastrointestinal (GI) cancer related files were downloaded, including oesophageal (ESCA), colorectal (COADREAD), liver (LIHC), pancreatic (PAAD), stomach (STAD) and bile duct (CHOL).

the output are in the **./step1_mapping/chol/mapping_gene_abstract.txt** with following format:

```html
gene1,abstract1;abstract2;......
gene2,abstract3;abstract4;......
... ...
```

the related code is **step1.1_get_GeneMapping.py**

#### get abstract-ontology mapping

First, we extracted the Gene Ontology terms corresponding concept unique identifiers (CUIs) defined in the Unified Medical Language System (UMLS) Metathesaurus (**./knol/UMLS_GO_CUI.txt**). Second, we processed the abstract by MetaMap, a highly configurable NLP tool. Third, we automatically identified the GO term entities in the abstract.

the results are in the **./step_mapping/chol/mapping_ontology_abstract.txt** with following format:

```html
GO1,abstract1;abstract2;......
GO2,abstract3;abstract4;......
... ...
```

the related code is **step1.2_get_OntologyMapping.py**

***

### step2. Generation of the GOF

#### calculate P-value

Given a gene and a GO term recognized in the corpus, we first calculated the statistical significance of the connectivity of each gene-GO pair using the unadjusted p-value generated from the hypergeometric test.

The **mathematical equation** is defined in paper, and the **code** is **calculatePvalue()** defined in **step2_OntologyFingerprint.py**

#### adjust P-value

Second, to convert the p-value into the adjusted p-value, we took into account several noised non-significant gene-GO relationships with the following three situations: i) abstracts that mentioned many genes; ii) abstracts that included general GO terms (e.g., “cell”, GO:0005623); and iii) hotspot genes or well-studied genes (e.g., TP53, and ERBB2) linked to massive biomedical studies. To eliminate these noised non-significant gene-GO relationships, we adopted the adjusted p-value defined by Tsoi LC.

The **mathematical equation** is defined in paper, and the **code** is **calculateAdjustedPvalue()** defined in **step2_OntologyFingerprint.py**

***

### step3. Calculate gene-gene similarity score (GGSS) 

The gene-gene similarity score (GGSS) was computed based on a modified inner product algorithm. The GGSS reflects the similarity between two genes, indicating their similar ontology fingerprints in a certain disease.

The **mathematical equation** is defined in paper, and the **code** is  **step3_GeneSimilarityScore.py**.

Furthermore, we determined the GGSS percentile threshold (**top 1%**) based on prior knowledge of the gene set enrichment analysis (GSEA) annotation file downloaded from https://www.gsea-msigdb.org/gsea/index.jsp (**./knol/c5.all.v6.2.symbols.gmt**).

To identify biologically justified subnetworks (highly dense interconnected clusters) in the GOF network,
we used the MCODE algorithm for network clustering and Cytoscape for visualization.



Descriptive statistics of genes in the GOF. For example, in the GOF of liver cancer, there are 5887 genes, with calculate the gene-gene similarity score (GGSS), we only found **5659** genes have shared the GO terms, means the GGSS more than 0.0, and after rank the gene-gene pairs with similarity score, the Top 1% GGSS has **4838** genes, after construct network and clustering, only **1768** genes are highly dense interconnected.

| **Cancer**           |   **GGSS** | **Top 1% of GGSS** | **Network clusters** |
| -------------------- | ---------: | -----------------: | -------------------: |
| Liver cancer         |      5,659 |              4,838 |                1,768 |
| Stomach Cancer       |      5,206 |              4,061 |                1,428 |
| Biliary Tract Cancer |      1,637 |              1,071 |                  430 |
| Esophagus Cancer     |      2,303 |              1,596 |                  697 |
| Pancreatic Cancer    |      4,332 |              3,486 |                1,377 |
| Colorectal Cancer    |      7,874 |              6,399 |                1,845 |
| **Total**            | **10,220** |          **8,774** |            **4,631** |



***

***

## **Evaluation of the discrimination ability of the identified gene biomarkers with machine learning**

To assess the identified genes that might be used for GI cancer classification, RNA-Seq data from the TCGA and ICGC databases were classified into two groups: the hepatobiliary tumours versus adjacent tissues group and the hepatobiliary tumours versus other GI cancer tissues group.

We generated six different machine learning models with the library ‘sklearn’ in Python3: 

* decision tree (DT), 
* k nearest neighbour (kNN), 
* logistic regression (LR), 
* naive Bayesian (NB), 
* random forest (RF) 
* and support vector machine (SVM). 

The performance of the model was evaluated using the expression of the identified genes as classifiers and the corresponding AUC of each model to determine the discrimination power of the two-classification tasks. 

We preprocessed the data with MinMaxScale, which scaled the gene expression value to a range between zero and one. Measurements such as accuracy, sensitivity, specificity and F1 score were also used to estimate the performance of the model. To avoid the batch effect of cross-platform gene expression data, we performed a 5-fold cross-validation on the TCGA and ICGC datasets, respectively.

The **code** and detailed information are in **./machine_learning_classifier/main.ipynb**.



## Contacts

zonghui@tongji.edu.cn, Tongji University, Shanghai, 200092, China

nadger_wang@139.com, Shanghai Eastern Hepatobiliary Surgery Hospital, Shanghai, 200438,
China



