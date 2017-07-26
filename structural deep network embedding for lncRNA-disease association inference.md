# Structural Deep Network Embedding for LncRNA-disease Association Inference

---
###### Author : Li Liu
###### Start Date : 2017-7-23

## Abstract
### Motivation
### Results

---

## Introduction


---

## Materials and Method
### Pipeline

![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\steps.jpg) 

> 利用lncRNA_Disease_Matrix.txt文件，a network based method for analysis of LncRNA-Disease associations,形成两个网络

![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\绘图1.gif) 

The **lncDN** and **DlncN**.
- (a) *In lncDN, each node corresponds to a distinct disease, colored based on the disease class [ICD-10] to which it belongs.* The names of 21 disease classes are shown on the right panel. A link between two diseases exists if they share at least one implicated lncRNA. The size of the node is proportional to the degree of the node in lncRNA-disease association network. 
- (b) *In DlncN, each node is a lncRNA, with two lncRNAs being connected if they are implicated in the same disease.* The size of each node is proportional to the number of diseases in which the lncRNA is implicated. The color of a node is based on the class of diseases in which the corresponding lncRNA implicated. Nodes are light purple if the corresponding lncRNAs are associated with more than one disease class. We label the lncRNAs implicated in more than five diseases by their names. 

**lncDN**

![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\disease159-ps.jpg) 

**DlncN**
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\lncRNA1171.png) 


### Structural Deep Network Embedding

#### for every vertex of node to compute embedding repersentation
```math
y^{(K)}_i
```
#### disease and lncRNA into common r-dimension
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\SDNE.jpg) 

#### 每一个顶点得到同一纬度的embedding向量表示
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\捕获3.JPG) 

### Similarity Learning

#### 根据embedding向量表示计算相似性得分

The similarity of two vertices u and v is calculated as follows,
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\公式1.jpg)

where d is the dimension, and ui, vi are the components of vector u and v respectively. 

#### 分别得到disease和lncRNA根据embedding向量表示计算的相似性得分

### Association Discovering

1. lncRNA-based similarity inference(LBSI)

    **LBSI** predicts a disease–lncRNA association , if a disease di is associated with a lncRNA that has a similar lncRNA lj. For a pair of (di,lj), a confidence score of the pair is calculated as,
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\公式3.jpg)
    
    where sim(lt,lj) is the similarity between lt and lj, and ait=1 if there is an existing association between li and lt otherwise ait=0.
2. disease-based similarity inference(DBSI)

    **DBSI** predicts a disease–lncRNA association,if a disease di is similar with a  disease that has an existing association with a lncRNA lj. For a pair of (di,lj), a confidence score of the pair is calculated as,
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\公式2.jpg)

    where sim(di,dt) is the similarity between di and dt; and atj=1 if there is an existing association between dt and lj otherwise atj=0. 

Operationally, for a disease di or a lncRNA lj as the input query, the DBSI and LBSI confidences are normalized as, 

![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\公式5.jpg)
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\公式4.jpg)

where Max(di,) is the maximum confidence and Min(di,) is the minimum confidence for di, and Max(,lj) is the maximum confidence and Min(,lj) is the minimum confidence for lj.

---

## Experiments
### Data preparation
- [ ] 文件
- lncRNA_Disease_Matrix.txt 
- lncRNA117.txt
- disease159.txt

NAME | NUMBER | FROM| SPECIES
---|---|---|---
lncRNA | 117| [==???==](http://XXX.com/)| ==???==
disease | 159| [==???==](http://XXX.com/)| 21(ICD-10[International Classification of Diseases])

- [ ] 疾病的分类文件
- disease_classificate.xlsx
- [ ] lncRNA的分类文件
- ==?????==
### Baseline Algorithms/Similaries
- [ ] 其他相似性文件--两两组合

其他相似性文件| 
---| 
disease_GIPSim_matrix.txt  | 
diseasesim_DOSEwang_matrix.txt | 
diseasesim_funSim_matrix.txt  | 
diseasesim_icodhprd.txt  | 
diseasesim_jaccard_go_matrix.txt  | 
  | 
lncRNA_GIPSim_matrix.txt  | 
LncRNA_SeqSim_Matrix.txt | 

- [ ] ==network embedding method==

network embedding method| 
---| 
DeepWalk  | 
LINE | 
GraRep | 
。。。  | 
- [ ] ==other lncRNA-disease association prediction method==

other lncRNA-disease association prediction method| 
---| 
LDAP  | 
RLSLD | 
RWR | 
。。。  | 
### Validation Metrics
10-fold cross-valiation

### Evaluation Metrics
#### precision@k
#### micro-F1 macro-F1
#### ROC curve / AUC
#### Recovered Fraction 

### Experiment Results
#### Network Reconstruction
precision@k

![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\NetworkReconstruction.jpg) 

#### Comparsion Experiment

##### Comparsion with other similarity method

###### ROC curve / AUC 
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\DBSI_ROC.png) 
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\LBSI_ROC.png) 


###### Recovered Fraction 


##### Comparsion with other network embedding method
==ROC curve / AUC 
Recovered Fraction==

##### Comparsion with other lncRNA-disease associations prediction method
==ROC curve / AUC Recovered Fraction==
![image](C:\Users\lecea1995\Desktop\lncRNA-prodige\2016-limin-bioinformatics-疾病lncRNA\result\其他方法对比柱状图_PS.png) 
#### Case Study
breast cancer
11/20   55%

Rank | lncRNA| Evidence| Description
---|---| ---|--- 
1 | BCAR4| LncRNADisease| BCAR4 is expressed in 27% of primary breast tumors. Forced expression of BCAR4 in human ZR-75-1 and MCF7 breast cancer cells resulted in cell proliferation in the absence of estrogen and in the presence of various antiestrogens.BCAR4 may be a good target for treating antiestrogen-resistant breast cancer.
2 | LSINCT5| LncRNADisease| Ovarian and breast tumours have also been associated with the expression of the LSINCT5 lncRNA; this transcript acts to target several other transcripts, including the antisense RNA NEAT-1 and the PSPC1 gene, which codes for a splicing regulatory factor
3 | MIR31HG| LncRNADisease| miR-31 and its host gene lncRNA LOC554202 (MIR31HG) are regulated by promoter hypermethylation in triple-negative breast cancer.Both miR-31 and the host gene LOC554202 are down-regulated in the TNBC cell lines of basal subtype and over-expressed in the luminal counterparts.
4 | PINC | LncRNADisease| In a finding of relevance to breast cancer pathogenesis, the mammary gland lncRNA PINC, whose genomic structure is substantially different between primates and rodents has been shown to function in both cell survival and cell cycle progression.
5 | GAS5| LncRNADisease| GAS5, a non-protein-coding RNA, controls apoptosis and is downregulated in breast cancer.
6 | SRA1| No lncRNA| 
7 | ZNFX1-AS1| LncRNADisease|SNORD-host RNA Zfas1 is a regulator of mammary development and a potential marker for breast cancer.ZFAS1 is highly expressed in the mammary gland and is down-regulated in breast tumors compared to normal tissue. ZFAS1 is a putative tumor suppressor gene. 
8 | HOTAIR| LncRNADisease| A subsequent study revealed that HOTAIR is overexpressed in approximately one quarter of human breast cancers.
9 | DSCAM-AS1| LncRNADisease|M41 (DSCAM-AS1) mRNA is expressed at a statistically significantly higher level in human breast cancer specimens than in normal human breast and benign lesions. 
10 | MEG3| LncRNADisease|MEG3 expression is lost. 
11 | IGF2-AS | No| 
12 | A130040M12Rik| No| 
13 | EPB41L4A-AS1| No| 
14 | LINC00032| No| 
15 | RRP1B| No lncRNA| 
16 | SPRY4-IT1| LncRNADisease|The long noncoding RNA SPRY4-IT1 increases the proliferation of human breast cancer cells by upregulating ZNF703 expression. 
17 | Yiya| No| 
18 | WRAP53| No| 
19 | PVT1| LncRNADisease| Amplification of PVT1 contributes to the pathophysiology of ovarian and breast cancer.
20 | WT1-AS| NO| 

---

## Disussion