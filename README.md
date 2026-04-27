# Graduation-Thesis

data and code for my thesis at Nankai University, about bioinformatics



我的毕业论文主题为“基于GO与STRING数据库的Spacedust保守基因簇预测工具的多物种可靠性评估——Reliability evaluation of Spacedust conserved gene clusters prediction tool based on GO and STRING databases”



本仓库展示了我的全部实验代码与原数据、输出结果。原定的10个物种分析，由于其中4个物种的数据注释情况与信息不符合要求，因此最终完整的output只有6个物种有效。


##code

仓库中展现了我的全部工作代码。

\------------ Introduction ------------

1 config：

##data

Species文件夹是整个实验用到的大部分原始数据，包括来自GO的STRING和NCBI的数据，以及Spacedust发掘保守基因簇所得到的保守矩阵数据。数据以物种为单位分类，每个物种内有5个子文件夹，以下提供介绍。共计6个物种，C.jejuni和S.meliloti被排除。


\------------ Introduction ------------

1 links：STRING的links\_detailed文件，以蛋白对protein pairs为单位，阐述蛋白对存在关系的概率大小。包含7种不同维度的特征，以及贝叶斯乘法算出的综合概率combined\_score。

2 aliases：STRING的aliases文件，以protein为单位，TaxonID + locus tag（如511145.b0001）为标准，记录蛋白质的别名以及对应的数据库来源。

3 goa：GO的annotation文件，以protein为单位，记录蛋白质在GO数据库的各类数据，如基因名称，被包含的GOTermID，evidence类型等。

4 feature：NCBI的feature\_table文件，以gene为单位，记录该基因组中基因的名称、位置、方向、genomic\_accession信息，及其对应编码的蛋白质的locus tag，refseqname等信息。该数据可以和氨基酸序列（.faa）整合，达到和General Feature Format（.gff）文件一样的效果。

5 conservation---querydb.lookup：Spacedust分析中产生的querydb的索引文件，包含每个基因对对应的genomic\_accession。

6 conservation---original\_conservation\_matrix：Spacedust发掘保守基因簇得到的最终的conservation\_matrix文件，每列为当前物种的一个同向基因对，每行为KEGG数据库中的一个reference物种。阐述当前物种的所有基因对的保守性情况，包括每个基因对在其他物种中的共现情况。

PS：conservation的文件来自于Spacedust的工作流程，是后续分析中得到的最重要的两个数据。该流程使用了来自NCBI官方的基因组氨基酸序列文件FASTA Amino Acid（.faa）以及feature\_table，以及Spacedust原论文中引用的reference数据库KEGG\_70。上述二者分别作为Spacedust的query和target。







Output文件夹是整个实验全部的分批次数据与图片输出，包括前中期的数据清洗与工程文件，以及后期和Spacedust的保守性分数的比对分析结果，共计跑通了6个物种（不包括C.jejuni和S.meliloti）。



\------------ Introduction ------------

1 output.log：代码的全部输出日志

2 match\_table.csv：GO的注释文件（GO-Annotation）加入了locus tag，用来和接下来的STRING的pairs数据匹配。拥有基因名、GOTermID、Evidence、Locustag

3 string\_links\_recalculated.csv：STRING的links文件加入了剔除channels之后重新计算的combined\_score。以protein pairs为单位。

4 go-string\_match\_pairs.csv：数据工程的最终文件，用来和Spacedust的保守性数据比对（以下统称PAIRS）。以proein pairs为单位，包含了该物种所有可能的pairs，及其被GO和STRING分别的注释情况。可能出现重复，因为一个pairs可以被多个GOTerm所注释。

5 query\_mapping.csv：Spacedust中的query.lookup加入了locus tag，用来和PAIRS比对。让原本没有ID，只有conservation score的数据拥有了ID，以此判断每个基因对的GO与STRING注释情况。

6 pairs overlap.png：描述GO和STRING注释的pairs占理论全组合的比例，与重叠情况。

7 high confidence pairs overlap.png：描述GO和STRING注释的pairs，以及其中的高置信度pairs占理论全组合的比例，与高置信度pairs的重叠情况。

8 spacedust pairs overlap.png：描述Spacedust的数据和GO/STRING注释的重叠情况

9 spacedust\_vs\_string\_combined\_scatter.png：描述Spacedust的conservation score和STRING的（剔除channels之后）combined score关系的散点图。

10 spacedust\_vs\_string\_original\_combined\_scatter.png：描述Spacedust的conservation score和STRING的combined score关系的散点图。

11 spacedust\_vs\_string\_XXX\_scatter.png：描述Spacedust的conservation score和STRING的某一个channel的score关系的散点图。

PS：以上的散点图数据来源于有STRING注释（score＞0）的pairs。

12 spacedust\_only\_distribution.png：描述仅在Spacedust的数据中存在，在PAIRS中不存在（GO和STRING均无注释）的pairs的conservation score分布。

13 intersection\_vs\_spacedust\_only\_boxplot.png：描述Spacedust和PAIRS的交集，与Spacedust中无注释的pairs的conservation score分布对比箱线图。

14 spacedust\_vs\_go\_strength.png：描述全集的conservation score分布，以不同GO evidence（High、Medium、Low、None）为标准来分箱。

15 spacedust\_tier\_go\_enrichment.png：描述全集中不同conservation scorequ区间内的GO evidence强度的构成比例柱状图。

16 high\_vs\_low\_conservation\_string\_boxplot.png：描述高保守性和低保守性pairs的STRING combined score分布对比箱线图，以150为conservation score的阈值。

17 spacedust\_tier\_string\_boxplot.png：描述全集中不同conservation score区间中的STRING combined score分布趋势箱线图。

18 go\_string\_confidence\_heatmap.png：描述全集中GO evidence X STRING combined score区间的样本分布与平均conservation score。

19 string\_go\_dual\_axis\_trend.png：双Y轴图，描述全集中STRING combined score和GO的注释比例随着conservation score改变的趋势。





\------------ Others ------------

faa是NCBI Assembly的官方基因组核苷酸序列文件，NCBI_ten_genomes是经过转换操作的核苷酸序列文件。其他信息不变，把Product\_Accesion换成Genomic\_Accesion，符合Spacedust预测工具的querydb的格式。

KEGG_70是Spacedust的reference数据库clusterdb，作为预测工具的target数据。通过KEGG_70和NCBI_ten_genomes的数据，运行runst.ipynb流程（在Linux预装好工具）即可获得输出。

以上提到的所有数据在我的Google Drive里展示：https://drive.google.com/drive/folders/1ArZGq40EQeVjJmZZ8ohVigx_3YmbPAtp?usp=drive_link
