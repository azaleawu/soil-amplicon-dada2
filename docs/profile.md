# Shuai Wu · Soil microbial ecology

I study how soil microbial communities change with depth and climate, with a particular interest in drylands and permafrost regions. A soil profile lets me ask whether two communities sampled at different depths in one soil profile differ more than communities from distant sites. Answering that question across studies requires knowing which observations came from the same profile, which depths were sampled and which parts of a marker gene the sequencing reads actually cover.

I use R and Linux workflows to connect public sequence data to field sampling records and environmental measurements. When I harmonize a study, I want to be able to trace a community estimate back to its physical sample, sequencing library and source accession. I retain the original depth intervals and distinguish field measurements from environmental variables extracted from a public grid.

For bacterial and archaeal 16S, eukaryotic 18S and fungal ITS data, this repository documents my approach to primer handling, run-specific error learning and marker-specific taxonomic assignment. I compare exact sequence variants across libraries with compatible sequence coverage. Where assays differ, I keep their separate feature tables and choose a level of comparison justified by the ecological question.

[Public metagenomics workflows](https://github.com/azaleawu/oligarchy-permafrost-microbiomes) · [GitHub profile](https://github.com/azaleawu)

## 中文简介

我是武帅，研究干旱区和永久冻土区的土壤微生物生态。我关心群落如何随土层深度、气候及土壤条件变化。垂向剖面研究里，我尤其想知道，同一剖面中不同深度的群落差异，在什么条件下会超过不同地点之间的差异。

我的工作包括收集公开扩增子和宏基因组数据、核对原始采样信息，以及在 R 和 Linux 中建立可复现的分析流程。我保留物理样品、土壤剖面、测序文库和公共数据库登录号之间的对应关系，记录原始深度区间，并区分野外实测值与公共栅格提取值。

这套流程整理了我处理细菌和古菌 16S、真核微生物 18S 与真菌 ITS 数据时采用的分析思路。比较前先核对引物、方向和实际覆盖区间，按测序批次学习误差，并保留每一步的读数。这样可以追溯群落表中的每条记录，也便于判断哪些跨研究差异具有可比性。
