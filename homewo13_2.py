import pandas as pd
from rpy2 import robjects
from rpy2.robjects import Formula
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
base = importr("base")
stats = importr("stats")
DESeq2 = importr("DESeq2")
# Прочитали файл 
a = pd.read_csv("colon_cancer_tumor_vs_normal_paired_counts.tsv", sep="\t", index_col=0)
print(a)

meta = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5}, index=a.columns)
meta["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta["Tissue"]), ref="Normal")
# Рассчитали коэффициенты нормализации
dds = DESeq2.DESeqDataSetFromMatrix(countData=a, colData=meta, design=Formula("~ Tissue"))
dds = DESeq2.DESeq(dds)
res = DESeq2.results(dds, name="Tissue_Tumor_vs_Normal")
res = DESeq2.lfcShrink(dds, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res = pd.DataFrame(base.as_data_frame(res))
res.index = counts.index
res = res.sort_values("padj")
res = res.loc[res["padj"] < 0.05]
res = res.loc[res["log2FoldChange"].abs() >= 1]
res.to_csv("Deseq_results_un.tsv", sep="\t")
un = res.loc[res["padj"] != 0, "padj"]

meta1 = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5}, index=a.columns)
meta1["Patient"] = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
meta1["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta["Tissue"]), ref="Normal")
# Рассчитали коэффициенты нормализации
dds1 = DESeq2.DESeqDataSetFromMatrix(countData=a, colData=meta1, design=Formula("~ Tissue + Patient"))
dds1 = DESeq2.DESeq(dds1)
res1 = DESeq2.results(dds1, name="Tissue_Tumor_vs_Normal")
res1 = DESeq2.lfcShrink(dds1, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res1 = pd.DataFrame(base.as_data_frame(res1))
res1.index = counts.index
res1 = res1.sort_values("padj")
res1 = res1.loc[res1["padj"] < 0.05]
res1 = res1.loc[res1["log2FoldChange"].abs() >= 1]
res1.to_csv("Deseq_results_p.tsv", sep="\t")
p = res1.loc[res1["padj"] != 0, "padj"]
print(len(un)) # Число дифференциально экспрессированных генов при непарном DESeq2
print(len(p)) # при парном DESeq2
print(un.head(10))# Топ-10 генов при непарном DESeq2
print(p.head(10))# при парном
# Вывод: 8 общих генов
Результат:
    
3698
3757
RP11-474D1.3    2.262371e-50
CDH3            7.650295e-40
MMP11           1.018072e-34
WNT2            7.607009e-33
ATG9B           9.604360e-33
CEMIP           1.686399e-30
SPTBN2          3.559111e-30
GYLTL1B         7.032899e-30
C2CD4A          9.433401e-29
TRIB3           2.085784e-27
Name: padj, dtype: float64
RP11-474D1.3    9.340975e-45
CDH3            8.460871e-43
GYLTL1B         4.097605e-36
MMP11           7.595253e-32
COMP            3.421126e-31
WNT2            5.619486e-31
ATG9B           8.273185e-30
CEMIP           2.323999e-29
C2CD4A          4.684854e-27
FOXQ1           2.513921e-26
Name: padj, dtype: float64
