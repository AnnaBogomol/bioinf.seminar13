import pandas as pd
from scipy.stats import ttest_rel, ttest_ind, rankdata, mannwhitneyu
df = pd.read_csv('colon_cancer_tumor_vs_normal_unpaired_FPKM.tsv', sep = '\t', index_col = 0)
df["unpaired"] = [ttest_ind(df.loc[gene].iloc[0:5], df.loc[gene].iloc[5:10])[1] for gene in df.index]
df_student = df["unpaired"]
df_student = df_student[df_student < 0.05].sort_values()
from rpy2 import robjects
from rpy2.robjects import Formula
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
base = importr("base")
stats = importr("stats")
DESeq2 = importr("DESeq2")
# Прочитали файл 
counts = pd.read_csv("colon_cancer_tumor_vs_normal_unpaired_counts.tsv", sep="\t", index_col=0)

meta = pd.DataFrame({"Tissue": ["Tumor"]*5 + ["Normal"]*5}, index=counts.columns)
meta["Tissue"] = stats.relevel(robjects.vectors.FactorVector(meta["Tissue"]), ref="Normal")
# Рассчитали коэффициенты нормализации
dds = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta, design=Formula("~ Tissue"))
dds = DESeq2.DESeq(dds)
res = DESeq2.results(dds, name="Tissue_Tumor_vs_Normal")
res = DESeq2.lfcShrink(dds, coef="Tissue_Tumor_vs_Normal", type="apeglm")
res = pd.DataFrame(base.as_data_frame(res))
res.index = counts.index
res = res.sort_values("padj")
res = res.loc[res["padj"] < 0.05]
res = res.loc[res["log2FoldChange"].abs() >= 1]
res.to_csv("DESeq2_results_unun.tsv", sep="\t")
df_deseq = res.loc[res["padj"] != 0, "padj"]
print(df_deseq.head(10))# Топ-10 генов при непарном DESeq2
df_deseq.head(10).to_csv("Deseq_top10.tsv", sep="\t")
df["mannwhitneyu"]=df.apply(lambda x: mannwhitneyu(x[:5], x[5:10])[1], axis = 1)
df_mw = df.loc[df["mannwhitneyu"] != 0, "mannwhitneyu"]
df_mw = df_mw.sort_values()
print(df_mw.head(10))# топ-10 генов при критерии Манна-Уитни
df_mw.head(10).to_csv("MannWh_top10.tsv", sep="\t")
print(df_student.head(10))# Топ-10 генов при непарном t-критерии Стьюдента
df_student.head(10).to_csv("student_top10.tsv", sep="\t")
# Общие гены при DESeq2 и критерии Манна-Уитни: 0
# Общие гены при критерии Манна-Уитни и t-критерии Стьюдента: 0 
# Общие гены при DESeq2 и t-критерии Стьюдента: IER5L, FUT1, C17orf96 

# Вывод: DESeq2 и t-критерий Стьюдента мощнее, так как дают меньшие значения p-value, 
# за счет того, что они делают предположение о виде распределения. 
# При использовании критерия Манна-Уитмана платим за то, что не делаем предположений о природе выборки 
# меньшей мощностью теста, из-за чего находим меньшее число дифференциальных генов

Результат:
FABP6       3.701666e-34
ETV4        1.082326e-30
IER5L       2.268427e-26
KRT80       4.786961e-26
FUT1        4.430398e-23
C17orf96    1.933375e-22
CLDN1       5.494094e-22
ATG9B       2.092626e-21
KIAA1257    2.687224e-20
SLC51B      2.687224e-20
Name: padj, dtype: float64
SFTA2            0.004850
RAET1L           0.005580
CTD-2147F2.1     0.005580
LINC00460        0.005580
AC007128.1       0.005963
RP5-884M6.1      0.005963
VAC14-AS1        0.005963
RP11-399O19.9    0.005963
CST1             0.005963
LINC00858        0.005963
Name: mannwhitneyu, dtype: float64
C17orf96    3.038037e-08
IER5L       8.053047e-08
FUT1        9.848163e-08
CDH3        5.472875e-07
FXYD5       8.408204e-07
ZNHIT2      1.611916e-06
CLCA4       2.305243e-06
ACADSB      2.840350e-06
MT1F        2.856975e-06
PIGN        3.159708e-06
Name: unpaired, dtype: float64
