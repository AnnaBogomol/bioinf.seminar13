import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
student = pd.read_csv("student_top10.tsv", sep="\t", index_col=0)
mw = pd.read_csv("MannWh_top10.tsv", sep="\t", index_col=0)
deseq = pd.read_csv("Deseq_top10.tsv", sep="\t", index_col=0)
deseq = deseq.melt(var_name="DESeq2", value_name="p-adj")
student = student.melt(var_name="Student", value_name="p-value") 
mw = mw.melt(var_name="MannWhitneyu", value_name="p.value")
sns.violinplot(x="DESeq2", y="p-adj", data=deseq) 
plt.tight_layout()
plt.savefig("DESeq2.jpg")
sns.violinplot(x="Student", y="p-value", data=student)
plt.tight_layout()
plt.savefig("student.jpg")
sns.violinplot(x="MannWhitneyu", y="p.value", data=mw)
plt.tight_layout()
plt.savefig("MannWhitneyu.jpg")
plt.close()
