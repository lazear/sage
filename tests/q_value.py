import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']

PATH = 'b1906_293T_proteinID_01A_QE3_122212.ms2.carina.pin'
df = pd.read_csv(PATH, sep='\t')
passing = df[df.q_value <= 0.01]
df.loc[:,'Target/Decoy'] = df.label.apply(lambda x: 'Target' if x == 1 else 'Decoy')
hyperscore_cutoff = passing.hyperscore.values[-1]
print(f"Peptides originally at q<=0.01 (hyperscore): {passing.shape[0]}")


sns.displot(data=df, x='hyperscore', hue='Target/Decoy', kde=True)
plt.vlines(hyperscore_cutoff, *plt.ylim(), color='black', linestyle='dashed')
plt.text(hyperscore_cutoff + 2, y=plt.ylim()[1] - 500, s="1% FDR", rotation=90)
plt.suptitle("PXD001468")
plt.xlabel("Hyperscore")
plt.ylabel("PSM count")
plt.show(block=True)
fig, ax = plt.subplots()


