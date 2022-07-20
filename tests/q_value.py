import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('LQSRPAAPPAPGPGQLTLR.ms2.carina.pin', sep='\t')

sns.displot(data=df, x='hyperscore', hue='label', bins=50, kde=True)
plt.show(block=True)

