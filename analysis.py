import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv('results.pin', sep='\t')
print(df)


clean = df.groupby('scannr')['q_value'].min()
clean = clean[clean < 0.05]
print(f"Initial PSMs with q < 0.05: {clean.shape[0]}")

sns.displot(data=df, x='deltascore', hue='label', bins=100)
plt.show()