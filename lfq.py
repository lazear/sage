import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm, gaussian_kde
from scipy.spatial.distance import cosine, mahalanobis, euclidean
from scipy.signal import find_peaks, savgol_filter

df = pd.read_csv("b1906_293T_proteinID_01A_QE3_122212.sage.pin", sep="\t")
df = df[(df.label == 1) & (df.posterior_error <= -2)]
co = df.groupby("peptide")["scannr"].count()
# df = df[df.peptide.isin(co[co > 10].index)]
# df = df.sort_values(by='rt')
quant = pd.read_csv("quant.csv")
# df = df[df.isotope_error != 0]
# df = df[df.peptide == df.iloc[1]['peptide']]
# print(df.charge)

# N = 9
sns.set()
# for N in range(100, 120):
sigmas = []
for peptide in np.random.choice(df.peptide.unique(), 5):
    # for peptide in df.peptide.unique()[5:10]:
    # for peptide in ['LQSRPAAPPAPGPGQLTLR']:

    scans = df[df.peptide == peptide]
    qq = quant[quant.scannr.isin(scans["scannr"])]

    N = 250
    # Create RT grid for binning/KDE
    grid_x = np.linspace(qq.rt.min() - 0.25, qq.rt.max() + 0.25, N)
    x_step = grid_x[1] - grid_x[0]

    # Bin intensities
    apex = np.zeros(N)
    for tup in qq.itertuples():
        i = int((tup.rt - grid_x[0]) / x_step)
        apex[i] = max(apex[i], tup.intensity)

    density = gaussian_kde(qq.rt, bw_method=0.1).pdf(grid_x)
    density = density / density.max()
    apex = apex / apex.max()

    f = norm(0, 1).pdf(np.linspace(-1, 1, 10))
    f /= f.sum()
    smoothed_vals = np.convolve(apex, f, mode="same")

    # Scale smoothed MS1 apex by ion density
    S = smoothed_vals * density

    # Peak picking
    S_max = S.max()
    S_peak = np.where(S == S_max)[0][0]
    S_max = smoothed_vals[S_peak]
    S_half = S_max / 2
    S_shol = [0, len(S)]
    i = S_peak
    while i > 0:
        if smoothed_vals[i] < S_half:
            S_shol[0] = i
            break
        i -= 1
    i = S_peak
    while i < len(S):
        if smoothed_vals[i] < S_half:
            S_shol[1] = i
            break
        i += 1

    # Cosine similarity between apex intensity & ion density
    dist = cosine(smoothed_vals, density)

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, sharex=True)

    # Plot 2D KDE of RT by mass for all isotopologues/charge states
    sns.kdeplot(data=qq, x="rt", y="mass", fill=True, common_norm=False, ax=ax0)
    ax1.plot(grid_x, smoothed_vals, label="apex intensity")
    ax1.plot(grid_x, S, label="smoothed area")
    ax1.hlines(0.5 * S_max, *grid_x[S_shol], colors="black", label="FWHM")
    ax1.set_ylabel("% of max intensity")
    ax1.set_title(dist)
    ax1.legend()

    ax2.plot(grid_x, density, color="red")
    ax2.vlines(scans["rt"], *ax2.get_ylim())
    ax2.set_ylabel("Density")
    plt.suptitle(f"{peptide}")
    plt.xlabel("RT")
    plt.show()
