from pyteomics import mzml, pylab_aux
import spectrum_utils.spectrum as sus
import matplotlib.pyplot as plt
from matplotlib import colors, patches
import pandas as pd
import numpy as np

sus.static_modification('C', 57.02146)
sus.static_modification('K', 229.1629)

scans = list(mzml.read('../Yeast_TMT_TKO_1_slice.mzML'))
# print(scans)
ms2 = [s for s in scans if 'scan=16260' in s['id']][0]
ms3 = [s for s in scans if 'scan=16269' in s['id']][0]
print(ms2)

mz = ms2['m/z array']
intensity = ms2['intensity array']

peptide = 'DAIITAEK'
# peptide = 'IGTGIGTGIGSR'

precursor = ms2['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]

selectedSPS = ms3['precursorList']['precursor']

sps = [
    p['selectedIonList']['selectedIon'][0]['selected ion m/z'] for p in selectedSPS[:-1]
]
print(sps)

fig, ax = plt.subplots()

pylab_aux.annotate_spectrum(
    ms2,
    peptide=peptide,
    # max_charge=precursor['charge state'] - 1,
    max_charge = 2,
    precursor_charge=precursor['charge state'],
    scaling="root",
    max_num_peaks=150,
    backend="spectrum_utils",
    ftol=1.0,
    modifications= {
        'N-term': 229.1629
    },
    ax=ax
)

sps_purity = 0.60024494

fragments = [974.601, 943.522, 903.5639, 790.39, 743.5, 677.34]
for sps_ion in sps:
    color='grey'
    for i in fragments:
        if abs(i - sps_ion) <= 1:
            color='green'
    rect = patches.Rectangle((sps_ion - 1, 0), 2, plt.ylim()[1], color=color, alpha=0.5)
    ax.add_patch(rect)

plt.suptitle(f"BTX420 Yeast TMT TKO: MS2 Scan 16260 - {peptide}. SPS purity score:  {sps_purity}")
plt.show(block=True)


