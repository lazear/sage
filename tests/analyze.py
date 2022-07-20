import os
import ppx
import json
from typing import Optional
from pyteomics import mzxml, pylab_aux
import spectrum_utils.spectrum as sus
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sus.static_modification('C', 57.02146)


class Pipeline:
    def __init__(self):
        if not os.path.exists("pxd/b1906_293T_proteinID_01A_QE3_122212.mzXML"):
            self.fetch()
        
        if not os.path.exists("b1906_293T_proteinID_01A_QE3_122212.ms2"):
            self.generate_ms2()

    def fetch(self):
        proj = ppx.find_project("PXD001468", local="pxd")
        proj.download("b1906_293T_proteinID_01A_QE3_122212.mzXML")

    def generate_ms2(self):
        with open("b1906_293T_proteinID_01A_QE3_122212.ms2", "w") as f:
            for scan in mzxml.read(
                "pxd/b1906_293T_proteinID_01A_QE3_122212.mzXML", huge_tree=True
            ):
                if scan["msLevel"] == 2:
                    scan_num = scan["num"]
                    rt = scan["retentionTime"]
                    precursor_mz = scan["precursorMz"][0]["precursorMz"]
                    precursor_charge = scan["precursorMz"][0]["precursorCharge"]
                    precursor_mz = (
                        (precursor_mz * precursor_charge)
                        - (precursor_charge * 1.007_276_4)
                        + 1.007_276_4
                    )
                    headers = f"S\t{scan_num}\nZ\t{precursor_charge}\t{precursor_mz}\nI\tRetTime {rt}\n"
                    f.write(headers)

                    for mz, intensity in zip(
                        scan["m/z array"].astype(np.float32),
                        scan["intensity array"].astype(int),
                    ):
                        f.write(f"{mz:0.4f}\t{intensity}\n")

class PostAnalysis:
    def __init__(self):
        self.results = pd.read_csv(
            "b1906_293T_proteinID_01A_QE3_122212.ms2.carina.pin", sep="\t"
        ).set_index("scannr")

    def dump_scan(self, peptide: str, scan):
        rt = scan["retentionTime"]
        precursor_mz = scan["precursorMz"][0]["precursorMz"]
        precursor_charge = scan["precursorMz"][0]["precursorCharge"]
        precursor_mz = (
            (precursor_mz * precursor_charge)
            - (precursor_charge * 1.007_276_4)
            + 1.007_276_4
        )
        df = pd.DataFrame(dict(mz=scan['m/z array'], intensity=scan['intensity array']))
        df.loc[:,'precursor_charge'] = precursor_charge
        df.loc[:,'precursor_mz'] = precursor_mz
        df.loc[:,'scan'] = scan['num']
        df.loc[:,'rt'] = rt
        print(df)
        df.to_csv(f"{peptide}.csv")



    def plot(self, specific_scan: Optional[int]):
        for scan in mzxml.read("pxd/b1906_293T_proteinID_01A_QE3_122212.mzXML", huge_tree=True):
            if scan["msLevel"] == 2:
                sn = int(scan["num"])
                if ss := specific_scan:
                    if ss != sn:
                        continue
                if sn in self.results.index:
                    row = self.results.loc[
                        sn,
                    ]
                    peptide = ''.join([c for c in row['peptide'] if c.isalpha()])
                    pylab_aux.annotate_spectrum(
                        scan,
                        peptide=peptide,
                        max_charge=row["charge"],
                        precursor_charge=row["charge"],
                        scaling="root",
                        max_num_peaks=100,
                        backend="spectrum_utils",
                    )
                    plt.suptitle(peptide)
                    plt.show()

                    if specific_scan:
                        self.dump_scan(peptide, scan)

_ = Pipeline()
p = PostAnalysis()
p.plot(30069)