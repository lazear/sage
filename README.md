# Proteomics Search Engine with Magical Performance

[![Rust](https://github.com/lazear/sage/actions/workflows/rust.yml/badge.svg)](https://github.com/lazear/sage/actions/workflows/rust.yml)

Check out the [blog post](https://lazear.github.io/sage/) for more information and full benchmarks!

<img src="figures/TMT_Panel.png" width="800">

I was inspired by the elegant data structure discussed in the [MSFragger paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5409104/), and decided to implement an (open source) version of it in Rust - with great results.

Sage has excellent performance characteristics (2-5x faster, 2-3x reduction in memory use compared to - the closed source - MSFragger), but does not sacrifice code quality or size to do so!
 
## Features & Anti-Features

- Incredible performance out of the box
- Fragment index search
- Assign chimeric spectra
- Effortlessly cross-platform, effortlessly parallel
- Small and simple codebase
- Configuration by JSON files
- X!Tandem hyperscore function
- Internal q-value/FDR calculation using a target-decoy competition approach
- Percolator/Mokapot compatible output
- Unsafe free
- Only uses mzML files
- Only Percolator PIN output
- Only outputs 1 protein ID even if the peptide is shared by multiple proteins


<img src="figures/chimera_27525.png" width="800">

Assign multiple peptides to complex spectra

# Usage 

Sage takes a single command line argument: a path to a JSON-encoded parameter file (see below). A new file (`results.json`) will be created that details input/output paths and all search parameters used for the search

Example usage: `sage tmt.json`

# Performance

To benchmark search performance versus MSFragger (closed source) and Comet (open source), I downloaded data from the paper [Benchmarking the Orbitrap Tribrid Eclipse for Next Generation
Multiplexed Proteomics](https://pubs.acs.org/doi/10.1021/acs.analchem.9b05685?goto=supporting-info).

Data repository: [PXD016766](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD016766)

**Sage has good peptide identity overlap**

<img src="figures/TMT_IDs.png" width="600">

Performance results: (c5ad.8xlarge, 32 vCPUs)

- ~11 seconds to process 12 files in narrow search, using less than ~3GB of RAM - 3s of this is generating the fragment index, so it's really ~750ms/mzML file
- Active scanning: ~50,000 scans/s for narrow window

### Search methods

- mzML files generated using the [ProteoWizard MSConvert tool](http://www.proteowizard.org/download.html)
- MSFragger and Comet were configured with analogous parameters (+/- 20 ppm precursor tolerance, +/- 10 ppm fragment tolerance - or for Comet setting `fragment_bin_tol` to 0.02 Da).
- [Mokapot](https://github.com/wfondrie/mokapot) was then used to refine FDR for all search results
- Parameter files for all engines can be found in the `figures/benchmark_params` folder!
- All searches for benchmarking were run on an c5ad.8xlarge (32 vCPU, 64 GB RAM, NVMe drives) EC2 instance

Sage search settings file:
```json
{
  "database": {
    "bucket_size": 8192,
    "fragment_min_mz": 75.0,
    "fragment_max_mz": 1500.0,
    "peptide_min_len": 5,
    "peptide_max_len": 50,
    "missed_cleavages": 1,
    "n_term_mod": 229.1629,
    "static_mods": {
      "K": 229.1629,
      "C": 57.0215
    },
    "decoy_prefix": "rev_",
    "fasta": "2022-07-23-decoys-reviewed-UP000005640.fas"
  },
  "precursor_tol": {
    "ppm": [-20, 20]
  },
  "fragment_tol": {
    "ppm": [-10.0, 10.0]
  },
  "isotope_errors": [
    -1,
    3
  ],
  "report_psms": 1,
  "chimera": false,
  "deisotope": true,
  "process_files_parallel": true,
  "mzml_paths": [
    "dq_00082_11cell_90min_hrMS2_A1.mzML",
    "dq_00083_11cell_90min_hrMS2_A3.mzML",
    "dq_00084_11cell_90min_hrMS2_A5.mzML",
    "dq_00085_11cell_90min_hrMS2_A7.mzML",
    "dq_00086_11cell_90min_hrMS2_A9.mzML",
    "dq_00087_11cell_90min_hrMS2_A11.mzML",
    "dq_00088_11cell_90min_hrMS2_B1.mzML",
    "dq_00089_11cell_90min_hrMS2_B3.mzML",
    "dq_00090_11cell_90min_hrMS2_B5.mzML",
    "dq_00091_11cell_90min_hrMS2_B7.mzML",
    "dq_00092_11cell_90min_hrMS2_B9.mzML",
    "dq_00093_11cell_90min_hrMS2_B11.mzML"
  ]
}
```

