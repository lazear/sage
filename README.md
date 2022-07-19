# Proteomics Search Engine in a Weekend!
<p align="center">
<img src="https://upload.wikimedia.org/wikipedia/commons/7/70/Carina_Nebula_by_Harel_Boren_%28151851961%2C_modified%29.jpg" width="450px" />
</p>

I wanted to see how far I could take a proteomics search engine in ~1000 lines of code, and spending a little more than a weekend on it. 

I was inspired by the elegant data structure discussed in the [MSFragger paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5409104/), and decided to implement a version of it in Rust - with great results.

Carina has excellent performance characteristics (>2x faster and >2x less memory usage than MSFragger), but does not sacrifice code quality or size to do so!
 
## Features

- Search by fragment, filter by precursor: memory layout using sorted tables for blazing fast performance (lots of binary searching!)
- Effortlessly cross-platform
- Configuration by JSON files
- Internal support for static/variable mods
- X!Tandem hyperscore function
- Internal q-value/FDR calculation using a target-decoy competition approach
- Unsafe free
- Percolator/mokapot compatible output
- 

## Limitations

- Only uses MS2 files
- Only percolator PIN output
- Only outputs 1 protein ID even if the peptide is shared by multiple proteins :)
- Probably has some calculation errors :)

# Usage 

Carina takes a single command line argument: a path to a JSON-encoded parameter file (see below). A new file (`results.json`) will be created that details input/output paths and all search parameters used for the search

Example usage: `carina tmt.json`

# Performance

To benchmark search performance versus MSFragger (closed source) and Comet (open source), I downloaded data from the paper [Benchmarking the Orbitrap Tribrid Eclipse for Next Generation
Multiplexed Proteomics](https://pubs.acs.org/doi/10.1021/acs.analchem.9b05685?goto=supporting-info).

Data repository: [PXD016766](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD016766)

***Carina* has good peptide identity overlap**

<img src="figures/TMT_IDs.png" width="600">

**And excellent performance relative to other search engines**

<img src="figures/TMT_Panel.png" width="800">

Performance results: (Intel i7-12700KF + 32GB ram)

- Active scanning: ~20,000 scans/s (can be tuned to use more ram and go up to 5x faster!)
- Amortized scanning: ~7,000 scans/s (most time used on fragment indexing, IO)


### Search methods

- MS2 files generated using the [ProteoWizard MSConvert tool](http://www.proteowizard.org/download.html)
- MSFragger and Comet were configured with analogous parameters (50ppm precursor tolerance, 10ppm fragment tolerance - or for Comet setting `fragment_bin_tol` to 0.02 Da).
- [Mokapot](https://github.com/wfondrie/mokapot) was then used to refine FDR for all search results

Carina search settings file:
```json
{
  "database": {
    "bucket_size": 8192,
    "fragment_min_mz": 75.0,
    "fragment_max_mz": 4000.0,
    "peptide_min_len": 5,
    "peptide_max_len": 50,
    "decoy": true,
    "missed_cleavages": 1,
    "n_term_mod": 229.1629,
    "static_mods": {
      "K": 229.1629,
      "C": 57.0215
    },
    "fasta": "UP000005640_9606.fasta"
  },
  "precursor_tol": {
    "ppm": 50.0
  },
  "fragment_tol": {
    "ppm": 10.0 
  },
  "report_psms": 1,
  "ms2_paths": [
    "./tmt_analysis/raw/dq_00082_11cell_90min_hrMS2_A1.ms2",
    "./tmt_analysis/raw/dq_00083_11cell_90min_hrMS2_A3.ms2",
    "./tmt_analysis/raw/dq_00084_11cell_90min_hrMS2_A5.ms2",
    "./tmt_analysis/raw/dq_00085_11cell_90min_hrMS2_A7.ms2",
    "./tmt_analysis/raw/dq_00086_11cell_90min_hrMS2_A9.ms2",
    "./tmt_analysis/raw/dq_00087_11cell_90min_hrMS2_A11.ms2",
    "./tmt_analysis/raw/dq_00088_11cell_90min_hrMS2_B1.ms2",
    "./tmt_analysis/raw/dq_00089_11cell_90min_hrMS2_B3.ms2",
    "./tmt_analysis/raw/dq_00090_11cell_90min_hrMS2_B5.ms2",
    "./tmt_analysis/raw/dq_00091_11cell_90min_hrMS2_B7.ms2",
    "./tmt_analysis/raw/dq_00092_11cell_90min_hrMS2_B9.ms2",
    "./tmt_analysis/raw/dq_00093_11cell_90min_hrMS2_B11.ms2"
  ]
}
```

