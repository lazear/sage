# Sage configuration documentation

This documentation covers the parameters in the JSON configuration file for the proteomics search engine. The configuration file contains information about the search engine's settings, including database, enzyme, modifications, and other settings. For a complete example of a configuration file, please see the bottom of [README.md](README.md#example-configuration-file). If any parameters or options are not clear, please open an issue!

## Database

- **bucket_size**: Integer. The number of fragments in each internal mass bucket (default: 8192). Tweaking this parameter can increase search performance for wide precursor or fragment searches.

### Enzyme

The enzyme section contains parameters related to the enzyme used for digestion. The default enzyme is trypsin, with the parameters specified below.

- **missed_cleavages**: Integer. The number of missed cleavages for tryptic digest (default: 1).
- **min_len**: Integer. The minimum amino acid (AA) length of peptides to search (default: 5).
- **max_len**: Integer. The maximum AA length of peptides to search (default: 50).
- **cleave_at**: String. Amino acids to cleave at (default: 'KR').
- **restrict**: Single character string. Do not cleave if this amino acid follows the cleavage site (default: 'P').
- **c_terminal**: Boolean. Cleave at the C-terminus of matching amino acids (default:true).

Example: 
```json
"database": {
  "enzyme": {
    "missed_cleavages": 1,
    "min_len": 5,
    "max_len": 50,
    "cleave_at": "KR",
    "restrict": "P",
    "c_terminal": true
  }
}
```

### Fragment Settings

- **fragment_min_mz**: Float. The minimum mass of fragments to search (default: 150.0).
- **fragment_max_mz**: Float. The maximum mass of fragments to search (default: 2000.0).
- **peptide_min_mass**: Float. The minimum monoisotopic mass of peptides to fragment *in silico* (default: 500.0).
- **peptide_max_mass**: Float. The maximum monoisotopic mass of peptides to fragment *in silico* (default: 5000.0).
- **ion_kinds**: List of strings. Which fragment ions to produce? Allowed values: "a", "b", "c", "x", "y", "z". (default: ["b", "y"])
- **min_ion_index**: Integer. Do not generate b1/bN/y1/yN ions for preliminary searching if `min_ion_index = N`. Does not affect full scoring of PSMs (default: 2).

Example:
```json
"database": {
  "fragment_min_mz": 150.0,
  "fragment_max_mz": 2000.0,
  "peptide_min_mass": 500.0,
  "peptide_max_mass": 5000.0,
  "ion_kinds": ["b", "y"],
  "min_ion_index": 2
}
```

### Modifications

#### Static Modifications

- **static_mods**: Dictionary with characters as keys and floats as values. Represents static modifications applied to amino acids or termini (default: {}). Static modifications are applied after variable modifications
  - Example: Apply a static modification of 304.207 to the N-terminus of the peptide and lysine, and 57.0215 to cysteine.
    ```json
    "database": {
      "static_mods": {
        "^": 304.207,
        "K": 304.207,
        "C": 57.0215
      }
    }
    ```

#### Variable Modifications

- **max_variable_mods**: Integer. Limit k-combinations of variable modifications (default: 2).
- **variable_mods**: Dictionary with characters as keys and list of floats (or single floats) as values. Represents variable modifications applied to amino acids or termini (default: {}).
  - Example: Apply a variable modification of 15.9949 to methionine, 49.2022 to the C-terminus of the peptide, 42.0 to the N-terminus of the protein, and 111.0 to the C-terminus of the protein.
    ```jsonc
    "database": {
      "variable_mods": {
        "M": [15.9949], 
        "^Q": [-17.026549],
        "^E": [-18.010565], // Applied to N-terminal glutamic acid
        "$": [49.2022],     // Applied to peptide C-terminus
        "[": 42.0,          // Applied to protein N-terminus
        "]": 111.0          // Applied to protein C-terminus
      }
    }
    ```
  - Syntax:
    "^X": Modification to be applied to amino acid X if it appears at the N-terminus of a peptide
    "$X": Modification to be applied to amino acid X if it appears at the C-terminus of a peptide
    "[X": Modification to be applied to amino acid X if it appears at the N-terminus of a protein
    "]X": Modification to be applied to amino acid X if it appears at the C-terminus of a protein

### Decoys

- **decoy_tag**: String. The tag used to identify decoy entries in the FASTA database (default: "rev_").
- **generate_decoys**: Boolean. If true, ignore decoys in the FASTA database matching `decoy_tag`, and generate internally reversed peptides (default: false).

### FASTA

- **fasta**: String. The path to the FASTA file, either a local path or s3 object URI.

## Quantification

The quant section is optional and should be specified only if TMT or LFQ is used.


- **tmt**: String. One of "Tmt6", "Tmt10", "Tmt11", "Tmt16", or "Tmt18" (default: null).
- **tmt_settings**: Object containing TMT-specific settings.
  - **level**: Integer. The MS-level to perform TMT quantification on (default: 3).
  - **sn**: Boolean. Use Signal/Noise instead of intensity for TMT quantification. Requires noise values in mzML (default: false).
- **lfq**: Boolean. Perform label-free quantification (default: null).
- **lfq_settings**: Object containing LFQ-specific settings.
  - **peak_scoring**: String. The method used for scoring peaks in LFQ, one of: "Hybrid", "RetentionTime", "SpectralAngle" (default: "Hybrid").
  - **integration**: String. The method used for integrating peak intensities, either "Sum" or "Max" (default: "Sum").
  - **spectral_angle**: Float. Threshold for the spectral angle similarity measure, ranging from 0 to 1 (default: 0.7).
  - **ppm_tolerance**: Float. Tolerance for matching MS1 ions in parts per million (default: 5.0).

Example: 
```json
 "quant": {
    "tmt": "Tmt16",
    "tmt_settings": {
      "level": 3,
      "sn": false
    },
    "lfq": true,
    "lfq_settings": {
      "peak_scoring": "Hybrid",
      "integration": "Sum",
      "spectral_angle": 0.7,
      "ppm_tolerance": 5.0
    }
  }
```


## Precursor Tolerance

- **precursor_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-500, 100] in daltons.
    ```json
    "precursor_tol": {
      "da": [-500, 100]
    }
    ```

## Fragment Tolerance

- **fragment_tol**: Dictionary with either "ppm" or "da" as keys, and lists of two integers as values (default: {}).
  - Example: Tolerance of [-10, 10] in parts per million.
    ```json
    "fragment_tol": {
      "ppm": [-10, 10]
    }
    ```

## Isotope Errors

- **isotope_errors**: List of two integers. The C13 isotopic envelope to consider for precursor (default: [0, 0]).
  - Example: Consider -1 and up to +3 C13 isotopes (-1/0/1/2/3).
    ```json
    "isotope_errors": [-1, 3]
    ```

**NOTE**: Searching with isotope errors is slower than searching with a wider precursor tolerance that encompasses the isotope errors, e.g. `"da": [-3.5, 1.25]`. Using the wider precursor tolerance will generally increase the number of confidently identified PSMs as well.

## Other Settings

Note on the settings below:

`predict_rt` is incompatible with `quant.lfq = true`. Setting `quant.lfq = true` will automatically turn on global retention time alignment and prediction, which are crucial for accurate direct ion current extraction.

- **deisotope**: Boolean. Perform deisotoping and charge state deconvolution on MS2 spectra (default: false). Recommended for high-resolution MS2 scans. This setting may interfere with TMT-MS2 quantification, use at your own risk.
- **chimera**: Boolean. Search for chimeric/co-fragmenting PSMs (default: false).
- **wide_window**: Boolean. Ignore `precursor_tol` and search spectra in wide-window/dynamic precursor tolerance mode (default: false).
- **predict_rt**: Boolean. Use retention time prediction model as a feature for LDA (default: false).
- **min_peaks**: Integer. Only process MS2 spectra with at least N peaks (default: 15).
- **max_peaks**: Integer. Take the top N most intense MS2 peaks to search (default: 150).
- **min_matched_peaks**: Integer. The minimum number of matched b+y ions to use for reporting PSMs (default: 4).
- **max_fragment_charge**: Integer. The maximum fragment ion charge states to consider (default: null - use precursor z-1).
- **report_psms**: Integer. The number of PSMs to report for each spectrum. Higher values might disrupt LDA (default: 1).
- **parallel**: Boolean. Parse and search files in parallel. For large numbers of files or low RAM, setting this to false can reduce memory usage at the cost of running slower (default: true).

## mzML Paths

- **mzml_paths**: List of strings. The paths to mzML (or gzipped-mzML) files for search. Paths are either local, or point to an S3 object. Files ended in ".gz" or ".gzip" are inferred to be compressed.
  - Example:
    ```json
    "mzml_paths": [
      "local/path.mzML",
      "s3://my-mass-spec-data/PXD0000001/foo.mzML.gz"
    ]
    ```
  
## Output directory:

- **output_directory**: Local directory, or S3 location where output files will be written. If the local directory does not already exist, it will be created. Write permissions are required for the directory or S3 path.
  - Possible output files are: "results.json", "results.sage.tsv", "lfq.tsv", and "tmt.tsv"
  - Example:
  ```json
  "output_directory": "s3://my-mass-spec-results/PXD003881/"
  ```

# Interpreting Sage Output

The "results.sage.tsv" file contains the following columns (headers):

- `peptide`: Peptide sequence, including modifications (e.g., NC\[+57.021\]HKGSFK).
- `proteins`: Proteins containing the peptide sequence.
- `num_proteins`: Number of proteins assigned to the peptide sequence.
- `filename`: File containing this PSM
- `scannr`: Spectrum identifier from mzML file.
- `rank`: Rank of the PSM. If `report_psms > 1`, then the best match will have rank = 1, the second best match will have rank = 2, etc. 
- `label`: Target/Decoy label (-1: decoy, 1: target).
- `expmass`: Experimental mass of the peptide.
- `calcmass`: Calculated mass of the peptide.
- `charge`: Reported precursor charge.
- `pepide_len`: Length of the peptide sequence.
- `missed_cleavages`: Number of missed cleavages.
- `isotope_error`: C13 isotope error.
- `precursor_ppm`: Difference between experimental mass and calculated mass, reported in parts-per-million.
- `fragment_ppm`: Average parts-per-million (delta mass) for matched fragment ions compared to theoretical ions.
- `hyperscore`: X!Tandem hyperscore for the PSM.
- `delta_next`: Difference between the hyperscore of this candidate and the next best candidate.
- `delta_bext`: Difference between the hyperscore of the best candidate (rank=1) and this candidate.
- `rt`: Retention time.
- `aligned_rt`: Globally aligned retention time.
- `predicted_rt`: Predicted retention time, if enabled.
- `delta_rt_model`: Difference between predicted and observed retention time.
- `matched_peaks`: Number of matched theoretical fragment ions.
- `longest_b`: Longest b-ion series.
- `longest_y`: Longest y-ion series.
- `longest_y_pct`: Longest y-ion series, divided by peptide length (as a percentage).
- `matched_intensity_pct`: Fraction of MS2 intensity explained by matched b- and y-ions (as a percentage of total MS2 intensity for this spectrum).
- `scored_candidates`: Number of scored candidates for this spectrum.
- `poisson`: Probability of matching exactly N peaks across all candidates (Pr(x=k)).
- `sage_discriminant_score`: Combined score from linear discriminant analysis, used for FDR (False Discovery Rate) calculation.
- `posterior_error`: Posterior error probability for this PSM / local FDR.
- `spectrum_fdr`: Assigned spectrum-level q-value.
- `peptide_fdr`: Assigned peptide-level q-value.
- `protein_fdr`: Assigned protein-level q-value.
- `ms1_intensity`: Intensity of the MS1 precursor ion
- `ms2_intensity`: Total intensity of MS2 spectrum

These columns provide comprehensive information about each candidate peptide spectrum match (PSM) identified by the Sage search engine, enabling users to assess the quality and characteristics of the results.