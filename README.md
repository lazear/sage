# Proteomics Database in a Weekend!
<p align="center">
<img src="https://upload.wikimedia.org/wikipedia/commons/7/70/Carina_Nebula_by_Harel_Boren_%28151851961%2C_modified%29.jpg" width="450px" />
</p>

I wanted to see how far I could take a proteomics search engine in ~1000 lines of code, and spending only one weekend on it. 

I was inspired by the elegant data structure discussed in the [MSFragger paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5409104/), and decided to implement a version of it in Rust
 
### Features

- Search by fragment, filter by precursor: memory layout using sorted tables for blazing fast performance (lots of binary searching!)
- Internal support for static/variable mods
- X!Tandem hyperscore function
- Interal q-value/FDR calculation using a target-decoy competition approach
- No unsafe
- Percolar output

### Limitations

- Only uses MS2 files
- Only percolator PIN output
- Lots of hardcoded parameters (however, you can supply static mods/missed cleavages, etc)
- Probably has some calculation errors :)

### Performance

Raw data extracted using RawConverter to MS2 format
Test data from [PXD001928](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001928) (Yeast_XP_Tricine_trypsin_147.ms2)

On my 6-year old 4-core desktop:

In silico digestion and generation of sorted fragment tables: ~12 s

Wide search with 300ppm fragment window and 5 Th precursor window:
- ~2300 scans/sec
- Total run time: 27.5s for a single raw file, < 2 Gb of peak memory usage
- Percolator/Mokapot results: 13,784 PSMs with q <= 0.05

Narrow search - 10ppm fragment window and 1 Th precursor window:
- ~19,000 scans/sec
- Total run time: 17.2s for a single raw file, < 2 Gb of peak memory usage
- Percolator/Mokapot results: 452 PSMs with q <= 0.05
