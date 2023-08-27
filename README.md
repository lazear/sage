<img src="figures/logo.png" width="300">

# Sage: proteomics searching so fast it seems like magic

[![Rust](https://github.com/lazear/sage/actions/workflows/rust.yml/badge.svg)](https://github.com/lazear/sage/actions/workflows/rust.yml) [![Anaconda-Server Badge](https://anaconda.org/bioconda/sage-proteomics/badges/version.svg)](https://anaconda.org/bioconda/sage-proteomics)


For more information please read [the online documentation!](https://sage-docs.vercel.app/docs)


# Introduction
 
Sage is, at it's core, a proteomics database search engine - 
    a tool that transforms raw mass spectra from proteomics experiments into peptide identifications 
    via database searching & spectral matching. 

However, Sage includes a variety of advanced features that make it a one-stop shop: retention time prediction, quantification (both isobaric & LFQ), peptide-spectrum match rescoring, and FDR control. You can directly use results from Sage without needing to use other tools for these tasks.

Additionally, Sage was designed with cloud computing in mind - massively parallel processing and the ability to directly stream compressed mass spectrometry data to/from AWS S3 enables unprecedented search speeds with minimal cost. 

 Sage also runs just as well reading local files from your Mac/PC/Linux device!

## Why use Sage instead of other tools?

Sage is **simple to configure**, **powerful** and **flexible**. 
It also happens to be well-tested, **mind-boggingly fast**, open-source (MIT-licensed) and free.

## Features

- Incredible performance out of the box
- [Effortlessly cross-platform](https://sage-docs.vercel.app/docs/started#download-the-latest-binary-release) (Linux/MacOS/Windows), effortlessly parallel (uses all of your CPU cores)
- [Fragment indexing strategy](https://sage-docs.vercel.app/docs/how_it_works) allows for blazing fast narrow and open searches (> 500 Da precursor tolerance)
- [Isobaric quantification](https://sage-docs.vercel.app/docs/how_it_works#tmt-based) (MS2/MS3-TMT, or custom reporter ions)
- [Label-free quantification](https://sage-docs.vercel.app/docs/how_it_works#label-free): consider all charge states & isotopologues *a la* FlashLFQ
- Capable of searching for [chimeric/co-fragmenting spectra](https://sage-docs.vercel.app/docs/configuration/additional)
- Wide-window (dynamic precursor tolerance) search mode - [enables WWA/PRM/DIA searches](https://sage-docs.vercel.app/docs/configuration/tolerance#wide-window-mode)
- Retention time prediction models fit to each LC/MS run
- [PSM rescoring](https://sage-docs.vercel.app/docs/how_it_works#machine-learning-for-psm-rescoring) using built-in linear discriminant analysis (LDA)
- PEP calculation using a non-parametric model (KDE)
- FDR calculation using target-decoy competition and picked-peptide & picked-protein approaches
- Percolator/Mokapot [compatible output](https://sage-docs.vercel.app/docs/configuration#env)
- Configuration by [JSON file](https://sage-docs.vercel.app/docs/configuration#file)
- Built-in support for reading gzipped-mzML files
- Support for reading/writing directly from [AWS S3](https://sage-docs.vercel.app/docs/configuration/aws)


<!-- <img src="figures/TMT_Panel.png" width="800"> -->

Check out the [blog post introducing Sage](https://lazear.github.io/sage/) for more information and full benchmarks!