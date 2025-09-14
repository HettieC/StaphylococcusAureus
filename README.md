# Enzyme-constrained genome-scale metabolic model of _S. aureus_ strain S0385.

This repository contains the julia code to build the ecGSMM, as well as the scripts to reproduce the results in chapter 3 of my doctoral thesis: 'Enzyme-constrained genome-scale metabolic model of _Staphylococcus aureus_ and its implications for drug target finding'. The strain has the NCBI accession number AM990992.1.

## Downloading the model

If only interested in downloading the model to use elsewhere, you can download the file ```data/s_aureus.json```. This is the model as of 14th September 2025.

## How to build the model 

To ensure that reactions are up-to-date with rhea, it is recommended to build the model from scratch. To do so, first clone the repo, and activate a new julia environment. 

To build the model, run:
```
using StaphylococcusAureus
model, reaction_isozymes = build_model()
```

Building the model uses  ```RheaReactions```, a small package for downloading reactions from [rhea](https://www.rhea-db.org/). This package works best once reactions have been cached locally. The first time the model is built, this will cache the reactions. Unfortunately, due to the high number of API requests, this first build often outputs errors. If this occurs, simply re-run ```build_model()``` again until it no longer errors. Each time, more reactions will be cached.

After the model has been built, the files in ```StaphylococcusAureus/scripts``` can be used to reproduce all figures from the thesis.

