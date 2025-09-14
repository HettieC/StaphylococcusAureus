# Enzyme-constrained genome-scale metabolic model of S. aureus strain S0385.

This repository contains the julia code to build the ecGSMM, as well as the scripts to reproduce the results in chapter 3 of my doctoral thesis: 'Enzyme-constrained genome-scale metabolic model of _Staphylococcus aureus_ and its implications for drug target finding'. The strain has the NCBI accession number AM990992.1.

## How to build the model 

Clone the repo, and activate a new environment. 

To build the model, run:
```
using StaphylococcusAureus
model, reaction_isozymes = build_model()
```

From there, the files in ```StaphylococcusAureus/scripts``` can be used to reproduce figures from the thesis.

