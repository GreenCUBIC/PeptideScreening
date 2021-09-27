# Peptide Screening

This repository provides the scripts used to produce the data
used in the paper *Computational design of peptide therapeutics: how useful are sequence-based protein-protein interaction predictors?* by Charih *et al.*

**The experiments were all performed in a Linux environment, with no expectation that this could be reproduced on Windows.**

## Authors

- Francois Charih <francois@charih.ca>
- Kyle K. Biggar <kyle_biggar@carleton.ca>
- James R. Green <jrgreen@sce.carleton.ca>

## Dependencies

The data preparation and analysis is written in Julia, so you will need
[Julia 1.6](https://julialang.org/downloads/) installed on your machine to run this. In addition to Julia, the following libraries are required:

- FASTX
- CSV
- ArgParse
- DataFrames
- Plots

To run this code, you will need to have the following softare installed on your path.

- [CD-HIT](http://cd-hit.org/)
- [blastp (BLAST+ toolkit)](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

The predictors used are available in their respective repositories.

- [D-SCRIPT](https://d-script.readthedocs.io/en/main/) and its [pre-trained weights](http://cb.csail.mit.edu/cb/dscript/data/models/human_v1.sav)
- [PIPR](https://github.com/muhaochen/seq_ppi)
- [SPRINT](https://github.com/lucian-ilie/SPRINT)

## Making predictions

To produce the data in the `data` folder, one can run the script `prepare_datasets.sh`.

The predictors can then be trained according to the instructions on the produced
training data and then used to predict the interactions between the therapeutic
peptides and the human proteome.

## Analysis

The code used to analyze the results is available in `notebooks/Plotting.ipynb`.

## License

Copyright (c) 2021 Carleton University Biomedical Informatics Collaboratory 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
