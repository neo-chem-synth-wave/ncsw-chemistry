# Overview
**Title:** RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and
Application<br>
**Authors:** Connor W. Coley, William H. Green, and Klavs F. Jensen<br>
**Publication Date:** 2019/06/13<br>
**Publication Links:** [ACS JCIM](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00286)<br>
**Alternative Links:** [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/60c741240f50dbe270395a6d) |
[DSpace@MIT](https://dspace.mit.edu/handle/1721.1/134762)<br>
**Code Links:** [GitHub](https://github.com/connorcoley/rdchiral) | [GitLab](https://gitlab.com/ljn917/rdchiral_cpp)


# Abstract
There is a renewed interest in computer-aided synthesis planning, where the vast majority of approaches require the
application of retrosynthetic reaction templates. Here we introduce RDChiral, an open-source Python wrapper for RDKit
designed to provide consistent handling of stereochemical information in applying retrosynthetic transformations
encoded as SMARTS strings. RDChiral is designed to enforce the introduction, destruction, retention, and inversion of
chiral tetrahedral centers as well as the cis/trans configuration of double bonds. We also introduce an open-source
implementation of a retrosynthetic template extraction algorithm to generate SMARTS patterns from atom-mapped reaction
SMILES strings. In this application note, we describe the implementation of these two pieces of code and illustrate
their use through many examples.


# Citation
```
@article{
    Coley2019,
    author = {Connor W. Coley and William H. Green and Klavs F. Jensen},
    title = {RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and Application},
    journal = {Journal of Chemical Information and Modeling},
    volume = {59},
    number = {6},
    pages = {2529-2537},
    year = {2019},
    doi = {https://doi.org/10.1021/acs.jcim.9b00286}
}
```
