# NeoChemSynthWave: Chemistry
![Static Badge](https://img.shields.io/badge/ncsw__chemistry-v.2024.06.1-%23E0457B?logo=github&style=flat)
![Static Badge](https://img.shields.io/badge/Institute%20of%20Science%20Tokyo-%231C3177)
![Static Badge](https://img.shields.io/badge/Elix%2C%20Inc.-%235EB6B3)


## Introduction
Over the last decade, computer-assisted chemical synthesis has re-emerged as a heavily researched subject in
Chemoinformatics. Even though the idea of utilizing computers for chemical synthesis has existed for as long as
computers themselves, the high level of reliability and innovation expected of such approaches has been repeatedly
proven difficult to achieve. In recent years, however, utilizing machine learning has proven particularly promising,
with novel approaches emerging frequently. Therefore, considering the interdisciplinary nature of the subject,
providing essential programming utilities to researchers from various backgrounds is paramount to maximizing success.
The main objective of the **NeoChemSynthWave: Chemistry** project is to achieve that by developing and maintaining a
comprehensive set of the most relevant computer-assisted chemical synthesis utilities within a single Python package.


## Installation
### Virtual Environment
A virtual environment for the [**ncsw_chemistry**](/ncsw_chemistry) package can be set up using
[**Conda**](https://docs.conda.io/en/latest) and [**pip**](https://pip.pypa.io/en/stable) as follows:

```shell
conda env create -f environment.yaml

conda activate ncsw-chemistry

pip install --no-build-isolation -e .
```


### Utilizing Code Snippets of Interest
The [**ncsw_chemistry**](/ncsw_chemistry) package is an organized collection of utility classes, which means that the
static methods are applicable outside the class context as long as all dependencies are satisfied. Therefore, utilizing
only code snippets of interest is an appropriate strategy to avoid the installation altogether.


## References
1. **RDKit: Open-source Cheminformatics**: https://www.rdkit.org. Accessed on: June 5th, 2024.
2. Coley, C.W., Green, W.H., and Jensen K.F. **RDChiral: An RDKit Wrapper for Handling Stereochemistry in 
   Retrosynthetic Template Extraction and Application**. _J. Chem. Inf. Model., 59, 6, 2529-2537, 2019_. DOI:
   https://doi.org/10.1021/acs.jcim.9b00286.


## License Information
The contents of this repository are published under the [**MIT**](/LICENSE) license. Please refer to the individual
sources for more details regarding the license information of external resources utilized within this repository.


## Contact
If you are interested in contributing to this repository by reporting bugs, submitting feedback or anything else that
might be beneficial, please feel free to do so via
[**GitHub Issues**](https://github.com/neo-chem-synth-wave/ncsw-chemistry/issues).
