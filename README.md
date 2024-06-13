# NeoChemSynthWave: Chemistry
![Static Badge](https://img.shields.io/badge/ncsw__chemistry-v.2024.06.1-%23E0457B?logo=github&style=flat)
![Static Badge](https://img.shields.io/badge/Institute%20of%20Science%20Tokyo-%231C3177?style=flat)
![Static Badge](https://img.shields.io/badge/Elix%2C%20Inc.-%235EB6B3?style=flat)


## Introduction
Over the last decade, computer-assisted chemical synthesis has re-emerged as a heavily researched subject in
Chemoinformatics. Even though the idea of utilizing computers for chemical synthesis has existed for as long as
computers themselves, the high level of reliability and innovation expected of such approaches has been repeatedly
proven difficult to achieve. In recent years, however, utilizing machine learning has proven particularly promising,
with novel approaches emerging frequently.

Therefore, considering the interdisciplinary nature of the subject, providing essential programming utilities to
researchers from various backgrounds is paramount to maximizing success. The main objective of the **NeoChemSynthWave:
Chemistry** project is to achieve that by developing, documenting, and maintaining a comprehensive and easy-to-use
Python package titled [ncsw_chemistry](/ncsw_chemistry).


## Installation
A minimal virtual environment can be created using [Conda](https://docs.conda.io/en/latest) as follows:

```shell
conda env create -f environment.yaml

conda activate ncsw-chemistry
```

The package can be locally installed using [pip](https://pip.pypa.io/en/stable) as follows:

```shell
pip install --no-build-isolation -e .
```


## What's Next?
The following updates are currently planned for version **v.2024.07**:

- [ ] Create the _/documentation_ directory.
- [ ] Create the _/notebooks_ directory.
- [ ] Create the _/scripts_ directory.
- [ ] Create the _/tests_ directory.
- [ ] Publish the package on [PyPI](https://pypi.org).


## License Information
The contents of this repository are published under the [MIT](/LICENSE) license. Please refer to individual
[references](#references) for more details regarding the license information of external resources utilized within this 
repository.


## Contact
If you are interested in contributing to this repository by reporting bugs, suggesting improvements, or submitting
feedback, feel free to use [GitHub Issues](https://github.com/neo-chem-synth-wave/ncsw-chemistry/issues).


## References
**[[1]](https://www.rdkit.org)** **RDKit: Open-source Cheminformatics**: https://www.rdkit.org. Accessed on: June
1st, 2024.

**[[2]](https://doi.org/10.1021/acs.jcim.9b00286)** Coley, C.W., Green, W.H., and Jensen K.F. **RDChiral: An RDKit
Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and Application**. _J. Chem. Inf. Model., 59,
6, 2529-2537, 2019_.
