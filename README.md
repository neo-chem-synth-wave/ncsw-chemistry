# NeoChemSynthWave: Chemistry
![Static Badge](https://img.shields.io/badge/ncsw__chemistry-2024.7.1-%23556DC8?logo=github&style=flat)
![Static Badge](https://img.shields.io/badge/Elix%2C%20Inc.-%235EB6B3?style=flat)
![Static Badge](https://img.shields.io/badge/Institute%20of%20Science%20Tokyo-%231C3177?style=flat)

Welcome to the **NeoChemSynthWave: Chemistry** project !!!

Over the last decade, computer-assisted chemical synthesis has re-emerged as a heavily researched subject in
Chemoinformatics. Even though the idea of utilizing computers to assist chemical synthesis has existed for nearly as
long as computers themselves, the expected blend of reliability and innovation has repeatedly been proven difficult to
achieve. Nevertheless, recent machine learning approaches have exhibited the potential to address these shortcomings.
Considering the interdisciplinary nature of such approaches, providing essential programming utilities to researchers
from various backgrounds is paramount to maximizing success. The main objective of the **NeoChemSynthWave: Chemistry**
project is to provide chemistry programming utilities that are easy to understand and easy to utilize regardless of
background and skill level.


## Installation
A standalone environment can be created using the [git](https://git-scm.com) and [conda](https://conda.io) commands as
follows:

```shell
git clone https://github.com/neo-chem-synth-wave/ncsw-chemistry.git

cd ncsw-chemistry

conda env create -f environment.yaml

conda activate ncsw-chemistry-env
```

The [ncsw_chemistry](/ncsw_chemistry) package can be installed using the [pip](https://pip.pypa.io) command as follows:

```shell
pip install --no-build-isolation -e .
```


## Utilization
The purpose of the [scripts](/scripts) directory is to illustrate how to utilize the [ncsw_chemistry](/ncsw_chemistry)
package in the following scenarios:

1. [Extraction of Chemical Reaction Retro Templates]()
2. [Application of Chemical Reaction Retro Templates]()
3. [Extraction of Chemical Reaction Reactive Sites and Synthons]()


### Scenario 1: Extraction of Chemical Reaction Retro Templates
The chemical reaction retro templates can be extracted as follows:

```shell
python scripts/extract_reaction_retro_template.py \
  --mapped_reactant_compounds_smiles "[CH3:9][CH2:8][O:7][C:6](=[O:10])[CH2:5][N:4]1[CH2:3][CH2:2][C:1](=O)[CH2:12][CH2:11]1.[CH2:18]([N:19]1[CH2:20][CH2:21][NH:22][CH2:23][CH2:24]1)[c:17]1[cH:16][cH:15][cH:14][cH:26][cH:25]1" \
  --mapped_product_compound_smiles "[CH3:9][CH2:8][O:7][C:6](=[O:10])[CH2:5][N:4]1[CH2:3][CH2:2][CH:1]([CH2:12][CH2:11]1)[N:22]1[CH2:21][CH2:20][N:19]([CH2:18][c:17]2[cH:16][cH:15][cH:14][cH:26][cH:25]2)[CH2:24][CH2:23]1"
```


### Scenario 2: Application of Chemical Reaction Retro Templates
The chemical reaction retro templates can be applied as follows:

```shell
python scripts/apply_reaction_retro_template.py \
  --reaction_retro_template_smarts "[#16:4]=[N;H0;D2;+0:5]-[C;H0;D3;+0:1](-[C;D1;H3:2])=[O;D1;H0:3]>>Cl-[C;H0;D3;+0:1](-[C;D1;H3:2])=[O;D1;H0:3].[#16:4]=[NH;D1;+0:5]" \
  --compound_smiles "CCOC(=O)CN1CCC(CC1)N1CCN(Cc2ccccc2)CC1"
```


### Scenario 3: Extraction of Chemical Reaction Reactive Sites and Synthons
The chemical reaction reactive sites and synthons can be extracted as follows:

```shell
python scripts/apply_reaction_retro_template.py \
  --mapped_reaction_smiles "[CH3:9][CH2:8][O:7][C:6](=[O:10])[CH2:5][N:4]1[CH2:3][CH2:2][C:1](=O)[CH2:12][CH2:11]1.[CH2:18]([N:19]1[CH2:20][CH2:21][NH:22][CH2:23][CH2:24]1)[c:17]1[cH:16][cH:15][cH:14][cH:26][cH:25]1>>[CH3:9][CH2:8][O:7][C:6](=[O:10])[CH2:5][N:4]1[CH2:3][CH2:2][CH:1]([CH2:12][CH2:11]1)[N:22]1[CH2:21][CH2:20][N:19]([CH2:18][c:17]2[cH:16][cH:15][cH:14][cH:26][cH:25]2)[CH2:24][CH2:23]1"
```


## License Information
The contents of this repository are published under the [MIT](/LICENSE) license. Please refer to individual references
for more details regarding the license information of external resources utilized within this repository.


## Contact
If you are interested in contributing to this repository by reporting bugs, suggesting improvements, or submitting
feedback, feel free to use [GitHub Issues](https://github.com/neo-chem-synth-wave/ncsw-chemistry/issues).
