""" The ``ncsw_chemistry.utility.reaction`` package ``typing`` module. """

from typing import List, Tuple

from rdkit.Chem.rdchem import Mol


ReactionCompoundsTuple = Tuple[
    List[Tuple[str, Mol, str, Mol]],
    List[Tuple[str, Mol, str, Mol]],
    List[Tuple[str, Mol, str, Mol]],
]
