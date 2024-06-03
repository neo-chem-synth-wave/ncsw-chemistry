""" The ``ncsw_chemistry.reaction`` package ``typing`` module. """

from typing import List, Optional, Tuple

from rdkit.Chem.rdchem import Mol


ReactionCompoundsTuple = Tuple[
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
]
