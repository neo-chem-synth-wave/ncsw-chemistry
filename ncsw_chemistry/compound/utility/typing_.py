""" The ``ncsw_chemistry.compound.utility`` package ``typing_`` module. """

from typing import FrozenSet, Tuple, Union

CompoundAtomPropertyIDTuple = Tuple[Union[bool, int, float, str], ...]

CompoundBondPropertyIDTuple = Tuple[
    FrozenSet[CompoundAtomPropertyIDTuple],
    Tuple[Union[bool, str], ...]
]

CompoundSubstructurePropertyIDTuple = Tuple[
    FrozenSet[Tuple[CompoundAtomPropertyIDTuple, int]],
    FrozenSet[Tuple[CompoundBondPropertyIDTuple, int]]
]
