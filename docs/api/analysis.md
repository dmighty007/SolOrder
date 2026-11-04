# Analysis classes

Generated from the docstrings. Conventions, units and validation are described
in [Core parameters](../order_parameters.md) and
[Extended parameters](../extended_parameters.md).

Every analysis class takes an MDAnalysis `Universe`, copies one frame, and
returns one value (or one row) per selected oxygen, in the order of
`atom_indices`. Move to another frame with `refresh(frame)`.

::: SolOrder.OrderParameter
    options:
      inherited_members: true
      filters: ["!^_", "!^single", "!^(dist|maping|neighbours|get_atoms|calculate_dihedral|arr)$"]

::: SolOrder.StructuralOrder
    options:
      members: false

::: SolOrder.eFingerprint
    options:
      members: ["Entropy", "localEntropy", "refresh"]
