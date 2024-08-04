import om

ATOM_TYPE_NAMES = {
    om.Atom.SALT: "SALT",
    om.Atom.AIR: "AIR",
    om.Atom.EARTH: "EARTH",
    om.Atom.FIRE: "FIRE",
    om.Atom.WATER: "WATER",
    om.Atom.QUICKSILVER: "QUICKSILVER",
    om.Atom.GOLD: "GOLD",
    om.Atom.SILVER: "SILVER",
    om.Atom.COPPER: "COPPER",
    om.Atom.IRON: "IRON",
    om.Atom.TIN: "TIN",
    om.Atom.LEAD: "LEAD",
    om.Atom.VITAE: "VITAE",
    om.Atom.MORS: "MORS",
    om.Atom.REPETITION_PLACEHOLDER: "REPETITION_PLACEHOLDER",
    om.Atom.QUINTESSENCE: "QUINTESSENCE",
}

ATOM_TYPE_SYMBOLS = {
    om.Atom.SALT: "Sa",
    om.Atom.AIR: "Ai",
    om.Atom.EARTH: "Ea",
    om.Atom.FIRE: "Fi",
    om.Atom.WATER: "Wa",
    om.Atom.QUICKSILVER: "QS",
    om.Atom.GOLD: "Au",
    om.Atom.SILVER: "Ag",
    om.Atom.COPPER: "Cu",
    om.Atom.IRON: "Fe",
    om.Atom.TIN: "Sn",
    om.Atom.LEAD: "Pb",
    om.Atom.VITAE: "Vi",
    om.Atom.MORS: "Mo",
    om.Atom.REPETITION_PLACEHOLDER: "><",
    om.Atom.QUINTESSENCE: "QE",
}


def atom_names(atom_types: list[int]):
    return [ATOM_TYPE_NAMES[atom] for atom in atom_types]


def print_puzzle_details(puzzle: om.Puzzle):
    print(f"puzzle: \"{puzzle.name.decode('utf-8')}\"")
    extra_parts = puzzle.nondef_parts_list()
    print(f"- extra parts: ({len(extra_parts)})")
    for i, part in enumerate(extra_parts):
        print(f"    [{i}]: {part.decode('utf-8')}")
    print(f"- reagents: ({len(puzzle.reagents)})")
    for i, mol in enumerate(puzzle.reagents):
        if len(mol.atoms) == 1:
            print(f"    [{i}]: 1 x {ATOM_TYPE_NAMES[mol.atoms[0].type]}")
        else:
            mol_str_list = molecule_str_list(mol)
            if len(mol_str_list) == 1:
                print(f"    [{i}]: {mol_str_list[0]}")
            else:
                print(f"    [{i}]:")
                for mol_str in mol_str_list:
                    print(mol_str)
    print(f"- products: ({len(puzzle.products)})")
    for i, mol in enumerate(puzzle.products):
        if len(mol.atoms) == 1:
            print(f"    [{i}]: 1 x {ATOM_TYPE_NAMES[mol.atoms[0].type]}")
        else:
            mol_str_list = molecule_str_list(mol)
            if len(mol_str_list) == 1:
                print(f"    [{i}]: {mol_str_list[0]}")
            else:
                print(f"    [{i}]:")
                for mol_str in mol_str_list:
                    print(mol_str)


def get_bond_type(pos1, pos2, bond_positions):
    matched_bonds = [
        bond for bond in bond_positions
        if (pos1 in bond and pos2 in bond)
    ]
    if len(matched_bonds) > 1:
        raise Exception("too many matched bonds")
    elif len(matched_bonds) == 0:
        return None
    else:
        return bond_positions[matched_bonds[0]]


def get_horiz_bond(bond_type):
    if bond_type is None or bond_type == 0:
        return "  "
    if bond_type == om.Bond.NORMAL:
        return "--"
    if bond_type == om.Bond.TRIPLEX:
        return "=="
    if bond_type == om.Bond.TRIPLEX_RED | om.Bond.TRIPLEX_BLACK:
        return "rb"
    if bond_type == om.Bond.TRIPLEX_RED | om.Bond.TRIPLEX_YELLOW:
        return "ry"
    if bond_type == om.Bond.TRIPLEX_BLACK | om.Bond.TRIPLEX_YELLOW:
        return "by"
    if bond_type == om.Bond.TRIPLEX_RED:
        return "r-"
    if bond_type == om.Bond.TRIPLEX_BLACK:
        return "b-"
    if bond_type == om.Bond.TRIPLEX_YELLOW:
        return "y-"
    else:
        return "ER"


def get_vert_bond(bond_type, diag_char):
    if bond_type is None or bond_type == 0:
        return "   "
    if bond_type == om.Bond.NORMAL:
        return " %s " % diag_char
    if bond_type & om.Bond.TRIPLEX:
        r = diag_char if bond_type & om.Bond.TRIPLEX_RED else "-"
        b = diag_char if bond_type & om.Bond.TRIPLEX_BLACK else "-"
        y = diag_char if bond_type & om.Bond.TRIPLEX_YELLOW else "-"
        return r + b + y
    else:
        return "ER"


def molecule_str_list(
    molecule,
    atom_repr=lambda i_at: ATOM_TYPE_SYMBOLS[i_at[1]]
):
    min_q = min(atom.position[0] for atom in molecule.atoms)
    max_q = max(atom.position[0] for atom in molecule.atoms)
    min_r = min(atom.position[1] for atom in molecule.atoms)
    max_r = max(atom.position[1] for atom in molecule.atoms)

    lookup = {atom.position: atom.type for atom in molecule.atoms}
    list_form = [
        [
            lookup[(q, r)] if (q, r) in lookup.keys() else None
            for q in range(min_q, max_q + 1)
        ]
        for r in range(max_r, min_r - 1, -1)
    ]

    atom_positions = {atom.position: atom.type for atom in molecule.atoms}
    bond_positions = {bond.positions: bond.type for bond in molecule.bonds}

    str_list = []

    for r in range(max_r, min_r - 1, -1):
        indent = (r - min_r) * '   '
        str_list.append(
            indent + ''.join(
                [
                    ("(%2s)%2s" if q == 0 and r == 0 else " %2s %2s") % (
                        '  '
                        if (q, r) not in atom_positions else
                        ATOM_TYPE_SYMBOLS[atom_positions[(q, r)]],
                        get_horiz_bond(
                            get_bond_type(
                                (q, r), (q + 1, r),
                                bond_positions
                            )
                        )
                    )
                    for q in range(min_q, max_q + 1)
                ]
            )
        )
        if r != min_r:
            indent = (len(indent) - 1) * ' '
            str_list.append(
                indent + ''.join(
                    '%s%s' % (
                        get_vert_bond(
                            get_bond_type(
                                (q, r), (q, r - 1),
                                bond_positions
                            ), "/"
                        ),
                        get_vert_bond(
                            get_bond_type(
                                (q, r), (q + 1, r - 1),
                                bond_positions
                            ), "\\"
                        )
                    )
                    for q in range(min_q, max_q + 1)
                )
            )
    return str_list
