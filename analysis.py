import om
import puzzledb

ELEMENTAL_ATOM_TYPES = {
    om.Atom.AIR,
    om.Atom.EARTH,
    om.Atom.FIRE,
    om.Atom.WATER,
}

SIMPLE_ELEMENTAL_ATOM_TYPES = {
    om.Atom.SALT,
    om.Atom.AIR,
    om.Atom.EARTH,
    om.Atom.FIRE,
    om.Atom.WATER,
}

COMPLEX_ELEMENTAL_ATOM_TYPES = {
    om.Atom.SALT,
    om.Atom.AIR,
    om.Atom.EARTH,
    om.Atom.FIRE,
    om.Atom.WATER,
    om.Atom.VITAE,
    om.Atom.MORS,
    om.Atom.QUINTESSENCE,
}

METAL_ATOM_TYPES = {
    om.Atom.QUICKSILVER,
    om.Atom.LEAD,
    om.Atom.TIN,
    om.Atom.IRON,
    om.Atom.COPPER,
    om.Atom.SILVER,
    om.Atom.GOLD,
}


def find_mixed_no_dupe_puzzles():
    return [
        puzzle
        for puzzle in puzzledb.get_all_puzzles()
        if om.Atom.SALT in puzzle.reagent_types() and
           om.Atom.SALT in puzzle.product_types() and
           any([
               atomtype in ELEMENTAL_ATOM_TYPES
               for atomtype in puzzle.reagent_types()
           ]) and
           any([
               atomtype in ELEMENTAL_ATOM_TYPES
               for atomtype in puzzle.product_types()
           ]) and
           all([
               atomtype in SIMPLE_ELEMENTAL_ATOM_TYPES
               for atomtype in puzzle.atom_types()
           ]) and
           not puzzle.are_parts_available(om.Puzzle.DUPLICATION)
    ]


def main_analyze_types():
    simple_elemental_only_puzzles = 0
    complex_elemental_only_puzzles = 0
    metal_only_puzzles = 0
    for puzzle in puzzledb.get_sample_puzzles():
        product_atoms = [
            atom.type
            for product in puzzle.products
            for atom in product.atoms
        ]
        reagent_atoms = [
            atom.type
            for reagent in puzzle.reagents
            for atom in reagent.atoms
        ]
        if all([
            atom in SIMPLE_ELEMENTAL_ATOM_TYPES
            for atom in (product_atoms + reagent_atoms)
        ]):
            simple_elemental_only_puzzles += 1
        elif all([
            atom in COMPLEX_ELEMENTAL_ATOM_TYPES
            for atom in (product_atoms + reagent_atoms)
        ]):
            complex_elemental_only_puzzles += 1

        if all([
            atom in METAL_ATOM_TYPES
            for atom in (product_atoms + reagent_atoms)
        ]):
            metal_only_puzzles += 1

    print(f"elemental_only_puzzles: {simple_elemental_only_puzzles}")
    print(f"complex_elemental_only_puzzles: {complex_elemental_only_puzzles}")
    print(f"metal_only_puzzles: {metal_only_puzzles}")
