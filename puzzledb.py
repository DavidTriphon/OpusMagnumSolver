import om
from pathlib import Path
import os
import puzzleparts


OUTPUT_PATH = Path("output\\24hour-1-test")
SAVEDATA_PATH = Path("savedata")

# === PUZZLE RESOURCES ===

def get_sample_puzzles():
    return [
        om.Puzzle("resources\\24hour-1-sample\\GEN%03d.puzzle" % i)
        for i in range(100)
    ]


def get_test_puzzles():
    return [
        om.Puzzle("resources\\24hour-1-test\\GEN%03d.puzzle" % i)
        for i in range(1000)
    ]


def get_base_puzzles():
    return scan_dir_for_puzzles(Path("resources\\base\\"))


def get_all_puzzles():
    return get_sample_puzzles() + get_test_puzzles() + get_base_puzzles()


def get_workshop_puzzles():
    return scan_dir_for_puzzles(Path("savedata\\workshop\\"))


def scan_dir_for_puzzles(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    puzzle_paths = [
        path / x
        for x in list(os.walk(path))[0][2]
        if x[-7:] == ".puzzle"
    ]
    return [om.Puzzle(str(path)) for path in puzzle_paths]


# === PUZZLE FILTERS ===


def filter_assembly(puzzles, func_any_or_all=all):
    return sorted(
        [
            p for p in puzzles
            if func_any_or_all(len(mol.atoms) == 1 for mol in p.reagents)
        ],
        key=lambda p: sum(len(m.atoms) for m in p.products)
    )


def filter_disassembly(puzzles, func_any_or_all=all):
    return sorted(
        [
            p for p in puzzles
            if func_any_or_all(len(mol.atoms) == 1 for mol in p.products)
        ],
        key=lambda p: sum(len(m.atoms) for m in p.reagents)
    )


def filter_bond_types_all(puzzles, bond_types):
    return filter_bond_types(puzzles, bond_types, all)


def filter_bond_types_any(puzzles, bond_types):
    return filter_bond_types(puzzles, bond_types, any)


def filter_bond_types(puzzles, bond_types, func_any_or_all):
    return [
        p for p in puzzles
        if func_any_or_all(
            b.type in bond_types
                for m in p.reagents + p.products
                for b in m.bonds
        )
    ]


def filter_atom_types(puzzles, types, func_any_or_all=all):
    return [
        p for p in puzzles
        if func_any_or_all(
            [
                type in types
                for type in p.atom_types()
            ]
        )
    ]


def filter_uses_parts(puzzles, parts_list, func_any_or_all=all):
    return [
        p for p in puzzles
        if func_any_or_all(
            part in puzzleparts.full_parts_list(p)
                for part in parts_list
        )
    ]


def filter_bans_parts(puzzles, parts_list, func_any_or_all=all):
    return [
        p for p in puzzles
        if func_any_or_all(
            part not in puzzleparts.full_parts_list(p)
                for part in parts_list
        )
    ]
