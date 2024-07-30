import om
from pathlib import Path
import os
import puzzleparts

OUTPUT_PATH = Path("output\\24hour-1-test")
SAVEDATA_PATH = Path("savedata")


# === PUZZLE RESOURCES ===

def get_sample_puzzles() -> list[om.Puzzle]:
    return [
        om.Puzzle("resources\\24hour-1-sample\\GEN%03d.puzzle" % i)
        for i in range(100)
    ]


def get_test_puzzles() -> list[om.Puzzle]:
    return [
        om.Puzzle("resources\\24hour-1-test\\GEN%03d.puzzle" % i)
        for i in range(1000)
    ]


def get_base_puzzles() -> list[om.Puzzle]:
    return scan_dir_for_puzzles(Path("resources\\base\\"))


def get_all_puzzles() -> list[om.Puzzle]:
    return get_sample_puzzles() + get_test_puzzles() + get_base_puzzles()


def get_workshop_puzzles() -> list[om.Puzzle]:
    return scan_dir_for_puzzles(Path("savedata\\workshop\\"))


def scan_dir_for_puzzles(path: Path) -> list[om.Puzzle]:
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    puzzle_paths = [
        path / x
        for x in list(os.walk(path))[0][2]
        if x[-7:] == ".puzzle"
    ]
    return [om.Puzzle(str(path)) for path in puzzle_paths]


# === PUZZLE FILTERS ===


def filter_assembly(puzzles: list[om.Puzzle],
        func_any_or_all=all) -> list[om.Puzzle]:
    def output_atom_count(p: om.Puzzle) -> int:
        return sum(len(m.atoms) for m in p.products)

    return sorted(
        [
            p for p in puzzles
            if func_any_or_all(len(mol.atoms) == 1 for mol in p.reagents)
        ],
        key=output_atom_count
    )


def filter_disassembly(puzzles: list[om.Puzzle],
        func_any_or_all=all) -> list[om.Puzzle]:
    def input_atom_count(p: om.Puzzle) -> int:
        return sum(len(m.atoms) for m in p.reagents)

    return sorted(
        [
            p for p in puzzles
            if func_any_or_all(len(mol.atoms) == 1 for mol in p.products)
        ],
        key=input_atom_count
    )


def filter_bond_types_all(puzzles: list[om.Puzzle],
        bond_types: list[int] | set[int]) -> list[om.Puzzle]:
    return filter_bond_types(puzzles, bond_types, all)


def filter_bond_types_any(puzzles: list[om.Puzzle],
        bond_types: list[int] | set[int]) -> list[om.Puzzle]:
    return filter_bond_types(puzzles, bond_types, any)


def filter_bond_types(puzzles: list[om.Puzzle],
        bond_types: list[int] | set[int], func_any_or_all) -> list[om.Puzzle]:
    return [
        p for p in puzzles
        if func_any_or_all(
            b.type in bond_types
                for m in p.reagents + p.products
                for b in m.bonds
        )
    ]


def filter_atom_types(puzzles: list[om.Puzzle], types: list[int] | set[int],
        func_any_or_all=all, negate=False) -> list[om.Puzzle]:
    return [
        p for p in puzzles
        if func_any_or_all(
            [
                type not in types if negate else type in types
                for type in p.atom_types()
            ]
        )
    ]


def filter_uses_parts(puzzles: list[om.Puzzle], parts_list: list[bytes],
        func_any_or_all=all) -> list[om.Puzzle]:
    return [
        p for p in puzzles
        if func_any_or_all(
            part in puzzleparts.full_parts_list(p)
                for part in parts_list
        )
    ]


def filter_bans_parts(puzzles: list[om.Puzzle], parts_list: list[bytes],
        func_any_or_all=all) -> list[om.Puzzle]:
    return [
        p for p in puzzles
        if func_any_or_all(
            part not in puzzleparts.full_parts_list(p)
                for part in parts_list
        )
    ]
