import om

AVAILABLE_FLAG_NAMES = {
    0: "ARM",
    1: "MULTIARM",
    2: "PISTON",
    3: "TRACK",
    8: "BONDER",
    9: "UNBONDER",
    10: "MULTIBONDER",
    11: "TRIPLEX",
    12: "CALCIFICATION",
    13: "DUPLICATION",
    14: "PROJECTION",
    15: "PURIFICATION",
    16: "ANIMISMUS",
    17: "DISPOSAL",
    18: "QUINTESSENCE",
    22: "GRAB_AND_ROTATE",
    23: "DROP",
    24: "RESET",
    25: "REPEAT",
    26: "PIVOT",
    28: "BERLO",
}

PARTS_PER_FLAG = {
    0: [om.Part.ARM1],
    1: [om.Part.ARM2, om.Part.ARM3, om.Part.ARM6],
    2: [om.Part.PISTON],
    3: [om.Part.TRACK],
    8: [om.Part.BONDER],
    9: [om.Part.UNBONDER],
    10: [om.Part.MULTIBONDER],
    11: [om.Part.TRIPLEX],
    12: [om.Part.CALCIFICATION],
    13: [om.Part.DUPLICATION],
    14: [om.Part.PROJECTION],
    15: [om.Part.PURIFICATION],
    16: [om.Part.ANIMISMUS],
    17: [om.Part.DISPOSAL],
    18: [om.Part.UNIFICATION, om.Part.DISPERSION],
    28: [om.Part.BERLO],
}


def flag_names(puzzle: om.Puzzle):
    return [
        AVAILABLE_FLAG_NAMES[i]
        if i in AVAILABLE_FLAG_NAMES
        else str(i)
        for i in all_puzzle_flags(puzzle)
    ]


def all_puzzle_flags(puzzle: om.Puzzle):
    return [i for i in range(32) if 1 << i & puzzle.parts_available]


def nondef_puzzle_flags(puzzle: om.Puzzle):
    return [i for i in range(32)
        if 1 << i & puzzle.parts_available & ~om.Puzzle.DEFAULT_PARTS_AVAILABLE]


def full_parts_list(puzzle: om.Puzzle) -> list[bytes]:
    return [
        part
        for flag in all_puzzle_flags(puzzle)
        if flag in PARTS_PER_FLAG
        for part in PARTS_PER_FLAG[flag]
    ]


def nondef_parts_list(puzzle: om.Puzzle) -> list[bytes]:
    return [
        part
        for flag in nondef_puzzle_flags(puzzle)
        if flag in PARTS_PER_FLAG
        for part in PARTS_PER_FLAG[flag]
    ]