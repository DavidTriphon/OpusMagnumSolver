import om
import puzzledb


def step_preconditions(step):
    if step[0] == "c":
        return []
    if step[0] == "b":
        return []


def available(chosen, steps):
    return [
        s for s in steps
        if s not in chosen and all(x in chosen for x in step_preconditions(s))
    ]


def generator_permutations(values, prefix=None):
    if prefix is None:
        prefix = []
    if len(values) == 0:
        yield tuple(prefix)
    for val in values:
        for x in generator_permutations(values - {val}, prefix + [val]):
            yield x


def step_order_generator(steps, prefix=None):
    if prefix is None:
        prefix = []
    if len(prefix) == len(steps):
        yield prefix
    for a in available(prefix, steps):
        for p in step_order_generator(steps, prefix + [a]):
            yield p


def moves_from_steps(all_steps):
    return [(int(i), s) for s in all_steps for i in s[1:]]


def find_puzzle_steps(puzzle: om.Puzzle):
    # assumes 1 reagent with 1 atom and 1 product for which every atom type
    # has only one manufacture path
    product = puzzle.products[0]
    reagent = puzzle.reagents[0]
    pr_mapping = {i: 0 for i in range(len(puzzle.products[0].atoms))}
    return [
        # assume all bonds come from bonder
        ("b%d%d" % (
            next((i for i, a in enumerate(product.atoms)
                if bond.positions[0] == a.position), None),
            next((i for i, a in enumerate(product.atoms)
                if bond.positions[1] == a.position), None)
        )) for bond in product.bonds
    ] + [
        "c%d" % p_i
        for p_i, p_atom in enumerate(puzzle.products[0].atoms)
        # assumes the only path is identity or calcification
        if p_atom.type != reagent.atoms[pr_mapping[p_i]].type
    ]


def main():
    puzzle = puzzledb.get_test_puzzles()[391]
    paths = step_order_generator(find_puzzle_steps(puzzle))
    for i in range(100):
        print(i, next(paths))


if __name__ == '__main__':
    main()
