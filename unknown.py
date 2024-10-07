import puzzledb
import step as step


def steps_for_atom(atom_i, steps):
    return [
        s for s in steps
        if str(atom_i) in s or s[0] == "p"
    ]


def available_next_atom_steps(atom_i, previous, steps):
    atom_steps = steps_for_atom(atom_i, steps)
    reagent_step = "r%d" % atom_i
    if reagent_step not in previous:
        return [reagent_step]
    non_prev_atom_steps = [s for s in atom_steps if s not in previous]
    if non_prev_atom_steps:
        return non_prev_atom_steps
    else:
        return "p"


def main():
    puzzle = puzzledb.get_test_puzzles()[391]
    steps = step.find_puzzle_steps(puzzle)
    product = puzzle.products[0]

    for i in range(len(product.atoms)):
        print(i, steps_for_atom(i, steps))

    avail = step.available([], steps)
    print(avail)

    for i in range(len(product.atoms)):
        print(i, steps_for_atom(i, avail))

    for step_1 in avail:
        atom_i = int(step_1[1])
        for step_2 in available_next_atom_steps(atom_i, [step_1], steps):
            print([step_1, step_2])


if __name__ == '__main__':
    main()
