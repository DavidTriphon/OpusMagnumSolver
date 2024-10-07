import puzzledb
import om
import math
import z3


def analyze_step1(puzzle):
    def reagent_counts(solve):
        return [
            solve.get("reagent_%d" % i, 0)
            for i in range(len(puzzle.reagents))
        ]

    def product_counts(solve):
        return [
            solve.get("product_%d" % i, 0)
            for i in range(len(puzzle.products))
        ]

    def operations(solve):
        return sum(v for k, v in solve.items() if "glyph-" != k[:6])

    def solve_parts(solve):
        return [
            part
            for part in om.Part.COSTS.keys()
            if any(part.decode("utf-8") == k[:len(part)] for k in solve.keys())
        ]

    def access_range(solve):
        parts = solve_parts(solve)
        r_accesses = sum(v > 0 for v in reagent_counts(solve))
        p_accesses = sum(v > 0 for v in product_counts(solve))
        part_accesses_min = sum(min(om.Part.ACCESSES[part]) for part in parts)
        part_accesses_max = sum(max(om.Part.ACCESSES[part]) for part in parts)
        access_min = r_accesses + p_accesses + part_accesses_min
        access_max = access_min + part_accesses_max - part_accesses_min
        return range(access_min, access_max + 1)

    def channel_min(solve):
        parts = solve_parts(solve)
        min_arm_length = 1 if all(k in om.Part.COMPATIBLE_MIN_ARM_LENGTH for
            k in parts) else 2
        part_channels = min([
            max([
                    om.Part.channels(part, length)
                    for part in parts
                ] + [1])
            for length in range(min_arm_length, 4)
        ])
        access_channels = math.ceil(min(access_range(solve)) / 6)
        return max(part_channels, access_channels)

    def cost(solve):
        recipe_parts_cost = sum(
            map(lambda k: om.Part.COSTS[k], solve_parts(solve)))
        arm_cost = 20
        channels = channel_min(solve)
        track_cost = (5 * channels if channels > 1 else 0)
        return arm_cost + recipe_parts_cost + track_cost

    solves = solve_step_counts(puzzle)
    solves = sorted(solves,
        key=lambda solve: tuple(
            [cost(solve)] +
            [solve.get(k, 0) for k in ["molecules", "bonds"]] +
            reagent_counts(solve)
        )
    )

    print("Solves (%d):" % (len(solves)))
    for solve in solves:
        access_r = access_range(solve)
        print("\treagents=%s operations=%d cost=%d parts=%s channels=%d "
              "access=%d-%d" % (
            reagent_counts(solve), operations(solve), cost(solve),
            solve_parts(solve), channel_min(solve),
            min(access_r), max(access_r)))
        print("\t" + str(solve))

    part_combos = {tuple(solve_parts(solve)) for solve in solves}
    print("Part Combos (%d):" % (len(part_combos)))
    for combo in part_combos:
        print("\t" + ", ".join(map(lambda b: b.decode("utf-8"), combo)))


def solve_step_counts(puzzle, n=100):
    import z3
    from collections import Counter
    import recipes

    def atom_type_key(type):
        return "atom_%s" % om.Atom.TYPE_NAMES[type]

    bonds = z3.Int("bonds")
    atoms = [z3.Int(atom_type_key(t)) for t in om.Atom.TYPES]
    molecules = z3.Int("molecules")
    count_vars = [bonds, molecules, *atoms]

    ###
    bond_actions = {
        om.Part.BONDER.decode("utf-8"): {"bonds": 1, "molecules": -1},
        om.Part.UNBONDER.decode("utf-8"): {"bonds": -1, "molecules": 1},
    }

    ###
    reagent_actions = {
        ("reagent_%d" % i): (
                {
                    "bonds": len(reagent.bonds),
                    "molecules": 1,
                } |
                Counter(map(
                    lambda a: atom_type_key(a.type),
                    reagent.atoms
                ))
        )
        for i, reagent in enumerate(puzzle.reagents)
    }

    ###
    product_actions = {
        ("product_%d" % i): (
                {
                    "bonds": len(product.bonds),
                    "molecules": 1,
                } |
                Counter(map(
                    lambda a: atom_type_key(a.type),
                    product.atoms
                ))
        )
        for i, product in enumerate(puzzle.products)
    }
    for p_action in product_actions.values():
        for k in p_action.keys():
            p_action[k] = -p_action[k]

    ###
    def recipe_name(recipe):
        name = recipe.part.decode("utf-8")
        if recipe.part == om.Part.CALCIFICATION:
            return name + "_" + om.Atom.TYPE_NAMES[recipe.conversions[0][0]]
        if recipe.part == om.Part.DUPLICATION:
            return name + "_" + om.Atom.TYPE_NAMES[recipe.conversions[0][1]]
        if recipe.part == om.Part.PURIFICATION:
            return name + "_" + om.Atom.TYPE_NAMES[recipe.produced[0]]
        if recipe.part == om.Part.PROJECTION:
            return name + "_" + om.Atom.TYPE_NAMES[recipe.conversions[0][1]]
        return name

    def recipe_offsets(recipe):
        offsets = Counter(map(atom_type_key, recipe.list_created()))
        offsets.subtract(Counter(map(atom_type_key, recipe.list_required())))
        offsets["molecules"] = offsets.total()
        return offsets

    recipe_actions = {
        recipe_name(recipe): recipe_offsets(recipe)
        for recipe in recipes.RECIPES
    }

    ###
    actions = bond_actions | recipe_actions | reagent_actions | product_actions

    ###
    action_vars = {action_name: z3.Int(action_name) for action_name in
        actions.keys()}

    ###
    all_vars = count_vars + list(action_vars.values())
    nonneg_constraints = [var >= 0 for var in all_vars]

    ###
    action_constraints = [
        var == sum(
            action.get(str(var), 0) * action_vars[name]
                for name, action in actions.items()
        )
        for var in count_vars
    ]

    ###
    puzzle_parts_list = puzzle.full_parts_list()
    availability_constraints = [
        z3.Or((recipe.part in puzzle_parts_list),
            action_vars[recipe_name(recipe)] == 0)
        for recipe in recipes.RECIPES
    ]

    ###
    product_constraints = [
        action_vars[product_key_name] > 0
        for product_key_name in product_actions.keys()
    ]

    ###
    solver = z3.Solver()

    solver.add(nonneg_constraints + action_constraints +
               availability_constraints + product_constraints)
    solves = find_n_solves(all_vars, solver, n)

    return [
        {
            str(variable): value
            for variable, value in zip(all_vars, solve)
            if not isinstance(value, int) or value != 0
        }
        for solve in solves
    ]


def find_n_solves(variables, solver: z3.Solver, n=100):
    results_table = []
    while solver.check() == z3.sat and len(results_table) < n:
        m = solver.model()
        results = tuple(m.evaluate(var).as_long() for var in variables)
        results_table.append(results)
        solver.add([z3.Or([
            val < result
            for val, result in zip(variables, results)
        ])])

    # enforce pareto, in case first solves weren't pareto
    # there is no guarantee in the z3 solver that all solutions are pareto,
    # just that they are pareto among all previous solutions.
    reversed_results = results_table[::-1]
    for results in reversed_results:
        if results in results_table:
            others = results_table[:results_table.index(results)]
            for other in others:
                if all(other[i] >= results[i] for i in range(len(results))):
                    results_table.remove(other)

    return results_table


def main_part_solve():
    from timeit import default_timer as timer

    puzzle = puzzledb.get_test_puzzles()[128]
    start = timer()
    analyze_step1(puzzle)
    end = timer()
    print("duration: %.2f secs" % (end - start))


if __name__ == '__main__':
    main_part_solve()
