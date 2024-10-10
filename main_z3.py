import puzzledb
import om
import math
import z3
from collections import Counter
import recipes


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

    puzzle = puzzledb.get_test_puzzles()[391]
    start = timer()
    solution = solve(puzzle)
    end = timer()
    print()
    print("duration: %.2f secs" % (end - start))
    if solution is not None:
        print("solution found")
    else:
        print("no solution")


def bounds(lower_bound, expr, upper_bound):
    return z3.And([
        lower_bound <= expr,
        expr < upper_bound
    ])


def solve(puzzle) -> om.Solution | None:
    solver = z3.Solver()

    puzzle_recipes = recipes.all_puzzle_recipes(puzzle)

    recipe_strs = {recipe: recipe.shorthand_str() for recipe in puzzle_recipes}
    atom_type_strs = {t: om.Atom.TYPE_NAMES[t] for t in om.Atom.TYPES}
    bond_type_strs = {t: om.Bond.type_str_id(t) for t in om.Bond.TYPES}

    # # types

    RecipeSort, RecipeSort_consts = z3.EnumSort("Recipe",
        list(recipe_strs.values()))
    AtomTypeSort, AtomTypeSort_consts = z3.EnumSort("AtomType",
        list(atom_type_strs.values()))
    BondTypeSort, BondTypeSort_consts = z3.EnumSort("BondType",
        list(bond_type_strs.values()))

    recipe_const_dict = {
        recipe: recipe_const
        for recipe, recipe_const in zip(puzzle_recipes, RecipeSort_consts)
    }
    atomtype_const_dict = {
        atomtype: atomtype_const
        for atomtype, atomtype_const in zip(om.Atom.TYPES, AtomTypeSort_consts)
    }
    bondtype_const_dict = {
        bondtype: bondtype_const
        for bondtype, bondtype_const in zip(om.Bond.TYPES, BondTypeSort_consts)
    }

    # # relationships and fields

    # structure =
    # array of atoms
    # array of bonds
    # dict from recipetype to array of recipes

    # bond/atom = {
    #   consumer: recipe reference
    #   producer: recipe reference
    #   array of conversions
    # }
    # conversion = recipe reference

    # recipe reference = {
    #   recipe
    #   recipe index
    #   material index
    # }

    # recipe = {
    #   array of atom/bond consumed
    #   array of atom/bond produced
    #   array of atom/bond converted
    #   array of atom/bond hovers
    # }


    # basic single ints

    f_recipe_count = z3.Function("recipe_count", RecipeSort, z3.IntSort())
    atom_count = z3.Int("atom_mass_count")
    bond_count = z3.Int("bond_mass_count")

    forall_recipe = z3.Const("forall_recipe", RecipeSort)
    solver.add([
        # bounded values for recipe_count
        z3.ForAll([forall_recipe], 0 <= f_recipe_count(forall_recipe)),
        0 <= atom_count,
        0 <= bond_count,
    ])

    # fields of recipes: list of atoms and bond consumed/produced/converted

    f_recipe_atom_consumed_count = z3.Function("recipe_atom_consumed_count",
        RecipeSort, z3.IntSort())
    f_recipe_atom_produced_count = z3.Function("recipe_atom_produced_count",
        RecipeSort, z3.IntSort())
    f_recipe_atom_conversion_count = z3.Function("recipe_atom_conversion_count",
        RecipeSort, z3.IntSort())
    f_recipe_bond_consumed_count = z3.Function("recipe_bond_consumed_count",
        RecipeSort, z3.IntSort())
    f_recipe_bond_produced_count = z3.Function("recipe_bond_produced_count",
        RecipeSort, z3.IntSort())
    f_recipe_bond_conversion_count = z3.Function("recipe_bond_conversion_count",
        RecipeSort, z3.IntSort())

    # bind funcs from above to data in the recipes
    solver.add([
        z3.And([
            f_recipe_atom_consumed_count(recipe_const_dict[recipe]) == len(
                recipe.consumed),
            f_recipe_atom_produced_count(recipe_const_dict[recipe]) == len(
                recipe.produced),
            f_recipe_atom_conversion_count(recipe_const_dict[recipe]) == len(
                recipe.conversions),
            f_recipe_bond_consumed_count(recipe_const_dict[recipe]) == len(
                recipe.bonds_consumed),
            f_recipe_bond_produced_count(recipe_const_dict[recipe]) == len(
                recipe.bonds_produced),
            f_recipe_bond_conversion_count(recipe_const_dict[recipe]) == len(
                recipe.bond_conversions),
        ])
        for recipe in puzzle_recipes
    ])

    # values in the recipe field lists: each atom and bond type for
    # consumptions/productions/conversions

    # the type of a material for a given recipe, field, and index
    # each func is a different field, so the signature is recipe and mat_index
    f_recipe_index_to_atom_consumed_type = z3.Function(
        "recipe_index_to_atom_consumed_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )
    f_recipe_index_to_atom_produced_type = z3.Function(
        "recipe_index_to_atom_produced_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )
    f_recipe_index_to_atom_conversion_in_type = z3.Function(
        "recipe_index_to_atom_conversion_in_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )
    f_recipe_index_to_atom_conversion_out_type = z3.Function(
        "recipe_index_to_atom_conversion_out_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )
    f_recipe_index_to_bond_consumed_type = z3.Function(
        "recipe_index_to_bond_consumed_type",
        RecipeSort, z3.IntSort(), BondTypeSort
    )
    f_recipe_index_to_bond_produced_type = z3.Function(
        "recipe_index_to_bond_produced_type",
        RecipeSort, z3.IntSort(), BondTypeSort
    )
    f_recipe_index_to_bond_conversion_in_type = z3.Function(
        "recipe_index_to_bond_conversion_in_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )
    f_recipe_index_to_bond_conversion_out_type = z3.Function(
        "recipe_index_to_bond_conversion_out_type",
        RecipeSort, z3.IntSort(), AtomTypeSort
    )

    def material_type_constraint(recipe_const,
            f_recipe_index_field_to_type, recipe_material_types):
        return z3.And([
            f_recipe_index_field_to_type(
                recipe_const, material_index) == material_type
            for material_index, material_type in
            enumerate(recipe_material_types)
        ])

    # bind the above funcs to the data in the recipes
    solver.add([
        z3.And([
            # atom consumed/produced types
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_atom_consumed_type,
                [atomtype_const_dict[t] for t in recipe.consumed]
            ),
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_atom_produced_type,
                [atomtype_const_dict[t] for t in recipe.produced]
            ),

            # atom conversion types
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_atom_conversion_in_type,
                [atomtype_const_dict[t] for (t, _) in recipe.conversions]
            ),
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_atom_conversion_out_type,
                [atomtype_const_dict[t] for (_, t) in recipe.conversions]
            ),

            # bond consumed/produced types
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_bond_consumed_type,
                [bondtype_const_dict[t] for t in recipe.bonds_consumed]
            ),
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_bond_produced_type,
                [bondtype_const_dict[t] for t in recipe.bonds_produced]
            ),

            # bond conversion types
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_bond_conversion_in_type,
                [bondtype_const_dict[t] for (t, _) in recipe.bond_conversions]
            ),
            material_type_constraint(
                recipe_const_dict[recipe],
                f_recipe_index_to_bond_conversion_out_type,
                [bondtype_const_dict[t] for (_, t) in recipe.bond_conversions]
            ),

        ])
        for recipe in puzzle_recipes
    ])

    # every application of a recipe acts on real material, that material has
    # an id.

    # atom ids
    f_recipe_id_atom_consumed_index_to_atomid = z3.Function(
        "recipe_id_atom_consumed_index_to_atomid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_recipe_id_atom_produced_index_to_atomid = z3.Function(
        "recipe_id_atom_produced_index_to_atomid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_recipe_id_atom_converted_index_to_atomid = z3.Function(
        "recipe_id_atom_converted_index_to_atomid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    # bond ids
    f_recipe_id_bond_consumed_index_to_bondid = z3.Function(
        "recipe_id_bond_consumed_index_to_bondid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_recipe_id_bond_produced_index_to_bondid = z3.Function(
        "recipe_id_bond_produced_index_to_bondid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_recipe_id_bond_converted_index_to_bondid = z3.Function(
        "recipe_id_bond_converted_index_to_bondid",
        RecipeSort, z3.IntSort(), z3.IntSort(), z3.IntSort()
    )

    forall_recipeid = z3.Int("forall_recipeid")
    forall_materialid = z3.Int("forall_materialid")

    # material count implies the existence of contiguous mass id assigns
    # the ids they return must be within bounds
    def material_id_constraint(f_field_count, f_field_id, material_count):
        return z3.ForAll([forall_recipe, forall_recipeid, forall_materialid],
            z3.Implies(
                z3.And([
                    bounds(0, forall_recipeid, f_recipe_count(forall_recipe)),
                    bounds(0, forall_materialid, f_field_count(forall_recipe)),
                ]),
                bounds(
                    0,
                    f_field_id(
                        forall_recipe, forall_recipeid, forall_materialid
                    ),
                    material_count)
            )
        )

    solver.add([
        # atoms consumed/produced/conversions
        material_id_constraint(
            f_recipe_atom_consumed_count,
            f_recipe_id_atom_consumed_index_to_atomid,
            atom_count
        ),
        material_id_constraint(
            f_recipe_atom_produced_count,
            f_recipe_id_atom_produced_index_to_atomid,
            atom_count
        ),
        material_id_constraint(
            f_recipe_atom_conversion_count,
            f_recipe_id_atom_converted_index_to_atomid,
            atom_count
        ),
        # bonds consumed/produced
        material_id_constraint(
            f_recipe_bond_consumed_count,
            f_recipe_id_bond_consumed_index_to_bondid,
            bond_count
        ),
        material_id_constraint(
            f_recipe_bond_produced_count,
            f_recipe_id_bond_produced_index_to_bondid,
            bond_count
        ),
        material_id_constraint(
            f_recipe_bond_conversion_count,
            f_recipe_id_bond_converted_index_to_bondid,
            bond_count
        ),
    ])

    # reverse relationship from atomids to the recipe, id, mat_index
    # consumer
    f_atomid_to_consumer_recipe = z3.Function("f_atomid_to_consumer_recipe",
        z3.IntSort(), RecipeSort)
    f_atomid_to_consumer_recipeid = z3.Function("f_atomid_to_consumer_recipeid",
        z3.IntSort(), z3.IntSort())
    f_atomid_to_consumer_matindex = z3.Function("f_atomid_to_consumer_matindex",
        z3.IntSort(), z3.IntSort())
    # producer
    f_atomid_to_producer_recipe = z3.Function("f_atomid_to_producer_recipe",
        z3.IntSort(), RecipeSort)
    f_atomid_to_producer_recipeid = z3.Function("f_atomid_to_producer_recipeid",
        z3.IntSort(), z3.IntSort())
    f_atomid_to_producer_matindex = z3.Function("f_atomid_to_producer_matindex",
        z3.IntSort(), z3.IntSort())

    # conversions
    f_atomid_to_conversion_count = z3.Function(
        "atomid_to_conversion_count",
        z3.IntSort(), z3.IntSort()
    )

    f_atomid_conversion_order_to_conversion_recipe = z3.Function(
        "atomid_conversion_order_to_conversion_recipe",
        z3.IntSort(), z3.IntSort(), RecipeSort
    )
    f_atomid_conversion_order_to_conversion_recipeid = z3.Function(
        "atomid_conversion_order_to_conversion_recipeid",
        z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_atomid_conversion_order_to_conversion_matindex = z3.Function(
        "atomid_conversion_order_to_conversion_matindex",
        z3.IntSort(), z3.IntSort(), z3.IntSort()
    )

    forall_atomid = z3.Int("forall_atomid")
    forall_conv_order = z3.Int("forall_conv_order")
    solver.add([
        z3.ForAll([forall_atomid],
            z3.Implies(
                bounds(0, forall_atomid, atom_count),
                z3.And([
                    # results of return functions are valid
                    bounds(0,
                        f_atomid_to_consumer_recipeid(forall_atomid),
                        f_recipe_count(
                            f_atomid_to_consumer_recipe(forall_atomid))
                    ),
                    bounds(0,
                        f_atomid_to_producer_recipeid(forall_atomid),
                        f_recipe_count(
                            f_atomid_to_producer_recipe(forall_atomid))
                    ),
                    bounds(0,
                        f_atomid_to_consumer_matindex(forall_atomid),
                        f_recipe_atom_consumed_count(
                            f_atomid_to_consumer_recipe(forall_atomid))
                    ),
                    bounds(0,
                        f_atomid_to_producer_matindex(forall_atomid),
                        f_recipe_atom_produced_count(
                            f_atomid_to_producer_recipe(forall_atomid))
                    ),
                    # consumed is inverse function
                    f_recipe_id_atom_consumed_index_to_atomid(
                        f_atomid_to_consumer_recipe(forall_atomid),
                        f_atomid_to_consumer_recipeid(forall_atomid),
                        f_atomid_to_consumer_matindex(forall_atomid)
                    ) == forall_atomid,
                    # produced is inverse function
                    f_recipe_id_atom_produced_index_to_atomid(
                        f_atomid_to_producer_recipe(forall_atomid),
                        f_atomid_to_producer_recipeid(forall_atomid),
                        f_atomid_to_producer_matindex(forall_atomid)
                    ) == forall_atomid,
                    # conversion count is positive
                    0 <= f_atomid_to_conversion_count(forall_atomid),
                ])

            )
        ),

        z3.ForAll([forall_atomid, forall_conv_order],
            z3.Implies(
                z3.And([
                    bounds(0, forall_atomid, atom_count),
                    bounds(0, forall_conv_order,
                        f_atomid_to_conversion_count(forall_atomid)),
                ]),
                f_recipe_id_atom_converted_index_to_atomid(
                    f_atomid_conversion_order_to_conversion_recipe(
                        forall_atomid, forall_conv_order),
                    f_atomid_conversion_order_to_conversion_recipeid(
                        forall_atomid, forall_conv_order),
                    f_atomid_conversion_order_to_conversion_matindex(
                        forall_atomid, forall_conv_order)
                ) == forall_atomid
            )
        )
    ])

    if False:
        forall_recipe = z3.Const("forall_recipe", RecipeSort)
        forall_recipeid = z3.Int("forall_recipeid")
        forall_matindex = z3.Int("forall_matindex")

        solver.add(
            z3.ForAll(
                [forall_atomid, forall_recipe, forall_recipeid,
                    forall_matindex],
                z3.Implies(
                    z3.And([
                        bounds(0, forall_atomid, atom_count),
                        bounds(0, forall_recipeid,
                            f_recipe_count(forall_recipe)),
                    ]),
                    z3.And([
                        z3.Implies(
                            bounds(0, forall_matindex,
                                f_recipe_atom_consumed_count(forall_recipe)),
                            # consumed is inverse function
                            (f_recipe_id_atom_consumed_index_to_atomid(
                                forall_recipe, forall_recipeid,
                                forall_matindex
                            ) == forall_atomid) == z3.And([
                                f_atomid_to_consumer_recipe(
                                    forall_atomid) == forall_recipe,
                                f_atomid_to_consumer_recipeid(
                                    forall_atomid) == forall_recipeid,
                                f_atomid_to_consumer_matindex(
                                    forall_atomid) == forall_matindex
                            ])
                        ),
                        z3.Implies(
                            bounds(0, forall_matindex,
                                f_recipe_atom_produced_count(forall_recipe)),
                            # produced is inverse function
                            (f_recipe_id_atom_produced_index_to_atomid(
                                forall_recipe, forall_recipeid,
                                forall_matindex
                            ) == forall_atomid) == z3.And([
                                f_atomid_to_producer_recipe(
                                    forall_atomid) == forall_recipe,
                                f_atomid_to_producer_recipeid(
                                    forall_atomid) == forall_recipeid,
                                f_atomid_to_producer_matindex(
                                    forall_atomid) == forall_matindex
                            ])
                        )
                    ])
                )
            ),

        )

    # reverse relationship from bondids to the recipe, id, mat_index
    # consumer
    f_bondid_to_consumer_recipe = z3.Function("f_bondid_to_consumer_recipe",
        z3.IntSort(), RecipeSort)
    f_bondid_to_consumer_recipeid = z3.Function("f_bondid_to_consumer_recipeid",
        z3.IntSort(), z3.IntSort())
    f_bondid_to_consumer_matindex = z3.Function("f_bondid_to_consumer_matindex",
        z3.IntSort(), z3.IntSort())
    # producer
    f_bondid_to_producer_recipe = z3.Function("f_bondid_to_producer_recipe",
        z3.IntSort(), RecipeSort)
    f_bondid_to_producer_recipeid = z3.Function("f_bondid_to_producer_recipeid",
        z3.IntSort(), z3.IntSort())
    f_bondid_to_producer_matindex = z3.Function("f_bondid_to_producer_matindex",
        z3.IntSort(), z3.IntSort())

    # conversions
    f_bondid_to_conversion_count = z3.Function(
        "bondid_to_conversion_count",
        z3.IntSort(), z3.IntSort()
    )

    f_bondid_conversion_order_to_conversion_recipe = z3.Function(
        "bondid_conversion_order_to_conversion_recipe",
        z3.IntSort(), z3.IntSort(), RecipeSort
    )
    f_bondid_conversion_order_to_conversion_recipeid = z3.Function(
        "bondid_conversion_order_to_conversion_recipeid",
        z3.IntSort(), z3.IntSort(), z3.IntSort()
    )
    f_bondid_conversion_order_to_conversion_matindex = z3.Function(
        "bondid_conversion_order_to_conversion_matindex",
        z3.IntSort(), z3.IntSort(), z3.IntSort()
    )

    forall_bondid = z3.Int("forall_bondid")
    forall_conv_order = z3.Int("forall_conv_order")
    solver.add([
        z3.ForAll([forall_bondid],
            z3.Implies(
                bounds(0, forall_bondid, bond_count),
                z3.And([
                    # results of return functions are valid
                    bounds(0,
                        f_bondid_to_consumer_recipeid(forall_bondid),
                        f_recipe_count(
                            f_bondid_to_consumer_recipe(forall_bondid))
                    ),
                    bounds(0,
                        f_bondid_to_producer_recipeid(forall_bondid),
                        f_recipe_count(
                            f_bondid_to_producer_recipe(forall_bondid))
                    ),
                    bounds(0,
                        f_bondid_to_consumer_matindex(forall_bondid),
                        f_recipe_bond_consumed_count(
                            f_bondid_to_consumer_recipe(forall_bondid))
                    ),
                    bounds(0,
                        f_bondid_to_producer_matindex(forall_bondid),
                        f_recipe_bond_produced_count(
                            f_bondid_to_producer_recipe(forall_bondid))
                    ),
                    # consumed is inverse function
                    f_recipe_id_bond_consumed_index_to_bondid(
                        f_bondid_to_consumer_recipe(forall_bondid),
                        f_bondid_to_consumer_recipeid(forall_bondid),
                        f_bondid_to_consumer_matindex(forall_bondid)
                    ) == forall_bondid,
                    # produced is inverse function
                    f_recipe_id_bond_produced_index_to_bondid(
                        f_bondid_to_producer_recipe(forall_bondid),
                        f_bondid_to_producer_recipeid(forall_bondid),
                        f_bondid_to_producer_matindex(forall_bondid)
                    ) == forall_bondid,
                    # conversion count is positive
                    0 <= f_bondid_to_conversion_count(forall_bondid),
                ])
            )
        ),
        z3.ForAll([forall_bondid, forall_conv_order],
            z3.Implies(
                z3.And([
                    bounds(0, forall_bondid, bond_count),
                    bounds(0, forall_conv_order,
                        f_bondid_to_conversion_count(forall_bondid)),
                ]),
                f_recipe_id_bond_converted_index_to_bondid(
                    f_bondid_conversion_order_to_conversion_recipe(
                        forall_bondid, forall_conv_order),
                    f_bondid_conversion_order_to_conversion_recipeid(
                        forall_bondid, forall_conv_order),
                    f_bondid_conversion_order_to_conversion_matindex(
                        forall_bondid, forall_conv_order)
                ) == forall_bondid
            )
        )
    ])

    # TODO: add functions defining relationships between data

    input_recipe = next(r for r in puzzle_recipes
        if isinstance(r.part, tuple) and r.part[0] == om.Part.INPUT)
    output_recipe = next(r for r in puzzle_recipes
        if isinstance(r.part, tuple) and r.part[0] == om.Part.OUTPUT_STANDARD)
    calc_recipe = next(r for r in puzzle_recipes
        if r.part == om.Part.CALCIFICATION)
    bond_recipe = next(r for r in puzzle_recipes
        if r.part == om.Part.BONDER)

    #solver.add(f_recipe_count(recipe_const_dict[input_recipe]) == 4)
    #solver.add(f_recipe_count(recipe_const_dict[output_recipe]) == 1)
    #solver.add(f_recipe_count(recipe_const_dict[calc_recipe]) == 3)
    #solver.add(f_recipe_count(recipe_const_dict[bond_recipe]) == 3)

    # the outputs must all be more than 0
    for recipe in puzzle_recipes:
        # TODO: fix this stopgap
        solver.add(f_recipe_count(recipe_const_dict[recipe]) < 5)
        if isinstance(recipe.part, tuple):
            if recipe.part[0] == om.Part.OUTPUT_STANDARD:
                solver.add(f_recipe_count(recipe_const_dict[recipe]) == 1)

    #solver.add(atom_count == 4)
    #solver.add(bond_count == 3)

    # value to optimize for
    """
    solver.minimize(
        sum(
            f_recipe_count(recipe_const_dict[recipe])
                for recipe in puzzle_recipes
        )
    )
    """

    print("thinking")
    if not solver.check() == z3.sat:
        print("gave up")
        return None
    print("done thinking")

    model = solver.model()
    solution = om.Solution()

    recipe_counts = {
        recipe: model.evaluate(f_recipe_count(recipe_const_dict[
            recipe])).as_long()
        for recipe in puzzle_recipes
    }
    atoms = model.evaluate(atom_count).as_long()
    bonds = model.evaluate(bond_count).as_long()

    print("atom count: %r" % atoms)
    print("bond count: %r" % bonds)
    print("recipe counts:")
    for recipe in puzzle_recipes:
        print("\t%s: %r" % (recipe.shorthand_str(), recipe_counts[recipe]))

    def print_recipe_ids(recipe, id, mat_list_names, func_for_id,
            prefix_str):
        recipe_const = recipe_const_dict[recipe]
        for matindex, mattype in enumerate(mat_list_names):
            atomid = model.evaluate(
                func_for_id(
                    recipe_const, id, matindex
                )
            ).as_long()
            print("\t\t%s %s: %d" % (prefix_str, mattype, atomid))

    print()
    print("recipe material ids")
    for recipe in puzzle_recipes:
        print(recipe.shorthand_str())
        for id in range(recipe_counts[recipe]):
            print("\t[%d]" % id)
            print_recipe_ids(recipe, id, [
                om.Atom.TYPE_NAMES[t] for t in recipe.consumed
            ],
                f_recipe_id_atom_consumed_index_to_atomid,
                "atom consumed")
            print_recipe_ids(recipe, id, [
                om.Atom.TYPE_NAMES[t] for t in recipe.produced
            ],
                f_recipe_id_atom_produced_index_to_atomid,
                "atom produced")
            print_recipe_ids(recipe, id, [
                om.Atom.TYPE_NAMES[t1] + "_to_" + om.Atom.TYPE_NAMES[t2]
                for (t1, t2) in recipe.conversions
            ],
                f_recipe_id_atom_converted_index_to_atomid,
                "atom converted")
            print_recipe_ids(recipe, id, [
                om.Bond.type_str_id(t) for t in recipe.bonds_consumed
            ],
                f_recipe_id_bond_consumed_index_to_bondid,
                "bond consumed")
            print_recipe_ids(recipe, id, [
                om.Bond.type_str_id(t) for t in recipe.bonds_produced
            ],
                f_recipe_id_bond_produced_index_to_bondid,
                "bond produced")
            print_recipe_ids(recipe, id, [
                om.Bond.type_str_id(t1) + "_to_" + om.Bond.type_str_id(t2)
                for (t1, t2) in recipe.bond_conversions
            ],
                f_recipe_id_bond_converted_index_to_bondid,
                "bond converted")

    print()
    print("atoms:")
    for i in range(atoms):
        print("\t[%d]" % i)
        producer_recipe = model.evaluate(f_atomid_to_producer_recipe(i))
        producer_recipeid = model.evaluate(f_atomid_to_producer_recipeid(i))
        producer_matindex = model.evaluate(f_atomid_to_producer_matindex(i))
        print("\t\tproduced by %r, id=%r, mat_index=%r" % (
            producer_recipe, producer_recipeid, producer_matindex
        ))
        consumer_recipe = model.evaluate(f_atomid_to_consumer_recipe(i))
        consumer_recipeid = model.evaluate(f_atomid_to_consumer_recipeid(i))
        consumer_matindex = model.evaluate(f_atomid_to_consumer_matindex(i))
        print("\t\tconsumed by %r, id=%r, mat_index=%r" % (
            consumer_recipe, consumer_recipeid, consumer_matindex
        ))
        conversions = model.evaluate(f_atomid_to_conversion_count(i))
        print("\t\tconversions = %r" % conversions)

    return solution


def test():
    solver = z3.Solver()
    left_count = z3.Int("left_count")
    right_count = z3.Int("right_count")
    right_int = z3.Int("right_int")
    left_int = z3.Int("left_int")
    left = z3.Function("left", z3.IntSort(), z3.IntSort())
    right = z3.Function("right", z3.IntSort(), z3.IntSort())

    rules = [
        z3.ForAll([right_int, left_int],
            z3.Implies(
                z3.And(0 < right_int, right_int <= right_count,
                    0 < left_int, left_int <= left_count),
                (right(left_int) == right_int) == (left(right_int) == left_int)
            )
        ),
        z3.ForAll(
            [right_int],
            z3.Implies(
                z3.And(0 < right_int, right_int <= right_count),
                z3.And(0 < left(right_int), left(right_int) <= left_count)
            )
        ),
        z3.ForAll(
            [left_int],
            z3.Implies(
                z3.And(0 < left_int, left_int <= left_count),
                z3.And(0 < right(left_int), right(left_int) <= right_count)
            )
        ),
        left_count == 3
    ]
    solver.add(rules)
    expressions = [
        left_count, right_count,
        left(1), left(2), left(3),
        right(1), right(2), right(3),
    ]
    solves = find_n_solves(expressions, solver)
    for solve in solves:
        print(solve)
        for expr, val in zip(expressions, solve):
            print(expr, val)


def print_out(solver, expressions):
    if solver.check() == z3.sat:
        m = solver.model()
        for var in expressions:
            print(str(var), m.evaluate(var))
    else:
        print("unsat")


if __name__ == '__main__':
    main_part_solve()
