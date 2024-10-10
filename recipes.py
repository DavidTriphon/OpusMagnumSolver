import bisect

import graph
import om

from collections import Counter


class Recipe:

    def __init__(self, part: bytes | tuple[bytes, int],
            consumed: list[int] = None,
            produced: list[int] = None,
            conversions: list[tuple[int, int]] = None,
            hovers: list[int] = None,
            bonds_consumed: list[int] = None,
            bonds_produced: list[int] = None,
    ):
        self.part: bytes = part
        self.consumed: list[int] = consumed or []
        self.produced: list[int] = produced or []
        self.conversions: list[tuple[int, int]] = conversions or []
        self.hovers: list[int] = hovers or []
        self.bonds_consumed: list[int] = bonds_consumed or []
        self.bonds_produced: list[int] = bonds_produced or []
        self.bond_conversions: list[tuple[int, int]] = []

    def creates(self, element):
        return (element in self.produced or
                element in (c[1] for c in self.conversions))

    def requires(self, element):
        return (element in self.consumed or
                element in (c[0] for c in self.conversions) or
                element in self.hovers)

    def list_created(self):
        return self.produced + list(c[1] for c in self.conversions)

    def list_required(self):
        return self.consumed + list(c for c, _ in self.conversions)

    def atom_metric(self, atom_type):
        return (self.list_created().count(atom_type)
                - self.list_required().count(atom_type))

    def atom_metrics(self):
        positive = Counter(self.list_created())
        negative = Counter(self.list_required())
        for k in negative:
            positive[k] = positive.get(k, 0) - negative[k]
        return positive

    def get(self, is_mass, is_consumed) -> list:
        if is_mass:
            if is_consumed:
                return self.consumed
            else:
                return self.produced
        else:
            if is_consumed:
                return self.bonds_consumed
            else:
                return self.bonds_produced

    def shorthand_str(self):
        if self.part == om.Part.CALCIFICATION:
            return "calc_%s" % om.Atom.TYPE_NAMES[self.conversions[0][0]]
        if self.part == om.Part.DUPLICATION:
            return "dupe_%s" % om.Atom.TYPE_NAMES[self.conversions[0][1]]
        if self.part == om.Part.PURIFICATION:
            return "purify_%s" % om.Atom.TYPE_NAMES[self.consumed[0]]
        if self.part == om.Part.PROJECTION:
            return "project_%s" % om.Atom.TYPE_NAMES[self.conversions[0][0]]
        if self.part == om.Part.ANIMISMUS:
            return "animismus"
        if self.part == om.Part.DISPERSION:
            return "dispersion"
        if self.part == om.Part.UNIFICATION:
            return "unification"
        if self.part == om.Part.BONDER:
            return "bond"
        if self.part == om.Part.UNBONDER:
            return "debond"
        if isinstance(self.part, tuple):
            return "%s_%d" % (self.part[0].decode("UTF-8"), self.part[1])
        return str(self)

    def __str__(self):
        str_list = []

        if isinstance(self.part, bytes):
            if self.part[:6] == b"glyph-":
                part_str = self.part[6:].decode("UTF-8")
            else:
                part_str = self.part.decode("UTF-8")
        else:
            part_str = "%s[%d]" % (self.part[0].decode("UTF-8"), self.part[1])

        if self.hovers:
            hovr_str = '(' + ", ".join([om.Atom.TYPE_NAMES[e] for e in
                self.hovers]) + ')'
            str_list.append(hovr_str)

        if self.consumed or self.produced:
            cons_str = " + ".join(
                [om.Atom.TYPE_NAMES[e] for e in self.consumed])
            prod_str = " + ".join(
                [om.Atom.TYPE_NAMES[e] for e in self.produced])

            comb_str = ""
            if len(self.consumed) > 0 and len(self.produced) > 0:
                comb_str = "%s => %s" % (cons_str, prod_str)
            elif len(self.produced) > 0:
                comb_str = "+" + prod_str
            elif len(self.consumed) > 0:
                comb_str = "-(%s)" % cons_str
            str_list.append(comb_str)

        if self.bonds_consumed:
            str_list.append(", ".join(
                "- %s Bond" % om.Bond.type_name(bond_type)
                    for bond_type in self.bonds_consumed
            ))

        if self.bonds_produced:
            str_list.append(", ".join(
                "+ %s Bond" % om.Bond.type_name(bond_type)
                    for bond_type in self.bonds_produced
            ))

        str_list += [
            "%s -> %s" % (om.Atom.TYPE_NAMES[r], om.Atom.TYPE_NAMES[p])
            for r, p in self.conversions]

        return "Recipe(%s: %s)" % (part_str, " ; ".join(str_list))

    def __repr__(self):
        d = {
            "part": self.part,
            "consumed": self.consumed,
            "produced": self.produced,
            "conversions": self.conversions,
            "hovers": self.hovers,
            "bonds_consumed": self.bonds_consumed,
            "bonds_produced": self.bonds_produced,
        }
        val_str = ", ".join((
            "%s=%s" % (k, v)
            for k, v in d.items()
            if v is not None
        ))

        return "Recipe(%s, %s)" % (self.part, val_str)

    def __hash__(self):
        return hash((
            self.part,
            tuple(self.produced), tuple(self.consumed),
            tuple(self.conversions), tuple(self.hovers),
            tuple(self.bonds_produced), tuple(self.bonds_consumed)
        ))

    def __eq__(self, other):
        return (
                self.part == other.part
                and self.consumed == other.consumed
                and self.produced == other.produced
                and self.conversions == other.conversions
                and self.hovers == other.hovers
                and self.bonds_consumed == other.bonds_consumed
                and self.bonds_produced == other.bonds_produced
        )


def recipe_calcify(element: int):
    return Recipe(part=om.Part.CALCIFICATION,
        conversions=[(element, om.Atom.SALT)])


def recipe_duplicate(element: int):
    return Recipe(part=om.Part.DUPLICATION,
        conversions=[(om.Atom.SALT, element)], hovers=[element])


def recipe_purify(metal: int, metal_up: int):
    return Recipe(part=om.Part.PURIFICATION,
        consumed=[metal, metal], produced=[metal_up])


def recipe_project(metal: int, metal_up: int):
    return Recipe(part=om.Part.PROJECTION,
        consumed=[om.Atom.QUICKSILVER], conversions=[(metal, metal_up)])


def recipe_reagent(molecule: om.Molecule, molecule_index: int = None):
    return Recipe(part=(om.Part.INPUT, molecule_index),
        produced=[atom.type for atom in molecule.atoms],
        bonds_produced=[bond.type for bond in molecule.bonds]
    )


def puzzles_reagent_recipes(puzzle):
    return [recipe_reagent(reagent, i) for i, reagent in enumerate(
        puzzle.reagents)]


def recipe_product(molecule: om.Molecule, molecule_index: int = None):
    return Recipe(part=(om.Part.OUTPUT_STANDARD, molecule_index),
        consumed=[atom.type for atom in molecule.atoms],
        bonds_consumed=[bond.type for bond in molecule.bonds]
    )


def puzzles_product_recipes(puzzle):
    return [recipe_product(product, i) for i, product in enumerate(
        puzzle.products)]


def recipe_disposal(atom_type: int):
    return Recipe(part=om.Part.DISPOSAL, consumed=[atom_type])


RECIPES: list[Recipe] = [
    # # elements
    # calcifications
    *[recipe_calcify(element) for element in om.Atom.TYPES_ELEMENTAL],
    # duplications
    *[recipe_duplicate(element) for element in om.Atom.TYPES_ELEMENTAL],
    # polarizations
    Recipe(om.Part.ANIMISMUS,
        consumed=[om.Atom.SALT, om.Atom.SALT],
        produced=[om.Atom.VITAE, om.Atom.MORS]),
    # quintessence
    Recipe(om.Part.UNIFICATION,
        consumed=[om.Atom.AIR, om.Atom.EARTH, om.Atom.FIRE, om.Atom.WATER],
        produced=[om.Atom.QUINTESSENCE]),
    Recipe(om.Part.DISPERSION,
        consumed=[om.Atom.QUINTESSENCE],
        produced=[om.Atom.AIR, om.Atom.EARTH, om.Atom.FIRE, om.Atom.WATER]),

    # # metals
    # purification
    recipe_purify(om.Atom.LEAD, om.Atom.TIN),
    recipe_purify(om.Atom.TIN, om.Atom.IRON),
    recipe_purify(om.Atom.IRON, om.Atom.COPPER),
    recipe_purify(om.Atom.COPPER, om.Atom.SILVER),
    recipe_purify(om.Atom.SILVER, om.Atom.GOLD),
    # projection
    recipe_project(om.Atom.LEAD, om.Atom.TIN),
    recipe_project(om.Atom.TIN, om.Atom.IRON),
    recipe_project(om.Atom.IRON, om.Atom.COPPER),
    recipe_project(om.Atom.COPPER, om.Atom.SILVER),
    recipe_project(om.Atom.SILVER, om.Atom.GOLD),

    # # etc
    # disposal
    *[recipe_disposal(atom_type) for atom_type in om.Atom.TYPES],

    # BOND RECIPES
    Recipe(part=om.Part.BONDER, bonds_produced=[om.Bond.NORMAL]),
    Recipe(part=om.Part.UNBONDER, bonds_consumed=[om.Bond.NORMAL]),
    Recipe(part=om.Part.TRIPLEX, bonds_produced=[om.Bond.TRIPLEX]),
]


def all_puzzle_recipes(puzzle):
    return (
            usable_recipes(puzzle.reagent_types(), puzzle.full_parts_list())
            + puzzles_reagent_recipes(puzzle)
            + puzzles_product_recipes(puzzle)
    )


def part_recipes(parts) -> list[Recipe]:
    return [r for r in RECIPES if r.part in parts]


def usable_recipes(given_types, parts) -> list[Recipe]:
    maximum_types = expand_possible_types(given_types, parts)
    available_recipes = part_recipes(parts)
    return [
        recipe for recipe in available_recipes
        if all(req in maximum_types for req in recipe.list_required())
    ]


def search_cheapest_partlist(puzzle: om.Puzzle) -> list[bytes] | None:
    glyph_parts = set(recipe.part for recipe in RECIPES) | {om.Part.BERLO}
    Node: type = tuple[bytes]

    def neighbors(selected=(), available=None) -> list[Node]:
        if available is None:
            available = glyph_parts
        else:
            available = [a for a in available if a in glyph_parts]

        nbs = []
        for part in available:
            if part not in selected:
                partlist = tuple(list(selected) + [part])
                nbs.append(partlist)
        return nbs

    product_types = puzzle.product_types()
    reagent_types = puzzle.reagent_types()
    available_parts = puzzle.full_parts_list()
    start: Node = tuple(search_required_partlist(puzzle))
    node_cost: dict[Node, int] = {start: 0}
    node_types: dict[Node, list[int]] = {
        start: expand_possible_types(
            reagent_types, start)
    }
    frontier: list[Node] = [start]
    explored: set[Node] = set()

    while frontier:
        current = frontier.pop()
        if current not in explored:
            explored.add(current)
            if all(p in node_types[current] for p in product_types):
                return list(current)
            for neighbor in neighbors(current, available_parts):
                node_cost[neighbor] = sum(om.Part.COSTS[p] for p in neighbor)
                node_types[neighbor] = expand_possible_types(
                    node_types[current], neighbor)
                bisect.insort_right(frontier, neighbor, key=lambda n:
                -node_cost[n])

    return None


def expand_possible_types(types: list[int],
        available_parts: tuple[bytes] = None):
    if available_parts is None:
        available_parts = ()

    hovering_types = []
    typelist = types.copy()
    type_checklist = types.copy()

    available_recipes = part_recipes(available_parts)

    if om.Part.BERLO in available_parts:
        hovering_types = om.Atom.TYPES_ELEMENTAL

    while type_checklist:
        next_type = type_checklist.pop()
        for recipe in available_recipes:
            if recipe.requires(next_type):
                if all([ing in typelist for ing in recipe.list_required()]):
                    if all([
                        hov in hovering_types or hov in typelist
                        for hov in recipe.hovers
                    ]):
                        for product in recipe.list_created():
                            if product not in typelist:
                                typelist.append(product)
                                type_checklist.append(product)

    return typelist


def search_required_partlist(puzzle: om.Puzzle) -> list[bytes] | None:
    required_parts = []
    available_parts = puzzle.full_parts_list()
    available_recipes = usable_recipes(puzzle.reagent_types(), available_parts)

    # if product size is bigger than all reagent sizes, bonder is required
    if is_bonder_required(puzzle) is True:
        required_parts.append(om.Part.BONDER)

    # if product size is smaller than all reagent sizes, debonder is required
    if is_debonder_required(puzzle) is True:
        required_parts.append(om.Part.UNBONDER)

    # if an atom type exists in the product and not in the reagents, and
    # there is only one available part that makes it, that part is required.
    required_types = [(t, ()) for t in puzzle.product_types()]
    reagent_types = puzzle.reagent_types()

    while required_types:
        p_type, prev_types = required_types.pop()
        if p_type not in reagent_types:
            p_recipes = [
                r for r in available_recipes
                if p_type in r.list_created()
                   and all(t not in prev_types for t in r.list_required())
            ]
            recipe_parts = set(r.part for r in p_recipes)
            if len(recipe_parts) == 1:
                if len(p_recipes) == 1:
                    for r_type in p_recipes[0].list_required():
                        # this needs some kind of recursion check to avoid
                        # infinite looping in the case unsatisfiable loops
                        required_types.append(
                            (r_type, (p_type,) + prev_types))
                    for hover_type in p_recipes[0].hovers:
                        # this is always true
                        if hover_type not in reagent_types:
                            if om.Part.BERLO not in required_parts:
                                required_parts.append(om.Part.BERLO)
                        else:
                            raise NotImplementedError(
                                "Somehow a hover was needed that was "
                                "already satisfied by the reagents...")
                part = list(recipe_parts)[0]
                if part not in required_parts:
                    required_parts.append(part)
            elif len(p_recipes) == 0:
                # invalid situation, atom type not satisfiable
                return None

    # for any given atom type in a product, if all the reagents that
    # contain that atom type have more atoms of that type than the
    # product, a debonder is required.
    # This is not true if a conversion recipe is added that creates more
    # reagents that could contain that atom type. So this rule only holds if
    # there are no more recipes to add that could possibly reveal this atom
    # type.

    # for any given atom type in a product, if all the reagents that
    # contain that atom type have fewer atoms of that type than the
    # product, a bonder is required.
    # This is not true if a conversion recipe is added that creates more
    # reagents that could contain that atom type. So this rule only holds if
    # there are no more recipes to add that could possibly reveal this atom
    # type.

    return required_parts


class CostConstraintNode(graph.BaseNode):

    def __init__(self, available_recipes, partlist: list[bytes] = None,
            type_solves: dict[int, list[Recipe]] = None,
            prev_types: dict[int, list[int]] = None):
        self.available_recipes = available_recipes
        self.partlist = [] if partlist is None else sorted(partlist)
        self.type_solves = type_solves or dict()
        self.prev_types = prev_types or dict()

    def copy(self):
        return CostConstraintNode(
            self.available_recipes,
            self.partlist.copy(),
            self.type_solves.copy(),
            self.prev_types.copy())

    def cost(self):
        return sum(om.Part.COSTS[part] for part in self.partlist)

    def unsolved(self):
        return [(k, v) for k, v in self.type_solves.items() if len(v) > 1]

    def satisfied(self):
        return all(len(v) == 1 for _, v in self.type_solves.items())

    def satisfiable(self):
        return all(len(v) > 0 for _, v in self.type_solves.items())

    def available_solves_for(self, type):
        return [
            r for r in self.available_recipes
            if type in r.list_created()
               and all(
                t not in self.prev_types.get(type, {})
                    for t in r.list_required()
            )
        ]

    def assign(self, type, recipe=None):
        if recipe is not None:
            if recipe not in self.type_solves[type]:
                raise ValueError(
                    "This recipe is not in the available options")
            self.type_solves[type] = [recipe]
        else:
            # recipe=None means we think there is already only 1 recipe.
            if len(self.type_solves[type]) != 1:
                raise ValueError(
                    "recipe is None but there is not only 1 option.")
            recipe = self.type_solves[type][0]
        self.add_part(recipe.part)
        if len(recipe.hovers) > 0:
            self.add_part(om.Part.BERLO)
        for new_type in recipe.list_required():
            if new_type in self.type_solves.keys():
                if new_type in self.prev_types[type]:
                    raise ValueError("%s loops" % type)
            # always update if key is in or not
            if not self.update_constraint(new_type, type):
                return False
        return self

    def update_constraint(self, type, prev_type=None):
        all_prev_types = [] if prev_type is None else \
            self.prev_types[prev_type] + [prev_type]
        if type not in self.prev_types:
            self.prev_types[type] = []
        self.prev_types[type] = self.prev_types[type] + all_prev_types

        self.type_solves[type] = self.available_solves_for(type)
        if len(self.type_solves[type]) == 0:
            return False
        if len(self.type_solves[type]) == 1:
            if not self.assign(type):
                return False
        return self

    def add_part(self, part):
        if part not in self.partlist:
            bisect.insort(self.partlist, part)
        return self

    def neighbors_edges(self):
        unsolved = sorted(self.unsolved(), key=lambda kv: len(kv[1]))
        neighbors = []
        for recipe in unsolved[0][1]:
            copy = self.copy()
            if copy.assign(unsolved[0][0], recipe):
                neighbors.append(copy)
        edges = [None for i in range(len(neighbors))]
        return neighbors, edges

    def __hash__(self):
        return hash((
            tuple(self.partlist),
            tuple(sorted(
                (k, tuple(v)) for k, v in self.type_solves.items()
            )),
            tuple(sorted(
                (k, tuple(v)) for k, v in self.prev_types.items()
            ))
        ))

    def __eq__(self, other):
        return (self.partlist == other.partlist and
                self.type_solves == other.type_solves and
                self.prev_types == other.prev_types)


def constraint_solve_min_cost_partlist(puzzle):
    # discover the smaller recipe list before we define a node since it
    # catches it
    available_parts = puzzle.full_parts_list()
    available_recipes = usable_recipes(puzzle.reagent_types(), available_parts)
    available_recipes += puzzles_reagent_recipes(puzzle)

    start = CostConstraintNode(available_recipes)
    for type in puzzle.product_types():
        start.update_constraint(type)
    if not start.satisfiable():
        return None

    end, cost = start.node_search(
        goal_condition=lambda node: node.satisfied(),
        cost_func=lambda node: node.cost(),
    )

    return end.partlist


def is_bonder_required(puzzle: om.Puzzle):
    if all(len(p.atoms) == 1 for p in puzzle.products):
        return False

    if any((
            all(len(r.atoms) < len(p.atoms) for r in puzzle.reagents)
            for p in puzzle.products
    )):
        return True

    # effectively "MAYBE / UNKNOWN / UNCERTAIN"
    return None


def is_debonder_required(puzzle: om.Puzzle):
    if all(len(r.atoms) == 1 for r in puzzle.reagents):
        return False

    if any((
            all(len(r.atoms) > len(p.atoms) for r in puzzle.reagents)
            for p in puzzle.products
    )):
        return True

    # effectively "MAYBE / UNKNOWN / UNCERTAIN"
    return None
