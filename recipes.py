import bisect

import om


class Recipe:

    def __init__(self, part: bytes,
            consumed: list[int] = None,
            produced: list[int] = None,
            conversions: list[tuple[int, int]] = None,
            hovers: list[int] = None):
        self.part: bytes = part
        self.consumed: list[int] = consumed or []
        self.produced: list[int] = produced or []
        self.conversions: list[tuple[int, int]] = conversions or []
        self.hovers: list[int] = hovers or []

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

    def __str__(self):
        str_list = []

        part_str = self.part[6:].decode("UTF-8")

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
                comb_str = "-" + cons_str
            str_list.append(comb_str)

        str_list += [
            "%s -> %s" % (om.Atom.TYPE_NAMES[r], om.Atom.TYPE_NAMES[p])
            for r, p in self.conversions]

        return "Recipe(%s: %s)" % (part_str, " ; ".join(str_list))

    def __repr__(self):
        d = {
            "consumed": self.consumed,
            "produced": self.produced,
            "conversions": self.conversions,
            "hovers": self.hovers
        }
        val_str = ", ".join("%s=%s" % (k, v) for k, v in d.items() if v)

        return "Recipe(%s, %s)" % (self.part, val_str)

    def __hash__(self):
        return hash((self.part, tuple(self.produced), tuple(self.consumed),
        tuple(self.conversions), tuple(self.hovers)))

    def __eq__(self, other):
        return (self.part == other.part
                and self.consumed == other.consumed
                and self.produced == other.produced
                and self.conversions == other.conversions
                and self.hovers == other.hovers)


def recipe_calcify(element: int):
    return Recipe(om.Part.CALCIFICATION,
        conversions=[(element, om.Atom.SALT)])


def recipe_duplicate(element: int):
    return Recipe(om.Part.DUPLICATION,
        conversions=[(om.Atom.SALT, element)], hovers=[element])


def recipe_purify(metal: int, metal_up: int):
    return Recipe(om.Part.PURIFICATION,
        consumed=[metal, metal], produced=[metal_up])


def recipe_project(metal: int, metal_up: int):
    return Recipe(om.Part.PROJECTION,
        consumed=[om.Atom.QUICKSILVER], conversions=[(metal, metal_up)])


def recipe_reagent(type):
    return Recipe(om.Part.INPUT, produced=[type])


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
]

glyph_parts = set(recipe.part for recipe in RECIPES) | {om.Part.BERLO}
Node = tuple[bytes]


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


def search_cheapest_partlist(puzzle: om.Puzzle) -> list[bytes] | None:
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


def expand_possible_types(types: list[int], available_parts: Node = None):
    if available_parts is None:
        available_parts = ()

    hovering_types = []
    typelist = types.copy()
    type_checklist = types.copy()

    available_recipes = [r for r in RECIPES
        if r.part in available_parts]

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
    available_recipes = [r for r in RECIPES if r.part in available_parts]

    # if product size is bigger than all reagent sizes, bonder is required
    if any(
            all(len(r.atoms) < len(p.atoms) for r in puzzle.reagents)
                    for p in puzzle.products
    ):
        required_parts.append(om.Part.BONDER)

    # if product size is smaller than all reagent sizes, debonder is required
    if any(
            all(len(r.atoms) > len(p.atoms) for r in puzzle.reagents)
                    for p in puzzle.products
    ):
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


def constraint_solve_min_cost_partlist(puzzle):
    # discover the smaller recipe list before we define a node since it
    # catches it
    available_parts = puzzle.full_parts_list()
    available_recipes = [r for r in RECIPES if r.part in available_parts]
    available_recipes += [recipe_reagent(t) for t in puzzle.reagent_types()]

    # define the Node
    class Node:

        def __init__(self, partlist: list[bytes] = None,
                type_solves: dict[int, list[Recipe]] = None,
                prev_types: dict[int, list[int]] = None):
            self.partlist = [] if partlist is None else sorted(partlist)
            self.type_solves = type_solves or dict()
            self.prev_types = prev_types or dict()

        def copy(self):
            return Node(self.partlist.copy(), self.type_solves.copy(),
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
                r for r in available_recipes
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

        def neighbors(self):
            unsolved = sorted(self.unsolved(), key=lambda kv: len(kv[1]))
            neighbors = []
            for recipe in unsolved[0][1]:
                copy = self.copy()
                if copy.assign(unsolved[0][0], recipe):
                    neighbors.append(copy)
            return neighbors

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

    # NOW begins the dijkstra solve
    start = Node()
    for type in puzzle.product_types():
        start.update_constraint(type)
    if not start.satisfiable():
        return None

    costs = {start: start.cost()}
    frontier = [start]
    explored = set()

    while frontier:
        current = frontier.pop()
        if current not in explored:
            explored.add(current)
            if current.satisfied():
                partlist = current.partlist
                partlist.remove(om.Part.INPUT)
                return partlist
            for neighbor in current.neighbors():
                costs[neighbor] = neighbor.cost()
                bisect.insort_right(frontier, neighbor, key=lambda n: -costs[n])

    return None
