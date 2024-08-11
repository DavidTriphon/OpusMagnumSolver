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


def search_min_partlist(puzzle: om.Puzzle) -> list[bytes] | None:
    product_types = puzzle.product_types()
    reagent_types = puzzle.reagent_types()
    available_parts = puzzle.full_parts_list()
    start: Node = ()
    node_cost: dict[Node, int] = {start: 0}
    node_types: dict[Node, list[int]] = {start: reagent_types}
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


def list_possible_theory_optimal_assembly_plans(puzzle):
    # WIP of a method for recursively searching for ways to construct the
    # product molecules
    max_partlist = puzzle.full_parts_list()
    available_recipes = [r for r in RECIPES if r.part in max_partlist]
    atomtype_recipes = dict()

    for product in puzzle.products:
        for atom in product.atoms:
            if atom.type not in atomtype_recipes:
                atomtype_recipes[atom.type] = [r for r in RECIPES
                    if atom.type in r.creates()]
            available_recipes = atomtype_recipes[atom.type]
