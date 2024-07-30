import analysis
import om


class Recipe:

    def __init__(self,
            ingredients: list[int],
            products: list[int],
            part: bytes,
            hovering: list[int] = None
    ):
        self.ingredients: list[int] = ingredients
        self.products: list[int] = products
        self.part: bytes = part
        if hovering is None:
            hovering: list[int] = []
        self.hovering = hovering


def recipe_calcify(element: int):
    return Recipe(
        [element],
        [om.Atom.SALT],
        om.Part.CALCIFICATION
    )


def recipe_duplicate(element: int):
    return Recipe(
        [om.Atom.SALT],
        [element],
        om.Part.DUPLICATION,
        [element]
    )


def recipe_purify(metal: int, metal_up: int):
    return Recipe(
        [metal, metal],
        [metal_up],
        om.Part.PURIFICATION
    )


def recipe_project(metal, metal_up):
    return Recipe(
        [metal, om.Atom.QUICKSILVER],
        [metal_up],
        om.Part.PROJECTION
    )


RECIPES: list[Recipe] = [
    # # elements
    # calcifications
    *[recipe_calcify(element) for element in analysis.ELEMENTAL_ATOM_TYPES],
    # duplications
    *[recipe_duplicate(element) for element in analysis.ELEMENTAL_ATOM_TYPES],
    # polarizations
    Recipe(
        [om.Atom.SALT, om.Atom.SALT],
        [om.Atom.VITAE, om.Atom.MORS],
        om.Part.ANIMISMUS
    ),
    # quintessence
    Recipe(
        [om.Atom.AIR, om.Atom.EARTH, om.Atom.FIRE,
            om.Atom.WATER],
        [om.Atom.QUINTESSENCE],
        om.Part.UNIFICATION
    ),
    Recipe(
        [om.Atom.QUINTESSENCE],
        [om.Atom.AIR, om.Atom.EARTH, om.Atom.FIRE, om.Atom.WATER],
        om.Part.DISPERSION
    ),

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
