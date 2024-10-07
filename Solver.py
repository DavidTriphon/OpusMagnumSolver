import sys
from collections import defaultdict
from enum import Enum

import graph
import om
import recipes

import logging

log = logging.getLogger()
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(
    logging.Formatter('[%(levelname)s][%(filename)s:%(lineno)s - '
                      '%(funcName)20s() ] %(message)s'))
handler.setLevel(logging.INFO)
log.addHandler(handler)


def constraint_satisfiable(constraints):
    return len(constraints) > 0


def constraint_assign(constraints, remaining_value):
    for value in constraints.copy():
        if value is not remaining_value:
            constraints.remove(value)
    return constraint_satisfiable(constraints)


def constraint_remove(constraints, value_to_remove):
    constraints.remove(value_to_remove)
    return constraint_satisfiable(constraints)


def list_range_inclusive(maximum):
    return list(range(0, maximum + 1))


PART_ACCESSES = {
    # these lists represent the default sets of possible numbers of hexes of a
    # part that require direct access from an arm.

    # bonders and debonders could be 0-access or full-access
    om.Part.BONDER: list_range_inclusive(2),
    om.Part.UNBONDER: list_range_inclusive(2),
    # a triplex has 3 hexes, but will never require direct access
    # to all 3 because it can rotate the molecule after the first bond.
    om.Part.TRIPLEX: list_range_inclusive(2),
    om.Part.CALCIFICATION: list_range_inclusive(1),
    # these parts will always require direct access to all hexes
    om.Part.DISPERSION: [5],
    om.Part.UNIFICATION: [5],
    om.Part.ANIMISMUS: [4],
    om.Part.PURIFICATION: [3],
    # disposal may take up more direct access hexes because of its size,
    # but it only ever requires 1 same exact hex to be accessible.
    om.Part.DISPOSAL: [1],
    # duplication could theoretically be done with a single atom, but most
    # times it will be unnecessary.
    om.Part.DUPLICATION: list_range_inclusive(2),
    # projection always consumes 1 QS, so it's never a 0-access
    om.Part.PROJECTION: [1, 2],
    # input and output parts are never0-access, but could warrant access to
    # separate hexes possibly? unsure.
    om.Part.INPUT: [1],
    om.Part.OUTPUT_STANDARD: [1],
    om.Part.OUTPUT_REPEATING: [1],
}


class TrackShape(Enum):
    """
    a 2 or 3 length arm config will require at most 4 tracks for access
    shapes, in the shape of a diamond.
    the animismus requires a 3 hex triangle, and the unification glyph
    requires a 3 hex C-shape. When added together, they form the
    animismus's required 4 hex diamond.
    \0\1\
     \2\3\

    this means there are 5 different track configurations for 2
    length arms:
    (1 hex single),
    (2 hex adjacent),
    (3 hex C-shape), (3 hex triangle),
    (4 hex diamond)
    where each configuration satisfies the other requirements of any
    config in a row above, but not a config in the same row or below.

    for length 3 arms: there are additional track arrangements for the
    unification and dispersion glyphs:
    (5 hourglass) (for the unification), (5, trapezoid) (for dispersal)
    (6 bent) (for both unification and dispersal)
    """
    ONE = (1, "single", 1)
    TWO = (2, "adjacent", 1)
    THREE_TRIANGLE = (3, "triangle", 2)
    THREE_CURVE = (3, "curve", 2)
    FOUR_DIAMOND = (4, "diamond", 2)
    FIVE_HOURGLASS = (5, "hourglass", 3)
    FIVE_TRAPEZOID = (5, "trapezoid", 3)
    SIX_BENT = (6, "bent", 3)

    def __init__(self, hexes, shape, min_arm_length):
        self.hexes = hexes
        self.shape = shape
        self.min_arm_length = min_arm_length

    def accessible_hexes(self, armlength):
        if armlength == 1:
            if self.hexes == 1:
                return 6
            elif self.hexes == 2:
                return 8
        if armlength == 2 or armlength == 3:
            return 6 * self.hexes
        raise ValueError(
            "invalid armlength (%r) for shape (%r)" % (armlength, self))

    def satisfies(self, other):
        # if self is bigger than other, we satisfy
        # NOT THE OTHER WAY AROUND
        if self.hexes > other.hexes:
            return True
        if self.hexes < other.hexes:
            return False
        return self.shape == other.shape

    def track_cost(self):
        return 0 if self.hexes <= 1 else 5 * self.hexes


ZER0_OR_ONE_ACCESS_TRACK_SHAPE = {
    1: TrackShape.ONE,
    2: TrackShape.ONE,
    3: TrackShape.ONE,
}
TWO_ACCESS_TRACK_SHAPE = {
    1: TrackShape.ONE,
    2: TrackShape.TWO,
    3: TrackShape.TWO,
}


def default_min_shapes(accesses):
    defaults = {
        0: ZER0_OR_ONE_ACCESS_TRACK_SHAPE,
        1: ZER0_OR_ONE_ACCESS_TRACK_SHAPE,
        2: TWO_ACCESS_TRACK_SHAPE,
    }
    return {
        access: defaults[access]
        for access in accesses
    }


# indexed by [part] [accesses] [arm_length]
# returns min track shape
PART_ACCESS_LENGTH_MIN_TRACK_SHAPES = {
    om.Part.BONDER: default_min_shapes(PART_ACCESSES[om.Part.BONDER]),
    om.Part.UNBONDER: default_min_shapes(PART_ACCESSES[om.Part.UNBONDER]),
    om.Part.TRIPLEX: default_min_shapes(PART_ACCESSES[om.Part.TRIPLEX]),
    om.Part.CALCIFICATION:
        default_min_shapes(PART_ACCESSES[om.Part.CALCIFICATION]),
    om.Part.DUPLICATION: default_min_shapes(PART_ACCESSES[om.Part.DUPLICATION]),
    om.Part.PROJECTION: default_min_shapes(PART_ACCESSES[om.Part.PROJECTION]),
    # now for one to one constraints
    om.Part.DISPERSION: {
        5: {
            1: None,
            2: TrackShape.FOUR_DIAMOND,
            3: TrackShape.FIVE_TRAPEZOID,
        }
    },
    om.Part.UNIFICATION: {
        5: {
            1: None,
            2: TrackShape.THREE_CURVE,
            3: TrackShape.FIVE_HOURGLASS,
        }
    },
    om.Part.ANIMISMUS: {
        4: {
            1: None,
            2: TrackShape.FOUR_DIAMOND,
            3: TrackShape.FOUR_DIAMOND,
        }
    },
    om.Part.PURIFICATION: {
        3: {
            1: None,
            2: TrackShape.THREE_TRIANGLE,
            3: TrackShape.THREE_TRIANGLE,
        }
    },
    # one to one constraints for 1-access only parts
    om.Part.DISPOSAL: {
        # requires only one access but cannot allow length 1 arms
        1: {
            1: None,
            2: TrackShape.ONE,
            3: TrackShape.ONE,
        }
    },
    # these are generic single access, the hex is variable
    om.Part.INPUT: default_min_shapes([1]),
    om.Part.OUTPUT_STANDARD: default_min_shapes([1]),
    om.Part.OUTPUT_REPEATING: default_min_shapes([1]),
}


class RecipeApplicationConstraints:

    def __init__(self, mass_produced, mass_consumed, bonds_produced,
            bonds_consumed):
        # key = atom_index
        # value = set of possible paired (recipe, iteration_id, atom_index)
        self.mass_produced = mass_produced or []
        self.mass_consumed = mass_consumed or []
        # key = bond_index
        # value = set of possible paired (recipe, iteration_id, bond_index)
        self.bonds_produced = bonds_produced or []
        self.bonds_consumed = bonds_consumed or []

    def get(self, is_mass, is_consumed):
        if is_mass:
            if is_consumed:
                return self.mass_consumed
            else:
                return self.mass_produced
        else:
            if is_consumed:
                return self.bonds_consumed
            else:
                return self.bonds_produced


class SolutionNode(graph.BaseNode):
    """
    each instance represents a set of possible solutions for the problem
    divide each set into separate but complete subsets in order to
    branch and bound the possible solutions to find an optimal solution
    """

    def __init__(self, puzzle):
        # constants, do not deep copy
        self.puzzle = puzzle
        self.applicable_recipes = (
                recipes.usable_recipes(puzzle.reagent_types(),
                    puzzle.full_parts_list()) +
                recipes.puzzles_reagent_recipes(puzzle) +
                recipes.puzzles_product_recipes(puzzle)
        )
        usable_parts = {
            part
            for part in (
                    set(puzzle.full_parts_list()) - om.Part.PARTS_ARMS -
                    {om.Part.MULTIBONDER, om.Part.EQUILIBRIUM}
            )
            if any((
                recipe in self.applicable_recipes
                for recipe in recipes.part_recipes(part)
            ))
        }

        ### this entire next section is dedicated to caching the simplest
        # paths between mass types, and the correponding available recipes
        # that are required for those simplest paths.

        # eg. if you choose that your fire atom comes from a water atom,
        # this declares that you need a duplicator and a calcifier.

        recipes_produce_mass = {
            mass_type: {
                (recipe, None, index)
                for recipe in self.applicable_recipes
                for index, matching_mass in enumerate(recipe.produced)
                if mass_type == matching_mass
            }
            for mass_type in om.Atom.TYPES
        }
        recipes_consume_mass = {
            mass_type: {
                (recipe, None, index)
                for recipe in self.applicable_recipes
                for index, matching_mass in enumerate(recipe.consumed)
                if mass_type == matching_mass
            }
            for mass_type in om.Atom.TYPES
        }
        self.valid_bond_starts = {
            bond_type: {
                (recipe, None, index)
                for recipe in self.applicable_recipes
                for index, matching_bond in enumerate(recipe.bonds_produced)
                if bond_type == matching_bond
            }
            for bond_type in om.Bond.TYPES
        }
        self.valid_bond_ends = {
            bond_type: {
                (recipe, None, index)
                for recipe in self.applicable_recipes
                for index, matching_bond in enumerate(recipe.bonds_consumed)
                if bond_type == matching_bond
            }
            for bond_type in om.Bond.TYPES
        }

        mass_conversion_edges = defaultdict(dict)
        for recipe in self.applicable_recipes:
            for start_mass, end_mass in recipe.conversions:
                mass_conversion_edges[start_mass][end_mass] = recipe
        self.mass_conversion_matrix = dict()
        for mass_type in om.Atom.TYPES:
            self.mass_conversion_matrix[mass_type] = {mass_type: []}
            to_check = {mass_type}
            while to_check:
                checking = to_check.pop()
                for end_type, recipe in mass_conversion_edges[checking].items():
                    if end_type not in self.mass_conversion_matrix[mass_type]:
                        self.mass_conversion_matrix[mass_type][end_type] = \
                            self.mass_conversion_matrix[mass_type][checking] + [
                                recipe]
                        to_check.add(end_type)
        reverse_mass_conversion_matrix = defaultdict(list)
        for start, ends in self.mass_conversion_matrix.items():
            for end in ends.keys():
                reverse_mass_conversion_matrix[end].append(start)

        self.valid_mass_ends = {
            start: {
                option
                for end in ends.keys()
                for option in recipes_consume_mass[end]
            }
            for start, ends in mass_conversion_edges.items()
        }
        self.valid_mass_starts = {
            end: {
                option
                for start in starts
                for option in recipes_produce_mass[start]
            }
            for end, starts in reverse_mass_conversion_matrix.items()
        }

        # back to your regularly scheduled constraints

        # recipe metrics
        self.byproduct_handling_method = [
            "wastechain", "disposer", "waste cells", "no byproduct"
        ]

        # recipe : [recipeApplicationConstraints obj]
        self.recipe_applications = {
            recipe: []
            for recipe in self.applicable_recipes
        }
        # recipe : recipeApplicationConstraints obj (not a list)
        self.recipes_received = {
            recipe: RecipeApplicationConstraints(
                [set() for _ in recipe.produced],
                [set() for _ in recipe.consumed],
                [set() for _ in recipe.bonds_produced],
                [set() for _ in recipe.bonds_consumed],
            )
            for recipe in self.applicable_recipes
        }
        # recipe : int
        self.recipe_count_bounds = {
            recipe: [0, None]
            for recipe in self.applicable_recipes
        }

        ## COST SPECIFIC CONSTRAINTS

        # PARTS
        self.reagents_used = {
            i: [True, False]
            for i in range(len(puzzle.reagents))
        }
        if len(self.reagents_used) == 1:
            self.reagents_used[0] = [True]

        self.parts_used = {
            part: [True, False]
            for part in usable_parts
        }
        # populated as parts are marked as required with the default values
        # for that part
        self.parts_accesses = dict()

        # set of valid armlength and trackshape combinations
        self.armlength_trackshapes = {
            (arm_length, shape)
            for arm_length in [1, 2, 3]
            for shape in list(TrackShape)
            if shape.min_arm_length <= arm_length
        }

        self.propagate_puzzle_to_partlist()
        for recipe in self.applicable_recipes:
            if recipe.part == om.Part.OUTPUT_STANDARD:
                if not self.assign_recipe_bound_min(recipe, 1):
                    raise ValueError("contradiction raised")

    def copy(self):
        # TODO copy
        pass

    def __hash__(self):
        # TODO hash
        pass

    # === methods of deterministic nature ===

    def min_cost(self):
        # for arm cost is 20
        cost = 20
        # min track cost
        cost += min(
            shape.track_cost() for _, shape in self.armlength_trackshapes)
        # add parts cost
        for part, part_constraints in self.parts_used.items():
            if part_constraints == [True]:
                cost += om.Part.COSTS[part]
        return cost

    def min_accesses(self):
        min_part_accesses = sum(
            min(accesses) for accesses in self.parts_accesses.values()
        )
        min_reagent_accesses = len([
            i for i, constraints in self.reagents_used.items()
            if constraints == [True]
        ]) or 1  # defaulting to zero is silly
        product_accesses = len(self.puzzle.products)
        return min_part_accesses + min_reagent_accesses + product_accesses

    def are_all_parts_constrained(self):
        parts_constrained = all(
            len(used) == 1
                for part, used in self.parts_used
        )
        return parts_constrained

    # === methods for checking how constraints propagate against each other ===

    def propagate_puzzle_to_partlist(self):
        # search for the minimum possible partlist
        required_partlist = recipes.search_required_partlist(self.puzzle)
        isValid = True
        for part in required_partlist:
            isValid &= self.assign_part_usage(part, True)
        # if bonder is activated, check if it's immediately obvious to be
        # 2-access
        if self.parts_used[om.Part.BONDER] == [True]:
            if all(len(reagent.atoms) == 1 for reagent in self.puzzle.reagents):
                self.assign_part_accesses(om.Part.BONDER, 2)

        return isValid

    def propagate_recipe_applications_counts(self, recipe):
        # call me when the minimum number of recipe applications increases
        recipe_bounds = self.recipe_count_bounds[recipe]
        while len(self.recipe_applications[recipe]) < recipe_bounds[0]:
            if not self.add_recipe_application_constraints(recipe):
                log.info("contradiction seen (recipe=%r", recipe)
                return False
        return True

    def propagate_recipe_application_options(self, recipe, id):
        # call me when a recipe application has been added or updated
        for is_mass in [True, False]:
            for is_consumed in [True, False]:
                log.debug("new loop %r %r %r", recipe, is_mass,
                    is_consumed)
                for i in range(len(recipe.get(is_mass, is_consumed))):
                    if not self.propagate_recipe_application_single(
                            recipe, id, is_mass, is_consumed, i):
                        log.info("contradiction seen (recipe=%r)",
                            recipe)
                        return False
        return True

    def propagate_recipe_application_single(self, recipe, id, is_mass,
            is_consumed, this_index):
        log.debug("propagate_recipe_application_single")
        # call this when the number of options for this has been reduced.
        options = self.recipe_applications[recipe][id].get(
            is_mass, is_consumed)[this_index]
        if len(options) == 1:
            # options are constrained, propagate:
            other_recipe, other_id, other_index = \
                next(iter(options))
            other_type = other_recipe.get(is_mass, not is_consumed)[other_index]
            this_type = recipe.get(is_mass, is_consumed)[this_index]
            if is_consumed:
                required_recipes = self.mass_conversion_matrix[other_type][
                    this_type]
            else:
                required_recipes = self.mass_conversion_matrix[this_type][
                    other_type]
            for required_recipe in required_recipes:
                if not self.assign_recipe_usage(required_recipe, True):
                    log.info("contradiction seen (recipe=%r, index=%d)",
                        recipe, this_index)
                    return False
            if other_id is None:
                # TODO: constrain the other side when the ID is known
                if not self.assign_recipe_received(other_recipe,
                        is_mass, not is_consumed, other_index,
                        (recipe, id, this_index)):
                    log.info("contradiction seen (recipe=%r, index=%r)",
                        recipe, this_index)
                    return False
        return True

    def propagate_min_access_to_track_count(self):
        # call me when min access increases, or when an access value is
        # eliminated

        # propagates min accesses from part access sums to limit the track
        # shapes based on
        min_accesses = self.min_accesses()
        invalid_armlength_trackshapes = [
            (armlength, trackshape)
            for armlength, trackshape, in self.armlength_trackshapes
            if trackshape.accessible_hexes(armlength) < min_accesses
        ]
        for armlength_trackshape in invalid_armlength_trackshapes:
            if not self.eliminate_armlength_trackshape(armlength_trackshape):
                return False
        return True

    def propagate_part_access_to_shapes(self, part):
        # call me when any part access increases, or when an access value is
        # eliminated

        # propagate a part and their access constraints to see if the track
        # shape is impacted
        part_access = self.parts_accesses[part]
        min_trackshape_for_armlength = \
            PART_ACCESS_LENGTH_MIN_TRACK_SHAPES[part][min(part_access)]
        invalid_armlength_trackshapes = [
            (armlength, trackshape)
            for armlength, trackshape in self.armlength_trackshapes
            if not trackshape.satisfies(min_trackshape_for_armlength[armlength])
        ]
        for armlength_trackshape in invalid_armlength_trackshapes:
            if not self.eliminate_armlength_trackshape(armlength_trackshape):
                return False
        return True

    # === methods for assigning and eliminating constraints according to
    # branching or to constraint propagation ===

    def add_recipe_application_constraints(self, recipe):
        new_application = RecipeApplicationConstraints(
            [
                {
                    (dest_recipe, id, index)
                    for (dest_recipe, id, index)
                    in self.valid_mass_ends[atom_type]
                    if self.recipe_count_bounds[dest_recipe][1] != 0
                }
                for atom_type in recipe.produced
            ],
            [
                {
                    (dest_recipe, id, index)
                    for (dest_recipe, id, index)
                    in
                    self.valid_mass_starts[atom_type]
                    if self.recipe_count_bounds[dest_recipe][1] != 0
                }
                for atom_type in recipe.consumed
            ],
            [
                {
                    (dest_recipe, id, index)
                    for (dest_recipe, id, index)
                    in
                    self.valid_bond_ends[bond_type]
                    if self.recipe_count_bounds[dest_recipe][1] != 0
                }
                for bond_type in recipe.bonds_produced
            ],
            [
                {
                    (dest_recipe, id, index)
                    for (dest_recipe, id, index)
                    in self.valid_bond_starts[bond_type]
                    if self.recipe_count_bounds[dest_recipe][1] != 0
                }
                for bond_type in recipe.bonds_consumed
            ]
        )
        self.recipe_applications[recipe].append(new_application)
        # check and propagate single options
        if not self.propagate_recipe_application_options(recipe,
                len(self.recipe_applications[recipe]) - 1):
            log.info("contradiction seen (recipe=%r)", recipe)
            return False

        return True

    def assign_recipe_usage(self, recipe, is_used):
        bounds = self.recipe_count_bounds[recipe]
        if is_used:
            # usages is at least 1 if used
            if bounds[1] == 0:
                log.info("contradiction raised for assigning recipe usage to "
                         "an unusable recipe")
                return False
            if bounds[0] == 0:
                bounds[0] = 1
            return True
        else:
            # usages is 0 if used
            if bounds[0] > 0:
                log.info("contradiction raised for disabling recipe usage on "
                         "a used recipe")
                return False
            if bounds[1] > 0 or bounds[1] is None:
                bounds[1] = 0
            return True

    def assign_recipe_bound_min(self, recipe, mini):
        current_min = self.recipe_count_bounds[recipe][0]
        if current_min == mini:
            return True
        if current_min > mini:
            log.info("contradiction raised for assigning a in incompatible "
                     "minimum")
            return False
        self.recipe_count_bounds[recipe][0] = mini
        if not self.propagate_recipe_applications_counts(recipe):
            log.info("contradiction seen (recipe=%r)", recipe)
            return False
        return True

    def assign_recipe_received(self, recipe, is_mass,
            is_consumed, index, sender):
        receiver = self.recipes_received[recipe].get(
            is_mass, is_consumed)[index]
        receiver |= {sender}
        if len(receiver) > self.recipe_count_bounds[recipe][0]:
            if not self.assign_recipe_bound_min(recipe, len(receiver)):
                log.info(
                    "contradiction seen (recipe=%r, sender=%r)",
                    recipe, sender)
                return False
        return True

    def assign_part_usage(self, part, used):
        isValid = constraint_assign(self.parts_used[part], used)
        if isValid and used:
            isValid &= self.add_part_access_constraints(part)
        return isValid

    def add_part_access_constraints(self, part):
        isValid = True
        self.parts_accesses[part] = PART_ACCESSES[part]
        if min(self.parts_accesses[part]) > 0:
            # propagation logic
            isValid &= self.propagate_part_access_to_shapes(part)
            isValid &= self.propagate_min_access_to_track_count()
        return isValid

    def assign_part_accesses(self, part, count):
        assert part in self.parts_accesses
        old_part_min_access = min(self.parts_accesses[part])
        isValid = True
        # always call this
        isValid &= constraint_assign(self.parts_accesses[part], count)
        # if the min access increases, check for propagation
        if isValid and count > old_part_min_access:
            # propagation logic
            isValid &= self.propagate_part_access_to_shapes(part)
            isValid &= self.propagate_min_access_to_track_count()
        return isValid

    def eliminate_part_accesses(self, part, count):
        assert part in self.parts_accesses
        old_part_min_access = min(self.parts_accesses[part])
        isValid = True
        # always call this
        isValid &= constraint_remove(self.parts_accesses[part], count)
        # if the min access is removed, check for propagation
        if isValid and count == old_part_min_access:
            # propagation logic
            isValid &= self.propagate_part_access_to_shapes(part)
            isValid &= self.propagate_min_access_to_track_count()
        return isValid

    # arm length + track shape

    def eliminate_armlength_trackshape(self, armlength_trackshape):
        return constraint_remove(self.armlength_trackshapes,
            armlength_trackshape)

    # TODO: complete all assign and elimination methods

    # === methods for checkpointing branch points ===

    def branch_and_bound(self):
        # first check parts
        if not self.are_all_parts_constrained():
            # target parts for next search
            unconstrained_parts = [
                part for part, is_used in self.parts_used.items()
                if len(is_used) > 1
            ]
            unconstrained_glyphs = [
                part for part in unconstrained_parts
                if part in om.Part.PARTS_GLYPHS
            ]

        # if all parts are constrained, create a list of possible recipe metrics

        # TODO: finish branch and bound method
