import cProfile
import math
import pstats
import traceback
from pathlib import Path
from timeit import default_timer as timer

import analysis
from frame import Frame
import puzzledb
import recipes
import puzzleparts
import display
import om


def not_a_solution(puzzle: om.Puzzle):
    next_pos = (0, 0)
    parts = []
    for i, reagent in enumerate(puzzle.reagents):
        xs = [atom.position[0] for atom in reagent.atoms]
        ys = [atom.position[1] for atom in reagent.atoms]
        lower = min(xs), min(ys)
        upper = max(xs), max(ys)
        size = tuple(u - l for u, l in zip(upper, lower))
        parts.append(om.Part(name=om.Part.INPUT, which_reagent_or_product=i,
            position=next_pos))
        next_pos = (next_pos[0] + size[0] + 2, next_pos[1] + size[1] + 2)
    for i, product in enumerate(puzzle.products):
        xs = [atom.position[0] for atom in product.atoms]
        ys = [atom.position[1] for atom in product.atoms]
        lower = min(xs), min(ys)
        upper = max(xs), max(ys)
        size = tuple(u - l for u, l in zip(upper, lower))
        parts.append(
            om.Part(name=om.Part.OUTPUT_STANDARD, which_reagent_or_product=i,
                position=next_pos))
        next_pos = (next_pos[0] + size[0] + 2, next_pos[1] + size[1] + 2)

    sol = om.Solution(parts=parts)
    return sol


def validate_workshop_puzzle(puzzle_id, solution_digit):
    puzzle_path = Path("savedata\\workshop\\{puzzle_name}.puzzle")
    if not puzzle_path.exists():
        raise FileNotFoundError(f"{puzzle_path} does not exist")
    puzzle = om.Puzzle(str(puzzle_path))
    solution_path = Path(f"savedata\\{puzzle.name}-{puzzle_id}"
                         f"-{solution_digit}.puzzle")
    if not solution_path.exists():
        raise FileNotFoundError(f"{solution_path} does not exist")
    solution = om.Solution(str(solution_path))
    om.Sim(puzzle, solution)


def possible_product_types(puzzle: om.Puzzle) -> list[int]:
    reagent_types = puzzle.reagent_types()
    product_types = reagent_types.copy()
    hovering_types = []
    type_checklist = product_types.copy()

    available_parts = puzzleparts.full_parts_list(puzzle)
    available_recipes = [r for r in recipes.RECIPES
        if r.part in available_parts]

    if om.Part.BERLO in available_parts:
        hovering_types = analysis.SIMPLE_ELEMENTAL_ATOM_TYPES

    while type_checklist:
        next_type = type_checklist.pop()
        for recipe in available_recipes:
            if next_type in recipe.ingredients:
                if all([ing in product_types for ing in recipe.ingredients]):
                    if all([
                        hov in hovering_types or hov in product_types
                        for hov in recipe.hovering
                    ]):
                        for product in recipe.products:
                            if product not in product_types:
                                product_types.append(product)
                                type_checklist.append(product)

    return product_types


def calc_io_ratio(puzzle: om.Puzzle) -> tuple[int, int]:
    input_types = puzzle.reagent_types()
    output_types = puzzle.product_types()

    # calculate what output types can be made from what input types
    # assume for elemental ones that ratio is 1:1 for now
    # use puzzle parts available flags to determine what parts are available
    # and from there what recipes are available to apply on the atom
    # types

    # once 1-step map is made of input to output types, figure out how many
    # inputs have to be repeated in order to construct 1 output.

    # later for metal puzzles it will be necessary to determine the cycles
    # ratio using an equation using the Power Log function. This will remove
    # linear time from the method in calculating the most efficient ratio in
    # puzzles that have both projection and purify.

    # the first match when considering satisfying the output should be
    # material as is, and then figuring out if alternative material is viable
    # to use. all recipes are reversible without loss, or not
    # reversible at all, so hopefully a greedy approach won't introduce too
    # many problems.

    # if remainders after 1 output, determine if the remainder is desirable to
    # use in the next output set or if it is to be discarded
    # if not to be discarded, continue trying to construct outputs to
    # determine if output ratio is not 1:n or n:1. this will inform rate for
    # this puzzle and is to be returned as the answer.
    raise NotImplementedError


def main_investigate_puzzle():
    puzzle = om.Puzzle("resources\\24hour-1-sample\\GEN000.puzzle")
    display.print_puzzle_details(puzzle)
    solution = not_a_solution(puzzle)
    solution.write_to_path("output\\24hour-1-sample\\GEN000-test.solution")


def omsim_record_save(puzzle, solution):
    solve_sim = om.Sim(puzzle=puzzle, solution=solution)

    try:
        sim_cost = solve_sim.metric("cost")
        sim_cycles = solve_sim.metric("cycles")
        sim_area = solve_sim.metric("area")
        solution_filename = "%s_%dg_%dc_%da.solution" % (
            puzzle.name.decode("utf-8"), sim_cost, sim_cycles, sim_area)
        solution.write_to_path(puzzledb.OUTPUT_PATH / solution_filename)
    except om.SimError as e:
        print("solution was found invalid by omsim:")
        print(e)


class Evaluator:

    def __init__(self, frame: Frame):
        self.frame = frame
        self._atoms = None
        self._bonds = None
        self._atoms_type = {}

    def atoms(self):
        if self._atoms is None:
            self._atoms = len(self.frame.atoms)
        return self._atoms

    def bonds(self):
        if self._bonds is None:
            self._bonds = len(self.frame.bonds)
        return self._bonds

    def atom_count_type(self, type):
        if type not in self._atoms_type:
            self._atoms_type[type] = len(
                [a for a in self.frame.atoms if a.type == type])
        return self._atoms_type[type]


def create_391_solve_partlist():
    arm_part = om.Part(name=om.Part.ARM1, position=(0, 0), rotation=0,
        length=1, arm_number=0)
    reagent_part = om.Part(name=om.Part.INPUT, position=(1, 0), rotation=0)
    product_part = om.Part(name=om.Part.OUTPUT_STANDARD, position=(-3, 0),
        rotation=-1)
    bonder_part = om.Part(name=om.Part.BONDER, position=(-1, 1), rotation=0)
    calcifier_part = om.Part(name=om.Part.CALCIFICATION, position=(-2, 1))

    return [arm_part, reagent_part, product_part, bonder_part,
        calcifier_part]


def main_create_cost_solve_391():
    # select puzzle
    puzzles = puzzledb.get_test_puzzles()
    puzzle = puzzles[391]

    # create solution parameters
    solution_name = "cost_prototype"

    #   create solution object
    solution = om.Solution(
        puzzle=puzzle.name,
        name=bytes(solution_name.replace("_", " "), "utf-8"),
        parts=create_391_solve_partlist()
    )

    # create initial start frame of the simulation
    start_frame = Frame(puzzle=puzzle, parts=create_391_solve_partlist())
    start_frame.iterate({})

    def overfitted_heuristic(frame: Frame):
        if frame.produced[0] > 0:
            return 0
        e = Evaluator(frame)
        return (  # 3 bonds take 1 cycle each
                max(0, 3 - e.bonds())
                # making a water takes 4 cycles each
                + 4 * max(0, 2 - e.atom_count_type(om.Atom.WATER))
                # making a salt takes 1+4 cycles each punish looping bonds
                + 5 * max(0, 3 - e.atom_count_type(om.Atom.SALT))
                # punish looping bonds
                + (math.inf if e.atoms() - e.bonds() <= 1
                   else 0)
                # punish extra bonds
                + (math.inf if e.bonds() > 3 else 0)
        )

    def output_heuristic(frame: Frame):
        if frame.produced[0] > 0:
            return 0
        e = Evaluator(frame)
        return (  # 3 bonds take 1 cycle each
                max(0, 3 - e.bonds())
                # making an atom takes 4 cycles each
                + 4 * max(0, 5 - e.atoms())
                # making a salt takes 1 cycle each
                + max(0, 3 - e.atom_count_type(om.Atom.SALT))
                # punish looping bonds
                + (math.inf if e.atoms() - e.bonds() <= 1 else 0)
                # punish extra bonds
                + (math.inf if e.bonds() > 3 else 0)
        )

    time_start = timer()
    # instruction path to product output
    product_instructions, output_frame, output_cost = start_frame.search(
        satisfy_condition=lambda f: f.produced[0] >= 1,
        heuristic=output_heuristic
    )
    time_midpoint = timer()
    # instruction path to reset
    return_instructions, end_frame, return_cost = output_frame.search(
        satisfy_condition=lambda f: f.parts == start_frame.parts
    )
    time_end = timer()

    print("output search duration=%.2f seconds" % (time_midpoint - time_start))
    print("return search duration=%.2f seconds" % (time_end - time_midpoint))
    print("cost = %d" % (output_cost + return_cost))

    # reformat instructions from search into solution format
    arm_instrs = {
        arm.arm_number: [
            om.Instruction(i, frame_instrs[arm.arm_number])
            for i, frame_instrs in
            enumerate(product_instructions + return_instructions)
            if arm.arm_number in frame_instrs
        ]
        for arm in start_frame.get_arm_parts()
    }

    # add instructions to arm
    solution_arm = solution.parts[0]
    solution_arm.instructions = arm_instrs[0]

    # save solution
    solution_filename = "%s_%s.solution" % (
        puzzle.name.decode("utf-8"), solution_name)
    solution.write_to_path(puzzledb.SAVEDATA_PATH / solution_filename)
    omsim_record_save(puzzle, solution)


def profile_code(func, name=None):
    if name is None:
        name = func.__name__
    # > snakeviz ./needs_profiling.prof
    with cProfile.Profile() as pr:
        try:
            func()
        except Exception as e:
            traceback.print_exception(e)
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.dump_stats(filename="profiles\\" + name + '.prof')


if __name__ == '__main__':
    main_create_cost_solve_391()
    # profile_code(main_create_cost_solve_391)
