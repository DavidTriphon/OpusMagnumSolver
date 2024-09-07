import cProfile
import math
import pstats
import traceback
from pathlib import Path
from timeit import default_timer as timer

import hexmath
from frame import Frame
import puzzledb
import display
import om


### TEST API FOR OM ###


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


def main_investigate_puzzle():
    puzzle = om.Puzzle("resources\\24hour-1-sample\\GEN000.puzzle")
    display.print_puzzle_details(puzzle)
    solution = not_a_solution(puzzle)
    solution.write_to_path("output\\24hour-1-sample\\GEN000-test.solution")


### UNUSED FUNCTIONS ###


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


def find_atom_steps(puzzle: om.Puzzle):
    # assumes 1 reagent with 1 atom and 1 product for which every atom type
    # has only one manufacture path
    product = puzzle.products[0]
    reagent = puzzle.reagents[0]
    pr_mapping = {i: 0 for i in range(len(puzzle.products[0].atoms))}
    return (
            {
                ("calcify", (p_i,))
                for p_i, p_atom in enumerate(puzzle.products[0].atoms)
                # assumes the only path is identity or calcification
                if p_atom.type != reagent.atoms[pr_mapping[p_i]].type
            } | {
                # assume all bonds come from bonder
                ("bond", (
                    next((i for i, a in enumerate(product.atoms)
                        if bond.positions[0] == a.position), None),
                    next((i for i, a in enumerate(product.atoms)
                        if bond.positions[1] == a.position), None)
                )) for bond in product.bonds
            }
    )


def generator_permutations(values, prefix=None):
    if prefix is None:
        prefix = []
    if len(values) == 0:
        yield tuple(prefix)
    for val in values:
        for x in generator_permutations(values - {val}, prefix + [val]):
            yield x


def find_calcifier_position_constraints(product: om.Molecule,
        armlength: int = 1):
    grab_pos = (armlength, 0)
    pos_dict_i_dict_req_list = {}
    for atom_grab_i in range(len(product.atoms)):
        for rot in range(6):
            rotated_product = product.copy().translate(
                hexmath.difference(grab_pos,
                    product.atoms[atom_grab_i].position)
            ).rotate(grab_pos, rot)
            collided_atoms_is = [i for i, atom in
                enumerate(rotated_product.atoms)
                if atom.position == (0, 0)]
            atom_hextants = [hexmath.hextant(atom.position)
                for atom in rotated_product.atoms]
            modified_positions = [
                hexmath.rotate(atom.position, -hextant)
                if hextant is not None else None
                for hextant, atom in zip(atom_hextants, rotated_product.atoms)
            ]

            for i, pos in enumerate(modified_positions):
                if pos is not None and product.atoms[i].type == om.Atom.SALT:
                    i_dict_req_list = pos_dict_i_dict_req_list.setdefault(
                        pos, {})
                    req_list = i_dict_req_list.setdefault(i, [])

                    required_atoms = (list(range(atom_grab_i, i + 1))
                                      if atom_grab_i <= i else
                                      list(range(i, atom_grab_i + 1)))
                    req_list.append({
                        "req": required_atoms,
                        "req_not": collided_atoms_is,
                        "arm_rot": (-atom_hextants[i]) % 6,
                        "molecule_rel_pivot": rot,
                        "grabbed": atom_grab_i
                    })

    return pos_dict_i_dict_req_list


def find_bonder_position_constraints():
    return [
        {
            "position": hexmath.pos_from_direction(1, dir),
            "rotation": (dir + 2) % 6,
            "occupies": {hexmath.pos_from_direction(1, dir),
                hexmath.pos_from_direction(1, dir + 1)}
        }
        for dir in range(6)
    ]


def find_reagent_positions():
    return [
        hexmath.pos_from_direction(1, dir)
        for dir in range(6)
    ]


def find_solution(puzzle: om.Puzzle):
    reagent_positions = find_reagent_positions()
    reagent_position = reagent_positions[0]
    bonder_constraints = find_bonder_position_constraints()
    calcifier_constraints = find_calcifier_position_constraints(
        puzzle.products[0])
    salt_count = len([atom for atom in puzzle.products[0].atoms
        if atom.type == om.Atom.SALT])
    calcifier_positions = [
        hexmath.rotate(calc_pos, rot)
        for calc_pos, constraint in calcifier_constraints.items()
        if len(constraint) == salt_count
        for rot in range(6)
    ]

    arm_part = om.Part(name=om.Part.ARM1, position=(0, 0),
        rotation=hexmath.direction_int(reagent_position),
        length=1, arm_number=0)
    reagent_part = om.Part(name=om.Part.INPUT,
        position=reagent_position, rotation=0)
    occupied = {arm_part.position, reagent_part.position}

    for bonder_details in bonder_constraints:
        if all(o not in occupied for o in bonder_details["occupies"]):
            occupied1 = occupied | bonder_details["occupies"]
            bonder_part = om.Part(name=om.Part.BONDER,
                position=bonder_details["position"],
                rotation=bonder_details["rotation"])
            for calcifier_pos in calcifier_positions:
                if calcifier_pos not in occupied1:
                    # occupied2 = occupied1 | {calcifier_pos}
                    calcifier_part = om.Part(name=om.Part.CALCIFICATION,
                        position=calcifier_pos)
                    parts = [arm_part, reagent_part, bonder_part,
                        calcifier_part]
                    start_frame = Frame(puzzle, parts=parts).iterate({})
                    instructions, path_frames, cost = find_391_program(
                        start_frame)
                    if instructions:
                        return create_solution(puzzle, "cost attempt",
                            instructions, parts)


### SAVE AND RECORD ALL GENERATED SOLUTIONS ###


def omsim_record_save(puzzle, solution):
    solve_sim = om.Sim(puzzle=puzzle, solution=solution)

    try:
        sim_cost = solve_sim.metric("cost")
        sim_cycles = solve_sim.metric("cycles")
        sim_area = solve_sim.metric("area")
        print("solution validated: cost=%d, cycles=%d, area=%d" % (
            sim_cost, sim_cycles, sim_area))
        solution_filename = "%s_%dg_%dc_%da.solution" % (
            puzzle.name.decode("utf-8"), sim_cost, sim_cycles, sim_area)
        solution.write_to_path(puzzledb.OUTPUT_PATH / solution_filename)
    except om.SimError as e:
        print("solution was found invalid by omsim:")
        print(e)


### SOLVE FOR 391 ###


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


def rotate_diff(x: int, y: int) -> int:
    return abs((y - x + 3) % 6 - 3)


def dist_from_new_atom(frame: Frame, e: Evaluator = None) -> int:
    e = e or Evaluator(frame)
    if e.atoms() >= 5:
        return 0

    # making an atom takes 3 cycles each (4- to offset wip atom)
    cost = 3 * max(0, 4 - e.atoms())

    arm = frame.get_arm_parts()[0]
    distance = rotate_diff(0, arm.rotation)
    # the distance is always a part of the cost for new atoms
    # but the distance can overlap with other tasks, like calc or bond
    # reduce addition to the cost by 1 for each of those unfinished
    # categories.
    # only care about reduction if the arm is grabbing
    """
    if arm.grabbing and arm.grabbed[0]:
        reduction = (e.atom_count_type(om.Atom.SALT < 3)) + (e.bonds() < 3)
        cost += max(0, distance - reduction)
    else:
        cost += distance
    """

    cost += (distance > 0)

    # whether the arm is grabbing or not determines whether a drop is
    # necessary or if a grab is necessary
    if distance == 0:
        if arm.grabbing:
            # ensure it's actually grabbed
            # grabbing an atom on top of the input means we have it
            # not grabbed means drop and then grab
            if not arm.grabbed[0]:
                cost += 2
        else:
            # still need to grab
            cost += 1
    else:
        if arm.grabbing:
            # grabbing a different atom means drop, move dist, grab
            cost += 2
        else:
            # not grabbing means move dist, grab
            cost += 1

    return cost


def bond_cost(frame: Frame = None, e: Evaluator = None) -> int:
    e = e or Evaluator(frame)
    cost = max(0, 3 - e.bonds())
    return max(cost, 0)


def salt_cost(frame: Frame = None, e: Evaluator = None) -> int:
    e = e or Evaluator(frame)
    cost = max(0, 3 - e.atom_count_type(om.Atom.SALT))
    return cost


def output_dist(frame: Frame, e: Evaluator = None) -> int:
    e = e or Evaluator(frame)
    if (e.atoms() == 5 and e.atom_count_type(om.Atom.SALT) == 3 and
            e.bonds() == 3):
        arm = frame.get_arm_parts()[0]
        cost = 1  # for the final drop
        if arm.rotation != 3:
            cost += 1
        return cost
    else:
        return 2


def bond_punish_cost(frame: Frame = None, e: Evaluator = None) -> int:
    e = e or Evaluator(frame)
    return (  # punish looping bonds
            (math.inf if e.atoms() - e.bonds() < 1 else 0)
            # punish extra bonds
            + (math.inf if e.bonds() > 3 else 0)
    )


def punish_bad_angles(frame: Frame, e: Evaluator = None) -> int | float:
    # no possibility of punishment if less than 2 bonds
    if e.bonds() < 2:
        return 0

    pos_nbs = dict()
    for bond in frame.bonds:
        pos1, pos2 = bond.positions
        pos_nbs[pos1] = pos_nbs.get(pos1, set()) | {pos2}
        pos_nbs[pos2] = pos_nbs.get(pos2, set()) | {pos1}

    for pos, nbs in pos_nbs.items():
        # ensure no atom has more than 2 bonds
        if len(nbs) > 2:
            return math.inf
        # ensure all bonds that share a pos are at angles of 2
        if len(nbs) == 2:
            pos1, pos2 = nbs
            pos1 = hexmath.difference(pos1, pos)
            pos2 = hexmath.difference(pos2, pos)
            dir1 = hexmath.direction_int(pos1)
            dir2 = hexmath.direction_int(pos2)
            if not hexmath.angle_directions(dir1, dir2) == 2:
                return math.inf

    return 0


def overfitted_heuristic(frame: Frame) -> int:
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
            + bond_punish_cost(e=e)
    )


def output_heuristic(frame: Frame) -> int:
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
            + bond_punish_cost(e=e)
    )


def consistent_heuristic(frame: Frame) -> int:
    if frame.produced[0] > 0:
        return 0
    e = Evaluator(frame)
    return (
            max(salt_cost(e=e), bond_cost(frame, e))
            + dist_from_new_atom(frame, e)
            + bond_punish_cost(e=e)
            + output_dist(frame, e)
            + punish_bad_angles(frame, e)
    )


def attach_instructions_to_solution(solution,
        instructions_dict_list):
    # reformat instructions from search into solution format
    # add instructions to arm
    for arm in solution.get_arm_parts():
        arm.instructions = [
            om.Instruction(i, frame_instrs[arm.arm_number])
            if arm.arm_number in frame_instrs else
            om.Instruction(i, om.Instruction.NOOP)
            for i, frame_instrs in enumerate(instructions_dict_list)
        ]


def find_391_program(start_frame):
    # instruction path to product output
    output_instructions, output_frames, output_cost = start_frame.search(
        goal_condition=lambda f: f.produced[0] >= 1,
        heuristic=consistent_heuristic
    )
    # instruction path to reset
    return_instructions, return_frames, return_cost = output_frames[-1].search(
        goal_condition=lambda f: f.parts == start_frame.parts
    )
    return (output_instructions + return_instructions,
    output_frames + return_frames, output_cost + return_cost)


def create_solution(puzzle, solution_name, instructions, partlist):
    # create and save solution
    solution = om.Solution(
        puzzle=puzzle.name,
        name=bytes(solution_name, "utf-8"),
        parts=partlist
    )
    # reformat instructions from search into solution format
    attach_instructions_to_solution(solution, instructions)
    return solution


def save_solution(puzzle, solution):
    # save solution
    solution_filename = "%s_%s.solution" % (
        puzzle.name.decode("utf-8"),
        solution.name.decode("utf-8").replace(" ", "_"))
    solution.write_to_path(puzzledb.SAVEDATA_PATH / solution_filename)
    omsim_record_save(puzzle, solution)


def main_create_cost_solve_391():
    # select puzzle
    puzzles = puzzledb.get_test_puzzles()
    puzzle = puzzles[391]

    # create initial start frame of the simulation
    start_frame = Frame(puzzle=puzzle, parts=create_391_solve_partlist())
    start_frame.iterate({})

    time_start = timer()
    instructions, output_frames, cost = find_391_program(start_frame)
    time_end = timer()

    print("search duration=%.2f seconds" % (time_end - time_start))
    print("cost = %d" % cost)
    print()
    print("heuristics:")
    for heuristic in [overfitted_heuristic, output_heuristic,
        consistent_heuristic]:
        print("  - %s:" % heuristic.__name__)
        print("    %r" % [heuristic(frame) + i
            for i, frame in enumerate([start_frame] + output_frames)])

    solution = create_solution(puzzle, "cost prototype", instructions,
        start_frame.parts)
    save_solution(puzzle, solution)


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
