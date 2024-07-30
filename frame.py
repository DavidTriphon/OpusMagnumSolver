import bisect
import math

import hexmath
import om

from sharedtypes import *


# programming for the spacial search of transformations to a molecule to put
# it in the right spot:

# node state format should possibly be entire grid state, including parts
# for now we just need a set of all molecules,
# where each molecule is just a set of atoms and bonds in positions
# and the molecules need to be testable for equality and overlap

# this could be done by modifying om.py or by implementing our own classes
# for objects

# positive rotations are CCW, negative rotations are CW.
# uses right hand system. redblob website uses left hand system.
# rotation = 0 is directly right
# negative and positive rotations are accepted values but are all (%= 6)'ed

# (q,r,s) is bottom right, top, bottom left.
# if q and r are x and y, then x is right, y is top right.


class Frame:

    def __init__(self,
            puzzle: om.Puzzle,
            atoms: list[om.Atom] = None,
            bonds: list[om.Bond] = None,
            parts: list[om.Part] = None):
        self.puzzle: om.Puzzle = puzzle
        self.atoms: list[om.Atom] = atoms or []
        self.bonds: list[om.Bond] = bonds or []
        self.parts: list[om.Part] = parts or []
        self.produced: list[int] = [0 for _ in puzzle.products]

    def copy(self):
        return Frame(
            self.puzzle,
            atoms=[atom.copy() for atom in self.atoms],
            bonds=[bond.copy() for bond in self.bonds],
            parts=[part.copy() if part.name in om.Part.PARTS_ARMS else part
                for part in self.parts]
        )

    # adders

    def add_molecule(self, molecule: om.Molecule, position: tuple = (0, 0),
            rotation: int = 0):
        # duplicate atoms
        new_atoms = [atom.copy() for atom in molecule.atoms]
        for atom in new_atoms:
            atom.translate(position, rotation)
        self.atoms.extend(new_atoms)
        # duplicate bonds
        new_bonds = [bond.copy() for bond in molecule.bonds]
        for bond in new_bonds:
            bond.translate(position, rotation)
        self.bonds.extend(new_bonds)

    def add_atom(self, atom: om.Atom):
        self.atoms.append(atom)

    def add_bond(self, bond: om.Bond):
        self.bonds.append(bond)

    def add_part(self, part: om.Part):
        self.parts.append(part.copy())

    # getters

    def get_atom(self, pos: Pos2D) -> om.Atom | None:
        return next((a for a in self.atoms if pos == a.position), None)

    def get_bond(self, pos0: Pos2D, pos1: Pos2D) -> om.Bond | None:
        pair = tuple(sorted((pos0, pos1)))
        return next((b for b in self.bonds if pair == b.positions), None)

    def get_arm_parts(self) -> list[om.Part]:
        return sorted(
            [part for part in self.parts
                if part.name in om.Part.PARTS_PROGRAMMABLE],
            key=lambda part: part.arm_number
        )

    def get_parts_in_set(self, name_set: set[bytes]) -> list[om.Part]:
        return [
            part for part in self.parts
            if part.name in name_set
        ]

    def all_grabbed_pos(self):
        return [
            grab_pos
            for arm in self.get_arm_parts()
            if arm.grabbing
            for grab_pos, grabbed in zip(arm.arm_grab_positions(), arm.grabbed)
            if grabbed
        ]

    def get_step4_parts(self) -> list[om.Part]:
        return [part for part in self.parts
            if part.name in om.Part.PARTS_POST_MOVEMENT]

    def get_step2_parts(self) -> list[om.Part]:
        return [part for part in self.parts
            if part.name in om.Part.PARTS_CONSUMERS]

    def get_molecules(self) -> list[om.Molecule]:
        linked_atoms = {
            atom.position: {atom}
            for atom in self.atoms
        }
        for bond in self.bonds:
            pos0 = bond.positions[0]
            pos1 = bond.positions[1]
            union = linked_atoms[pos0] | linked_atoms[pos1]
            for atom in union:
                linked_atoms[atom.position] = union
        groups = set(frozenset(group) for group in linked_atoms.values())
        molecules = [
            om.Molecule(
                atoms=list(group),
                bonds=[
                    bond for bond in self.bonds
                    if any(atom.position in bond.positions for atom in group)
                ]
            )
            for group in groups
        ]
        assert sum(len(molecule.bonds) for molecule in molecules) == len(
            self.bonds)
        assert sum(len(molecule.atoms) for molecule in molecules) == len(
            self.atoms)

        return molecules

    def get_simplified_rep(self):
        mol_count = len(self.get_molecules())

        atom_type_count_tuples = tuple(
            ("Atom:" + om.Atom.TYPE_NAMES[a_type], count)
            for (a_type, count) in self.count_atom_types().items()
        )
        bond_type_count_tuples = tuple(
            ("Bond:" + om.Bond.type_name(b_type), count)
            for (b_type, count) in self.count_bond_types().items()
        )
        return tuple(sorted(
            atom_type_count_tuples +
            bond_type_count_tuples +
            (("Molecules", mol_count),)
        ))

    def count_looping_bonds(self):
        return len(self.get_molecules()) - (len(self.atoms) - len(self.bonds))

    def count_atom_types(self):
        counts = {}
        for atom in self.atoms:
            counts[atom.type] = counts.get(atom.type, 0) + 1
        return counts

    def count_bond_types(self):
        counts = {}
        for bond in self.bonds:
            counts[bond.type] = counts.get(bond.type, 0) + 1
        return counts

    # evaluators

    def overlaps_track(self, pos: Pos2D):
        track_parts = self.get_parts_in_set({om.Part.TRACK})
        all_track_poses = {pos for part in track_parts
            for pos in part.track_hexes}
        return pos in all_track_poses

    def is_valid(self):
        # atoms cannot coexist in the same space
        atom_pos_set = set(atom.position for atom in self.atoms)
        if len(atom_pos_set) != len(self.atoms):
            return False
        # bonds must be attached to atoms
        bond_position_set = set(pos for bond in self.bonds for pos
            in bond.positions)
        if any(pos not in atom_pos_set for pos in bond_position_set):
            raise ValueError("Broken State!!!")
        # arms cannot exist in the same space as atoms
        arm_parts = self.get_arm_parts()
        arm_pos_set = {part.position for part in arm_parts}
        if any(arm_pos in atom_pos_set for arm_pos in arm_pos_set):
            return False
        # arms cannot coexist in the same space
        if len(arm_pos_set) != len(arm_parts):
            raise ValueError("Broken State!!!")
        # arm numbers should be sequential with no gaps
        if list(range(len(arm_parts))) != [arm.arm_number for arm in arm_parts]:
            raise ValueError("Broken State!!!")
        # if an arm has grabbed an atom, it should be in the spot of the hand
        for arm in arm_parts:
            if arm.grabbing and arm.grabbed:
                for grabbed_pos, grabbed \
                        in zip(arm.arm_grab_positions(), arm.grabbed):
                    if grabbed:
                        if grabbed_pos not in [
                            atom.position for atom in self.atoms
                        ]:
                            raise ValueError("Broken State!!!")

        # after everything else, it's valid
        return True

    # complicated and important methods

    def search(self, match_frame=None, satisfy_condition=None, heuristic=None):
        if match_frame is not None:
            def satisfy_condition(frame):
                return frame == match_frame
        if satisfy_condition is None:
            raise ValueError("match_frame or satisfy_condition must be set")

        h = heuristic or (lambda x: 0)

        frame_cost: dict[Frame, int] = {self: 0}
        frontier_queue: list[tuple[int, Frame]] = [(h(self), self)]
        explored: set[Frame] = set()
        explored_simplified: set[tuple] = set()
        came_from: dict[Frame, tuple[Frame, dict[int, bytes]]] = dict()

        def instruction_path(frame):
            instrs = []
            while frame in came_from:
                instrs.insert(0, came_from[frame][1])
                frame = came_from[frame][0]
            return instrs

        # fScore: dict[Frame, int] = {self: 0}
        past_count = 0
        # iterate until exit is found or queue is empty
        while frontier_queue:
            new_count = len(frame_cost) % 1000
            if past_count > new_count:
                print("Frames explored = %d, frontier size = %d, depth = %d" % (
                    len(frame_cost), len(frontier_queue), frame_cost[
                        frontier_queue[-1][1]]))
            past_count = new_count

            current_frame = frontier_queue.pop()[1]
            if current_frame not in explored:
                explored.add(current_frame)
                simp_rep = current_frame.get_simplified_rep()
                if simp_rep not in explored_simplified:
                    explored_simplified.add(simp_rep)
                # exit when the goal is found
                if satisfy_condition(current_frame):
                    return (instruction_path(current_frame),
                    current_frame, frame_cost[current_frame])
                for neighbor_frame, neighbot_instrs in \
                        current_frame.neighbor_frames():
                    neighbor_cost = frame_cost[current_frame] + 1
                    if neighbor_cost < frame_cost.get(neighbor_frame, math.inf):
                        frame_cost[neighbor_frame] = neighbor_cost
                        came_from[neighbor_frame] = (
                            current_frame, neighbot_instrs)
                        if neighbor_frame in explored:
                            explored.remove(neighbor_frame)
                        h_cost = neighbor_cost + h(neighbor_frame)
                        bisect.insort_left(frontier_queue,
                            (h_cost, neighbor_frame), key=lambda t: -t[0])

        # frontier is empty and goal was never found
        return None, None, None

    def neighbor_frames(self):
        arm_parts = self.get_arm_parts()
        arm_instructions = [arm.valid_instructions(
            overlaps_track=self.overlaps_track(arm.position))
            for arm in arm_parts]

        def combos(list_of_lists_of_instrs):
            if len(list_of_lists_of_instrs) == 1:
                return [[x] for x in list_of_lists_of_instrs[0]]
            mid_i = len(list_of_lists_of_instrs) // 2
            a_s = combos(list_of_lists_of_instrs[:mid_i])
            b_s = combos(list_of_lists_of_instrs[mid_i:])
            return [a + b for a in a_s for b in b_s]

        instruction_combos = [{i: instr for i, instr in enumerate(combo) if
            instr is not None}
            for combo in combos(arm_instructions)]

        neighbors = []
        for combo in instruction_combos:
            copy = self.copy()
            if copy.iterate(combo):
                neighbors.append((copy, combo))

        return neighbors

    def iterate(self, arm_instructions: dict[int, bytes] = None):
        arm_instructions = arm_instructions or {}
        # instructions must index to existing arms
        arm_parts = self.get_arm_parts()
        if any(key >= len(arm_parts) for key in arm_instructions.keys()):
            raise ValueError(
                "arm instruction keys must be within the number "
                "of arms. len(arms)=%d, max(key)=%d" % (
                    len(arm_parts), max(arm_instructions.keys()))
            )

        # identify what atoms are linked to each other and what bonds are
        # involved
        molecules = self.get_molecules()

        # identify all rotation and slide movements of molecules
        queue_mol_rotate: set[tuple[om.Molecule, Pos2D, Angle]] = set()
        queue_arm_rotate: list[tuple[om.Part, Angle]] = []

        # slide_queue = set()

        def queue_rotate(arm, theta, pivoting=False):
            if arm.grabbing and arm.grabbed:
                for grabbed_pos, grabbed in zip(arm.arm_grab_positions(),
                        arm.grabbed):
                    if grabbed:
                        grabbed_atom = self.get_atom(grabbed_pos)
                        queue_mol_rotate.add((
                            next(m for m in molecules
                                if grabbed_atom in m.atoms
                            ),
                            (grabbed_atom.position
                             if pivoting else
                             arm.position),
                            theta
                        ))
            # update the arm state afterwards
            if not pivoting:
                queue_arm_rotate.append((arm, theta))

        for arm in arm_parts:
            # nullable instruction, key not guaranteed, use get()
            instruction = arm_instructions.get(arm.arm_number)

            if instruction == om.Instruction.ROTATE_CCW:
                queue_rotate(arm, 1)
            elif instruction == om.Instruction.ROTATE_CW:
                queue_rotate(arm, -1)
            elif instruction == om.Instruction.PIVOT_CCW:
                queue_rotate(arm, 1, True)
            elif instruction == om.Instruction.PIVOT_CW:
                queue_rotate(arm, -1, True)
            elif instruction == om.Instruction.GRAB:
                # arms cannot grab if they are already grabbing
                if not arm.grabbing:
                    grab_positions = arm.arm_grab_positions()
                    arm.grabbing = True
                    arm.grabbed = [self.get_atom(grab_pos)
                        for grab_pos in grab_positions]
                queue_rotate(arm, 0)
            elif instruction == om.Instruction.DROP:
                arm.grabbing = False
                arm.grabbed = None
            elif instruction is None:
                if arm.grabbing and arm.grabbed:
                    queue_rotate(arm, 0, False)
            else:
                raise NotImplementedError("Instruction not implemented for "
                                          "iterate_frame: (%r)" % instruction)

        # TODO: go through the rotate queue and check for molecules grabbed
        #  by different pivot points that could crash the simulation

        # TODO: check if any of the rotations would collide with each other
        #  mid-swing

        # TODO: add pre-movement interaction checks
        # handle standard outputs
        outputs = self.get_parts_in_set({om.Part.OUTPUT_STANDARD})
        for out in outputs:
            out_mol = self.puzzle.products[out.which_reagent_or_product]
            moved_mol = out_mol.copy()
            moved_mol.translate(out.position, out.rotation)
            matching_mols = [m for m in self.get_molecules() if
                m == moved_mol]
            if matching_mols:
                matched = matching_mols[0]
                grabbed_poses = self.all_grabbed_pos()
                if not any(atom.position in grabbed_poses
                        for atom in matched.atoms):
                    self.produced[out.which_reagent_or_product] += 1
                    for atom in matched.atoms:
                        self.atoms.remove(atom)
                    for bond in matched.bonds:
                        self.bonds.remove(bond)

        # update the positions of the linked atoms and bonds
        for grab_mol, pivot_pos, angle in queue_mol_rotate:
            grab_mol.rotate(pivot_pos, angle)
        for arm, theta in queue_arm_rotate:
            arm.rotation = (arm.rotation + theta) % 6

        # TODO: complete post-movement interaction checklist
        # handle bonders
        bonders = self.get_parts_in_set({om.Part.BONDER})
        for bonder in bonders:
            pos0 = bonder.position
            pos1 = hexmath.summate(bonder.position, hexmath.ROTATION_VECTORS[
                bonder.rotation])

            atom0 = self.get_atom(pos0)
            atom1 = self.get_atom(pos1)
            bond_opt = self.get_bond(pos0, pos1)

            if atom0 and atom1 and not bond_opt:
                self.bonds.append(om.Bond(om.Bond.NORMAL, (pos0, pos1)))

        # handle inputs
        inputs = self.get_parts_in_set({om.Part.INPUT})
        for inp in inputs:
            inp_mol = self.puzzle.reagents[inp.which_reagent_or_product]
            inp_positions = {
                hexmath.translate(a.position, inp.position, inp.rotation)
                for a in inp_mol.atoms
            }
            covering_atoms = [a for a in self.atoms if
                a.position in inp_positions]
            if not covering_atoms:
                self.add_molecule(inp_mol, inp.position, inp.rotation)

        # handle calcification
        calcifiers = self.get_parts_in_set({om.Part.CALCIFICATION})
        for calcifier in calcifiers:
            overlapping_atom = self.get_atom(calcifier.position)
            if (overlapping_atom and overlapping_atom.type in
                    om.Atom.TYPES_ELEMENTAL):
                overlapping_atom.type = om.Atom.SALT

        # maybe perform a validity check? might as well if we're checking
        #  iteration validity already
        if not self.is_valid():
            return False

        return self

    # default methods

    def __eq__(self, other):
        if not isinstance(other, Frame):
            return False
        return (
                self.puzzle is other.puzzle and
                all(a == b for a, b in zip(self.atoms, other.atoms)) and
                all(a == b for a, b in zip(self.bonds, other.bonds)) and
                all(a == b for a, b in zip(self.parts, other.parts)) and
                self.produced == other.produced
        )

    def __hash__(self):
        return hash((
            self.puzzle, tuple(self.atoms), tuple(self.bonds),
            tuple(self.parts), tuple(self.produced)
        ))
