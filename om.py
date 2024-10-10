# CHANGELOG
#
#  2024-06-06 (panic): properly convert from memoryview to bytes when decoding strings
#  2024-05-06 (panic): add constants for part names
#  2024-05-06 (panic): properly use the `instructions` parameter in om.Solution.__init__() (thanks GuiltyBystander for the report)
#  2024-05-06 (panic): fix error in om.Instruction.REPEAT and om.Instruction.RESET (thanks gtw for the report)
#  2024-05-01 (panic): update om.Sim error message to link to the omsim repository
#  2024-04-29 (panic): fix parts_available constants
#  2024-04-24 (panic): move lcm import into the function that uses it for compatibility
#  2024-04-24 (panic): remove use of dataclasses for compatibility
#  2024-04-24 (panic): initial release

# EXAMPLE CODE
#
# === add an equilibrium glyph to a solution ===
#
#  import om
#  sol = om.Solution('/path/to/file.solution')
#  sol.parts.append(om.Part(name=om.Part.EQUILIBRIUM, position=(5, -2)))
#  sol.write_to_path('/path/to/file.solution')

# === change every fire to water in a puzzle's reagents ===
#
#  import om
#  puzzle = om.Puzzle('/path/to/file.puzzle')
#  for reagent in puzzle.reagents:
#      for atom in reagent.atoms:
#          if atom.type == om.Atom.FIRE:
#              atom.type = om.Atom.WATER
#  puzzle.write_to_path('/path/to/file.puzzle')

# === see if a solution completes with any of its parts removed ===
#
#  import om
#  sol = om.Solution('/path/to/file.solution')
#  for i in range(len(sol.parts)):
#      modified = om.Solution(sol)
#      del modified.parts[i]
#      try:
#          cycles = om.Sim('/path/to/file.puzzle', modified).metric('cycles')
#          print(f'solution completed in {cycles} cycles with part {i} removed')
#      except om.SimError:
#          pass

# API REFERENCE
#
# === om.Puzzle ===
#  om.Puzzle represents a puzzle file
#
#  puzzle = om.Puzzle()
#     creates a new, empty puzzle
#
#  puzzle = om.Puzzle('/path/to/file.puzzle')
#     loads a puzzle from a path
#
#  puzzle = om.Puzzle(b'...')
#     creates a puzzle from bytes in puzzle file format
#
#  data = puzzle.to_bytes()
#     encodes a puzzle as a byte array in puzzle file format and returns it
#
#  puzzle.write_to_path('/path/to/file.puzzle')
#     writes a puzzle to a path
#
#  puzzle_copy = om.Puzzle(puzzle)
#     creates a "deep copy" of a puzzle -- equivalent to om.Puzzle(puzzle.to_bytes())
#
#  puzzle = om.Puzzle(name=b'SAMPLE', reagents=[om.Molecule(...)])
#     creates a new, empty puzzle with some attributes set
#
#  puzzle = om.Puzzle()
#  puzzle.name = b'SAMPLE'
#  puzzle.reagents=[om.Molecule(...)]
#     creates a new, empty puzzle and sets some attributes afterward
#
#     here is the full list of attributes on om.Puzzle with some example values:
#  puzzle = om.Puzzle()
#  puzzle.name = b'puzzle name'
#  puzzle.creator = 0  (steam creator ID)
#  puzzle.parts_available = om.Puzzle.DEFAULT_PARTS_AVAILABLE | om.Puzzle.ANIMISMUS  (bitmask of which parts are enabled in the editor)
#  puzzle.reagents = [om.Molecule(...), ...]
#  puzzle.products = [om.Molecule(...), ...]
#  puzzle.output_scale = 1  (a scale of 1 means 6 outputs are required for completion)
#  puzzle.production_info = om.ProductionInfo(...)  (or None, if the puzzle is not a production puzzle)
#
#  true_or_false = puzzle.are_parts_available(om.Puzzle.QUINTESSENCE | om.Puzzle.DISPOSAL)
#     checks if all the given parts are available for use in solutions to the puzzle
#     different types of parts can be combined using the | operator
#
#  here is the full list of parts_available constants:
#     om.Puzzle.ARM
#     om.Puzzle.MULTIARM
#     om.Puzzle.PISTON
#     om.Puzzle.TRACK
#     om.Puzzle.BONDER
#     om.Puzzle.UNBONDER
#     om.Puzzle.MULTIBONDER
#     om.Puzzle.TRIPLEX
#     om.Puzzle.CALCIFICATION
#     om.Puzzle.DUPLICATION
#     om.Puzzle.PROJECTION
#     om.Puzzle.PURIFICATION
#     om.Puzzle.ANIMISMUS
#     om.Puzzle.DISPOSAL
#     om.Puzzle.QUINTESSENCE
#     om.Puzzle.GRAB_AND_ROTATE
#     om.Puzzle.DROP
#     om.Puzzle.RESET
#     om.Puzzle.REPEAT
#     om.Puzzle.PIVOT
#     om.Puzzle.BERLO
#     om.Puzzle.DEFAULT_PARTS_AVAILABLE

# === om.Molecule ===
#  om.Molecule is a collection of atoms and bonds representing a reagent or product
#
#  m = om.Molecule()
#     creates a new, empty molecule
#
#  m = om.Molecule(atoms=[om.Atom(...)], bonds=[om.Bond(...)])
#     creates a new molecule with the given atoms and bonds
#
#  m = om.Molecule()
#  m.atoms = [om.Atom(...)]
#  m.bonds = [om.Bond(...)]
#     creates a new molecule and sets the atoms and bonds afterward
#     these are the only two attributes on om.Molecule

# === om.Atom ===
#  om.Atom represents a single atom in a molecule
#
#  a = om.Atom(om.Atom.WATER, (0, 1))
#     creates a water atom located at position (0, 1) within a molecule
#
#  print(a.type, a.position)
#     print the numeric type and position (e.g., water has numeric code 5)
#
#  here is the full list of atom types:
#     om.Atom.SALT
#     om.Atom.AIR
#     om.Atom.EARTH
#     om.Atom.FIRE
#     om.Atom.WATER
#     om.Atom.QUICKSILVER
#     om.Atom.GOLD
#     om.Atom.SILVER
#     om.Atom.COPPER
#     om.Atom.IRON
#     om.Atom.TIN
#     om.Atom.LEAD
#     om.Atom.VITAE
#     om.Atom.MORS
#     om.Atom.REPETITION_PLACEHOLDER
#     om.Atom.QUINTESSENCE

# === om.Bond ===
#  om.Bond represents a bond between two atoms in a molecule
#
#  b = om.Bond(om.Bond.NORMAL, ((0, 0), (0, 1)))
#     creates a bond between positions (0, 0) and (0, 1) of a molecule
#
#  print(b.type, b.positions)
#     print the numeric bond type and the positions tuple (e.g., a normal bond has numeric code 1)
#
#  here is the full list of bond types:
#     om.Bond.NORMAL
#     om.Bond.TRIPLEX_RED
#     om.Bond.TRIPLEX_BLACK
#     om.Bond.TRIPLEX_YELLOW
#     om.Bond.TRIPLEX  (all three types of triplex combined)
#  bond types can be combined using the | operator, e.g. om.Bond.TRIPLEX_RED | om.Bond.TRIPLEX_BLACK

# === om.ProductionInfo ===
#  om.ProductionInfo stores information about a production puzzle
#
#  pi = om.ProductionInfo()
#     creates a new, empty ProductionInfo
#
#  pi = om.ProductionInfo(name=b'glyph-marker', position=(10, -5))
#     creates a new, empty ProductionInfo with some attributes set
#
#     here is the full list of attributes on om.ProductionInfo with some example values:
#  pi = om.ProductionInfo()
#  pi.shrink_left = False
#  pi.shrink_right = False
#  pi.isolate_inputs_from_outputs = False
#  pi.cabinets = [om.Cabinet(...), ...]
#  pi.conduits = [om.Conduit(...), ...]
#  pi.vials = [om.Vial(...), ...]

# === om.Cabinet ===
#  om.Cabinet represents a production cabinet
#
#  c = om.Cabinet(b'MediumWide', (0, 0))
#     creates a cabinet of type MediumWide at position (0, 0)
#
#  print(c.type, c.position)
#     prints the type and position of the cabinet

# === om.Conduit ===
#  om.Conduit represents a conduit in a production puzzle
#  note that this is different from conduits in *solutions*, which are represented as om.Parts
#
#  c = om.Conduit((-2, 0), (2, 0), hexes=[(0, 0)])
#     creates a new conduit with initial positions at (-2, 0) and (2, 0) and a single hex
#
#     here is the full list of attributes on om.Conduit with some example values:
#  c.starting_position_a = (0, 0)
#  c.starting_position_b = (4, 0)
#  c.hexes = [(0, 0), (1, 0)]

# === om.Vial ===
#  om.Vial represents a vial in a production puzzle
#
#  v = om.Vial(position=(2, 0), vial_count=2, flip_vertically=False)
#     creates a vial at position (2, 0), with 2 vials shown visually, positioned with the wider side downward
#
#  print(v.position, v.vial_count, v.flip_vertically)
#     accesses the properties of the vial to print them

# === om.Solution ===
#  om.Solution represents a solution file
#
#  sol = om.Solution()
#     creates a new, empty solution
#
#  sol = om.Solution('/path/to/file.solution')
#     loads a solution from a path
#
#  sol = om.Solution(b'...')
#     creates a solution from bytes in solution file format
#
#  data = sol.to_bytes()
#     encodes a solution as a byte array in solution file format and returns it
#
#  sol.write_to_path('/path/to/file.solution')
#     writes a solution to a path
#
#  sol_copy = om.Solution(sol)
#     creates a "deep copy" of a solution -- equivalent to om.Solution(sol.to_bytes())
#
#  sol = om.Solution(puzzle=b'P007', name=b'NEW SOLUTION 1')
#     creates a new, empty solution with some attributes set
#
#  sol = om.Solution()
#  sol.puzzle = b'P007'
#  sol.name = b'NEW SOLUTION 1'
#     creates a new, empty solution and sets some attributes afterward
#
#     here is the full list of attributes on om.Solution with some example values:
#  sol = om.Solution()
#  sol.puzzle = b'puzzle file name'
#  sol.name = b'solution name'
#  sol.solved = True  (or False)
#  sol.cycles = 123  (or None, if solved is False)
#  sol.cost = 250  (or None, if solved is False)
#  sol.area = 45  (or None, if solved is False)
#  sol.instructions = 51  (or None, if solved is False)
#  sol.parts = [om.Part(...), ...]

# === om.Part ===
#  om.Part represents an individual part in a solution file
#
#  part = om.Part()
#     creates a new, empty part
#
#  part = om.Part(name=om.Part.EQUILIBRIUM, position=(10, -5))
#     creates a new, empty part with some attributes set
#
#     here is the full list of attributes on om.Part with some example values:
#  part = om.Part()
#  part.name = om.Part.ANIMISMUS
#  part.position = (1, 2)
#  part.length = 2  (for arms)
#  part.rotation = 4
#  part.which_reagent_or_product = 0  (for inputs and outputs)
#  part.instructions = [om.Instruction(0, b'g'), om.Instruction(1, b'r'), om.Instruction(2, b'G'), om.Instruction(3, b'R')]
#  part.track_hexes = [(0, 0), (1, 0)]  (for parts named b'track')
#  part.arm_number = 0  (for arms)
#  part.conduit_id = 100  (for parts named b'pipe')
#  part.conduit_hexes = [(0, 0), (1, 0)]  (for parts named b'pipe')
#
#  here is the full list of part name constants:
#     om.Part.ARM1
#     om.Part.ARM2
#     om.Part.ARM3
#     om.Part.ARM6
#     om.Part.PISTON
#     om.Part.TRACK
#     om.Part.BERLO
#     om.Part.BONDER
#     om.Part.UNBONDER
#     om.Part.TRIPLEX
#     om.Part.MULTIBONDER
#     om.Part.CALCIFICATION
#     om.Part.DISPERSION
#     om.Part.DISPOSAL
#     om.Part.DUPLICATION
#     om.Part.ANIMISMUS
#     om.Part.EQUILIBRIUM
#     om.Part.PROJECTION
#     om.Part.PURIFICATION
#     om.Part.UNIFICATION
#     om.Part.INPUT
#     om.Part.OUTPUT_STANDARD
#     om.Part.OUTPUT_REPEATING
#     om.Part.CONDUIT

# === om.Instruction ===
#  om.Instruction represents an instruction
#
#  inst = om.Instruction(0, b'g')
#     creates an instruction at the given index with the given encoding byte
#
#  inst = om.Instruction(0, om.Instruction.DROP)
#     named constants are available as well
#
#  print(inst.index, inst.instruction)
#     prints the index and the encoding byte
#
#  here is the full list of instruction constants:
#     om.Instruction.ROTATE_CW
#     om.Instruction.ROTATE_CCW
#     om.Instruction.EXTEND
#     om.Instruction.RETRACT
#     om.Instruction.GRAB
#     om.Instruction.DROP
#     om.Instruction.PIVOT_CW
#     om.Instruction.PIVOT_CCW
#     om.Instruction.TRACK_PLUS
#     om.Instruction.TRACK_MINUS
#     om.Instruction.REPEAT
#     om.Instruction.RESET
#     om.Instruction.NOOP

# === om.Sim ===
#  om.Sim is a convenient wrapper around omsim's libverify
#
#  sim = om.Sim('/path/to/file.puzzle', '/path/to/file.solution')
#     creates a new sim by reading the puzzle and solution from files
#     requires either libverify.so or libverify.dll to be in the library search path (otherwise raises RuntimeError)
#     other errors are raised as om.SimError
#
#  sim = om.Sim(b'...', b'...')
#     creates a new sim by reading the puzzle and solution from bytes in puzzle and solution file formats
#
#  sim = om.Sim(om.Puzzle(...), om.Solution(...))
#     creates a new sim from om.Puzzle and om.Solution objects
#
#  sim = om.Sim('/path/to/file.puzzle', om.Solution(...))
#     these forms can be mixed and matched
#
#  result = sim.metric('cycles')
#     measures a metric (as listed on http://events.critelli.technology/static/metrics.html)
#     collisions and other errors are raised as om.SimError
#
#  result = sim.approximate_metric('cycles')
#     like sim.metric(), but supports "approximate" metrics like "per repetition^2 area"
#
#  result = sim.rate()
#     computes the rate of the solution as either a python Fraction value or a floating-point infinity value
#
#  growth_order, growth_rate = sim.area_at_infinity()
#     computes the area at infinity of a solution as a (growth_order, growth_rate) pair
#     for example, a solution that remains at 32 area forever will have a pair of (0, 32)
#     while a solution that grows by 6 area per output will have a pair of (1, 6)
#
#  intervals = sim.output_intervals()
#     measures the intervals between output drops (aka lexicographic cycles)
#     returns the intervals as an om.OutputIntervals object which can be compared with other such objects

# === om.SimError ===
#  om.SimError tracks the message, cycle, and location of a collision or error reported by om.Sim
#
#  sim = om.Sim('/path/to/file.puzzle', '/path/to/file.solution')
#  try:
#     result = sim.metric('cycles')
#  except om.SimError as err:
#     print(err.message, err.cycle, err.location[0], err.location[1])

# SPECIAL THANKS TO
#
#   F43nd1r, for documenting the OM puzzle and solution formats at <https://github.com/F43nd1r/omsp/blob/master/Formats.md>

import ctypes
from fractions import Fraction
import re
import struct
import pathlib

import numpy as np
import hexmath

from sharedtypes import *


class Decoder:

    def __init__(self, initial_bytes):
        self.bytes = memoryview(initial_bytes)

    def read_struct_format(self, struct_format):
        n = struct.calcsize(struct_format)
        if n > len(self.bytes):
            raise ValueError('not enough bytes left in file to parse value')
        result = struct.unpack_from(struct_format, self.bytes)
        self.bytes = self.bytes[n:]
        return result

    def read_string(self):
        n = 0
        shift = 0
        while True:
            byte, = self.read_struct_format('<B')
            n |= (byte & 0x7f) << shift
            shift += 7
            if byte & 0x80 == 0:
                break
        if n > len(self.bytes):
            raise ValueError('not enough bytes left in file to parse value')
        result = bytes(self.bytes[:n])
        self.bytes = self.bytes[n:]
        return result


class Encoder:

    def __init__(self):
        self.bytes = bytearray()

    def write_struct_format(self, struct_format, *values):
        offset = len(self.bytes)
        self.bytes.extend(bytes(struct.calcsize(struct_format)))
        struct.pack_into(struct_format, self.bytes, offset, *values)

    def write_string(self, string):
        n = len(string)
        while True:
            if n < 0x80:
                self.write_struct_format('<B', n)
                break
            else:
                self.write_struct_format('<B', 0x80 | (n & 0x7f))
                n >>= 7
        self.bytes.extend(string)


class Molecule:

    def __init__(self, *, decoder=None, atoms=None, bonds=None):
        self.atoms: list[Atom] = atoms or []
        self.bonds: list[Bond] = bonds or []
        self._hash = None
        if decoder is None:
            return
        natoms, = decoder.read_struct_format('<I')
        for i in range(natoms):
            t, u, v = decoder.read_struct_format('<Bbb')
            self.atoms.append(Atom(t, (u, v)))
        nbonds, = decoder.read_struct_format('<I')
        for i in range(nbonds):
            t, u0, v0, u1, v1 = decoder.read_struct_format('<Bbbbb')
            self.bonds.append(Bond(t, ((u0, v0), (u1, v1))))

    def copy(self):
        return Molecule(
            atoms=[atom.copy() for atom in self.atoms],
            bonds=[bond.copy() for bond in self.bonds]
        )

    def encode(self, encoder):
        encoder.write_struct_format('<I', len(self.atoms))
        for atom in self.atoms:
            encoder.write_struct_format(
                '<Bbb', atom.type, atom.position[0],
                atom.position[1]
            )
        encoder.write_struct_format('<I', len(self.bonds))
        for bond in self.bonds:
            encoder.write_struct_format(
                '<Bbbbb', bond.type,
                bond.positions[0][0], bond.positions[0][1],
                bond.positions[1][0], bond.positions[1][1]
            )

    def atom_types(self):
        return sorted(list({atom.type for atom in self.atoms}))

    def bond_types(self):
        return sorted(list({bond.type for bond in self.bonds}))

    def translate(self, position: Pos2D = (0, 0), rotation: Angle = 0):
        for atom in self.atoms:
            atom.translate(position, rotation)
        for bond in self.bonds:
            bond.translate(position, rotation)
        self._hash = None
        return self

    def move_atom_to(self, position: Pos2D, atom_i: int = 0):
        position = hexmath.difference(position, self.atoms[atom_i].position)
        for atom in self.atoms:
            atom.translate(position)
        for bond in self.bonds:
            bond.translate(position)
        self._hash = None
        return self

    def rotate(self, rotation: Angle, pivot_pos: Pos2D = (0, 0)):
        for atom in self.atoms:
            atom.rotate(rotation, pivot_pos)
        for bond in self.bonds:
            bond.rotate(rotation, pivot_pos)
        self._hash = None
        return self

    def translate_center(self):
        min_x = min(atom.position[0] for atom in self.atoms)
        min_y = min(atom.position[1] for atom in self.atoms)
        return self.translate(position=(-min_x, -min_y))

    def matches(self, other):
        rotated_copies = [
            self.copy().rotate(rot).translate_center()
            for rot in range(6)
        ]
        centered_other_copy = other.copy().translate_center()
        return any(centered_other_copy == rotated_copy
            for rotated_copy in rotated_copies)

    def angle_agnostic(self):
        rotated_copies = [
            self.copy().rotate(rot).translate_center()
            for rot in range(6)
        ]
        return min(rotated_copies)

    def symmetry_angle(self):
        centered_copy = self.translate_center()
        rotated_copies = [
            self.copy().rotate(rot).translate_center()
            for rot in range(6)
        ]
        equalities = [centered_copy == rot_copy for rot_copy in rotated_copies]
        return 6 // sum(equalities)

    def overlap(self, space_s):
        # convert space_s to set of tuples for various types of input
        if isinstance(space_s, Molecule):
            space_s = set(atom.position for atom in space_s.atoms)
        elif isinstance(space_s, np.ndarray):
            if len(space_s.shape) > 1:
                space_s = set(tuple(x) for x in space_s.tolist())
            else:
                space_s = {tuple(space_s.tolist())}
        elif isinstance(space_s, list):
            if isinstance(space_s[0], tuple):
                space_s = set(space_s)
            elif isinstance(space_s[0], list):
                space_s = set(tuple(x) for x in space_s)
            elif isinstance(space_s[0], int):
                space_s = set(tuple(space_s))
            else:
                raise TypeError(
                    "space_s is a list of something I don't "
                    "recognize: %r" % space_s
                )
        elif isinstance(space_s, tuple):
            if isinstance(space_s[0], int):
                return any(atom.position == space_s for atom in self.atoms)
            else:
                raise TypeError("space_s is a tuple but not of int")
        else:
            raise TypeError(
                "space_s must be a Molecule, numpy array, list, "
                "or tuple: space_s=%r" % space_s
            )
        # check for overlap
        return any(atom.position in space_s for atom in self.atoms)

    def find_disconnected(self):
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
        if len(groups) == 1:
            return [self]
        molecules = [
            Molecule(
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

    def __eq__(self, other):
        return (
                isinstance(other, Molecule)
                and set(self.atoms) == set(other.atoms)
                and set(self.bonds) == set(other.bonds)
        )

    def __lt__(self, other):
        assert isinstance(other, Molecule)
        if len(self.atoms) != len(other.atoms):
            return len(self.atoms) < len(other.atoms)
        if len(self.bonds) != len(other.bonds):
            return len(self.bonds) < len(other.bonds)
        s_atoms_sorted = sorted(self.atoms, key=lambda atom: atom.position)
        o_atoms_sorted = sorted(other.atoms, key=lambda atom: atom.position)
        for s_atom, o_atom in zip(s_atoms_sorted, o_atoms_sorted):
            if s_atom != o_atom:
                return s_atom < o_atom
        s_bonds_sorted = sorted(self.bonds, key=lambda bond: bond.positions)
        o_bonds_sorted = sorted(other.bonds, key=lambda bond: bond.positions)
        for s_bond, o_bond in zip(s_bonds_sorted, o_bonds_sorted):
            if s_bond != o_bond:
                return s_bond < o_bond
        return False

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((tuple(self.atoms), tuple(self.bonds)))
        return self._hash

    def __repr__(self):
        return "Molecule(%r, %r)" % (self.atoms, self.bonds)

    def __str__(self):
        return "Molecule([%s], [%s])" % (', '.join(str(a) for a in self.atoms),
        ', '.join(str(b) for b in self.bonds))


class Atom:
    SALT = 1
    AIR = 2
    EARTH = 3
    FIRE = 4
    WATER = 5
    QUICKSILVER = 6
    GOLD = 7
    SILVER = 8
    COPPER = 9
    IRON = 10
    TIN = 11
    LEAD = 12
    VITAE = 13
    MORS = 14
    REPETITION_PLACEHOLDER = 15
    QUINTESSENCE = 16

    TYPE_NAMES = {
        SALT: "SALT",
        AIR: "AIR",
        EARTH: "EARTH",
        FIRE: "FIRE",
        WATER: "WATER",
        QUICKSILVER: "QUICKSILVER",
        GOLD: "GOLD",
        SILVER: "SILVER",
        COPPER: "COPPER",
        IRON: "IRON",
        TIN: "TIN",
        LEAD: "LEAD",
        VITAE: "VITAE",
        MORS: "MORS",
        REPETITION_PLACEHOLDER: "REPETITION_PLACEHOLDER",
        QUINTESSENCE: "QUINTESSENCE",
    }

    TYPES = list(range(1, 15)) + [16]  # skip repetition symbol
    TYPES_ELEMENTAL = {AIR, EARTH, FIRE, WATER}
    TYPES_METAL = {LEAD, TIN, IRON, COPPER, SILVER, GOLD}

    TYPES_CATEGORY_ELEMENTAL = TYPES_ELEMENTAL | {VITAE, MORS, SALT,
        QUINTESSENCE}
    TYPES_CATEGORY_METAL = TYPES_METAL | {QUICKSILVER}

    def __init__(self, type, position: Pos2D):
        self.type: int = type
        self.position: Pos2D = tuple(position)
        self._hash = None

    def copy(self):
        return Atom(self.type, self.position)

    def translate(self, position=(0, 0), rotation=0):
        self.position = hexmath.translate(self.position, position, rotation)
        self._hash = None

    def rotate(self, rotation: Angle, pivot: Pos2D = (0, 0)):
        self.position = hexmath.rotate(self.position, rotation, pivot)
        self._hash = None

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return False
        return self.type == other.type and self.position == other.position

    def __lt__(self, other):
        assert isinstance(other, Atom)
        return (self.type, self.position) < (other.type, other.position)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((self.type, self.position))
        return self._hash

    def __repr__(self):
        return "Atom(%r, %r)" % (self.type, self.position)

    def __str__(self):
        return "Atom(%s, %s)" % (Atom.TYPE_NAMES[self.type], self.position)


class Bond:
    NORMAL = 1 << 0
    TRIPLEX_RED = 1 << 1
    TRIPLEX_BLACK = 1 << 2
    TRIPLEX_YELLOW = 1 << 3
    TRIPLEX = TRIPLEX_RED | TRIPLEX_BLACK | TRIPLEX_YELLOW
    TYPES = {
        # only 1 type of normal
        NORMAL,
        # 7 types of triplex (2^3 - 1 (the none one))
        # all 3 together
        TRIPLEX,
        # just 1
        TRIPLEX_RED,
        TRIPLEX_BLACK,
        TRIPLEX_YELLOW,
        # just 2
        TRIPLEX_RED | TRIPLEX_BLACK,
        TRIPLEX_BLACK | TRIPLEX_YELLOW,
        TRIPLEX_YELLOW | TRIPLEX_RED,
    }

    @staticmethod
    def type_name(bond_type):
        types = []
        if bond_type & Bond.NORMAL:
            types.append("Normal")
        if bond_type & Bond.TRIPLEX:
            types.append("Triplex")
        else:
            if bond_type & Bond.TRIPLEX_RED:
                types.append("Red")
            if bond_type & Bond.TRIPLEX_BLACK:
                types.append("Black")
            if bond_type & Bond.TRIPLEX_YELLOW:
                types.append("Yellow")

        if bond_type & ~(Bond.NORMAL | Bond.TRIPLEX):
            types.append(
                "Unknown(%d)" % (
                        bond_type & ~(Bond.NORMAL | Bond.TRIPLEX))
            )

        if len(types) > 1:
            type_str = "[%s]" % ", ".join(types)
        else:
            type_str = types[0]
        return type_str

    @staticmethod
    def type_str_id(bond_type):
        types = []
        if bond_type & Bond.NORMAL:
            types.append("normal")
        if bond_type & Bond.TRIPLEX:
            types.append("triplex")
        else:
            if bond_type & Bond.TRIPLEX_RED:
                types.append("red")
            if bond_type & Bond.TRIPLEX_BLACK:
                types.append("black")
            if bond_type & Bond.TRIPLEX_YELLOW:
                types.append("yellow")

        unknown_flags = bond_type & ~(Bond.NORMAL | Bond.TRIPLEX)
        if unknown_flags:
            types.append(
                "Unknown(%d)" % (unknown_flags)
            )

        if len(types) > 1:
            type_str = "_".join(types)
        else:
            type_str = types[0]
        return type_str

    def __init__(self, type: int, positions: tuple[tuple, tuple]):
        self.type: int = type
        self.positions: tuple[tuple[int, int], tuple[int, int]] = \
            tuple(sorted(positions))
        self._hash = None

    def copy(self):
        return Bond(self.type, self.positions)

    def translate(self, position=(0, 0), rotation=0):
        self.positions = tuple(sorted((
            hexmath.translate(self.positions[0], position, rotation),
            hexmath.translate(self.positions[1], position, rotation)
        )))
        self._hash = None

    def rotate(self, rotation: Angle, pivot: Pos2D = (0, 0)):
        self.positions = tuple(sorted((
            hexmath.rotate(self.positions[0], rotation, pivot),
            hexmath.rotate(self.positions[1], rotation, pivot)
        )))
        self._hash = None

    def other(self, position):
        if self.positions[0] == position:
            return self.positions[1]
        else:
            return self.positions[0]

    def __eq__(self, other):
        if not isinstance(other, Bond):
            return False
        return (self.type == other.type
                and self.positions == other.positions)

    def __lt__(self, other):
        assert isinstance(other, Bond)
        return (self.type, self.positions) < (other.type, other.positions)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((self.type, self.positions))
        return self._hash

    def __repr__(self):
        return "Bond(%r, %r)" % (self.type, self.positions)

    def __str__(self):

        return "Bond(%s, %s-%s)" % (Bond.type_name(self.type),
        self.positions[0], self.positions[1])


class ProductionInfo:

    def __init__(
            self, *, decoder=None, shrink_left=False, shrink_right=False,
            isolate_inputs_from_outputs=False, cabinets=None, conduits=None,
            vials=None
    ):
        self.shrink_left = shrink_left
        self.shrink_right = shrink_right
        self.isolate_inputs_from_outputs = isolate_inputs_from_outputs
        self.cabinets = cabinets or []
        self.conduits = conduits or []
        self.vials = vials or []
        if decoder is None:
            return
        self.shrink_left, self.shrink_right, self.isolate_inputs_from_outputs, ncabinets = decoder.read_struct_format(
            '<???I'
        )
        for i in range(ncabinets):
            self.cabinets.append(
                Cabinet(
                    position=decoder.read_struct_format('<bb'),
                    type=decoder.read_string()
                )
            )
        nconduits, = decoder.read_struct_format('<I')
        for i in range(nconduits):
            self.conduits.append(Conduit(decoder=decoder))
        nvials, = decoder.read_struct_format('<I')
        for i in range(nvials):
            u, v, flip, count = decoder.read_struct_format('<bb?I')
            self.vials.append(Vial((u, v), count, flip))

    def encode(self, encoder):
        encoder.write_struct_format(
            '<???I', self.shrink_left,
            self.shrink_right, self.isolate_inputs_from_outputs,
            len(self.cabinets)
        )
        for cabinet in self.cabinets:
            encoder.write_struct_format('<bb', *cabinet.position)
            encoder.write_string(cabinet.type)
        encoder.write_struct_format('<I', len(self.conduits))
        for conduit in self.conduits:
            conduit.encode(encoder)
        encoder.write_struct_format('<I', len(self.vials))
        for vial in self.vials:
            encoder.write_struct_format(
                '<bb?I', vial.position[0],
                vial.position[1], vial.flip_vertically, vial.vial_count
            )


class Cabinet:

    def __init__(self, type, position):
        self.type = type
        self.position = position


class Conduit:

    def __init__(
            self, starting_position_a=(0, 0), starting_position_b=(0, 0),
            *, decoder=None, hexes=None
    ):
        self.starting_position_a = starting_position_a
        self.starting_position_b = starting_position_b
        self.hexes = hexes or []
        if decoder is None:
            return
        u0, v0, u1, v1, nhexes = decoder.read_struct_format('<bbbbI')
        self.starting_position_a = (u0, v0)
        self.starting_position_b = (u1, v1)
        for i in range(nhexes):
            self.hexes.append(decoder.read_struct_format('<bb'))

    def encode(self, encoder):
        encoder.write_struct_format(
            '<bbbbI', self.starting_position_a[0],
            self.starting_position_a[1], self.starting_position_b[0],
            self.starting_position_b[1], len(self.hexes)
        )
        for h in self.hexes:
            encoder.write_struct_format('<bb', *h)


class Vial:

    def __init__(self, position, vial_count, flip_vertically):
        self.position = position
        self.vial_count = vial_count
        self.flip_vertically = flip_vertically


class Solution:

    def __init__(
            self, file=None, *, decoder=None, puzzle=b'', name=b'',
            solved=False, cycles=None, cost=None, area=None, instructions=None,
            parts=None
    ):
        self.puzzle = puzzle
        self.name = name
        self.solved = solved
        self.cycles = cycles
        self.cost = cost
        self.area = area
        self.instructions = instructions
        self.parts = parts or []
        if isinstance(file, str):
            with open(file, 'rb') as f:
                decoder = Decoder(f.read())
        elif isinstance(file, bytes):
            decoder = Decoder(file)
        elif isinstance(file, Solution):
            decoder = Decoder(file.to_bytes())
        elif file is not None:
            raise ValueError(
                'Solution() may only take a filename, bytes, or another Solution as a positional parameter'
            )
        if decoder is None:
            return
        if decoder.read_struct_format('<I') != (7,):
            raise ValueError('unknown version number in solution file')
        self.puzzle = decoder.read_string()
        self.name = decoder.read_string()
        nmetrics, = decoder.read_struct_format('<I')
        if nmetrics == 4:
            self.solved = True
            zero, self.cycles, one, self.cost, two, self.area, three, self.instructions = decoder.read_struct_format(
                '<IIIIIIII'
            )
            if (zero, one, two, three) != (0, 1, 2, 3):
                raise ValueError(
                    'unexpected numbering of metrics in solution file'
                )
        elif nmetrics != 0:
            raise ValueError(
                'wrong number of metrics in solution file (expecting 0 or 4)'
            )
        nparts, = decoder.read_struct_format('<I')
        for i in range(nparts):
            self.parts.append(Part(decoder=decoder))

    def get_arm_parts(self) -> list:
        return sorted(
            [part for part in self.parts
                if part.name in Part.PARTS_PROGRAMMABLE],
            key=lambda part: part.arm_number
        )

    def encode(self, encoder):
        encoder.write_struct_format('<I', 7)
        encoder.write_string(self.puzzle)
        encoder.write_string(self.name)
        if self.solved:
            encoder.write_struct_format(
                '<IIIIIIIII', 4, 0, self.cycles, 1,
                self.cost, 2, self.area, 3, self.instructions
            )
        else:
            encoder.write_struct_format('<I', 0)
        encoder.write_struct_format('<I', len(self.parts))
        for part in self.parts:
            part.encode(encoder)

    def to_bytes(self):
        encoder = Encoder()
        self.encode(encoder)
        return encoder.bytes

    def write_to_path(self, filename):
        with open(filename, 'wb') as f:
            f.write(self.to_bytes())


class Part:
    ARM1 = b'arm1'
    ARM2 = b'arm2'
    ARM3 = b'arm3'
    ARM6 = b'arm6'
    PISTON = b'piston'
    TRACK = b'track'
    BERLO = b'baron'
    BONDER = b'bonder'
    UNBONDER = b'unbonder'
    TRIPLEX = b'bonder-prisma'
    MULTIBONDER = b'bonder-speed'
    CALCIFICATION = b'glyph-calcification'
    DISPERSION = b'glyph-dispersion'
    DISPOSAL = b'glyph-disposal'
    DUPLICATION = b'glyph-duplication'
    ANIMISMUS = b'glyph-life-and-death'
    EQUILIBRIUM = b'glyph-marker'
    PROJECTION = b'glyph-projection'
    PURIFICATION = b'glyph-purification'
    UNIFICATION = b'glyph-unification'
    INPUT = b'input'
    OUTPUT_STANDARD = b'out-std'
    OUTPUT_REPEATING = b'out-rep'
    CONDUIT = b'pipe'

    PARTS_PROGRAMMABLE = {ARM1, ARM2, ARM3, ARM6, PISTON, BERLO}
    PARTS_ARMS = {ARM1, ARM2, ARM3, ARM6, PISTON}
    PARTS_POST_MOVEMENT = {BONDER, TRIPLEX, MULTIBONDER, UNBONDER,
        CALCIFICATION, PROJECTION, DUPLICATION, INPUT, OUTPUT_REPEATING}
    PARTS_CONSUMERS = {PROJECTION, PURIFICATION, UNIFICATION, DISPERSION,
        ANIMISMUS, DISPOSAL, OUTPUT_STANDARD, CONDUIT}
    PARTS_OUTPUTS = {OUTPUT_STANDARD, OUTPUT_REPEATING}
    PARTS_GLYPHS = {CALCIFICATION, DISPERSION, DUPLICATION, ANIMISMUS,
        PROJECTION, PURIFICATION, UNIFICATION, BERLO, DISPOSAL}
    PARTS_BONDS = {BONDER, UNBONDER, MULTIBONDER, TRIPLEX}

    COSTS = {
        ARM1: 20,
        ARM2: 30,
        ARM3: 30,
        ARM6: 30,
        PISTON: 40,
        TRACK: 5,
        BERLO: 30,
        BONDER: 10,
        UNBONDER: 10,
        TRIPLEX: 20,
        MULTIBONDER: 30,
        CALCIFICATION: 10,
        DISPERSION: 20,
        DISPOSAL: 0,
        DUPLICATION: 20,
        ANIMISMUS: 20,
        EQUILIBRIUM: 0,
        PROJECTION: 20,
        PURIFICATION: 20,
        UNIFICATION: 20,
        INPUT: 0,
        OUTPUT_STANDARD: 0,
        OUTPUT_REPEATING: 0,
        CONDUIT: 0,
    }

    AREAS = {
        ARM1: 1,
        ARM2: 1,
        ARM3: 1,
        ARM6: 1,
        PISTON: 1,
        TRACK: 1,
        BERLO: 1,
        BONDER: 2,
        UNBONDER: 2,
        TRIPLEX: 3,
        MULTIBONDER: 4,
        CALCIFICATION: 1,
        DISPERSION: 5,
        DISPOSAL: 7,
        DUPLICATION: 2,
        ANIMISMUS: 4,
        EQUILIBRIUM: 1,
        PROJECTION: 2,
        PURIFICATION: 3,
        UNIFICATION: 5,
        INPUT: 1,
        OUTPUT_STANDARD: 1,
        OUTPUT_REPEATING: 1,
        CONDUIT: 1,
    }

    COMPATIBLE_MIN_ARM_LENGTH = {
        BONDER, UNBONDER, CALCIFICATION, DUPLICATION, PROJECTION,
        INPUT, OUTPUT_STANDARD, OUTPUT_REPEATING, CONDUIT,
    }

    ACCESSES = {
        ARM1: [0],
        ARM2: [0],
        ARM3: [0],
        ARM6: [0],
        PISTON: [0],
        TRACK: [0],
        BERLO: [0],
        BONDER: list(range(3)),
        UNBONDER: list(range(3)),
        TRIPLEX: list(range(4)),
        MULTIBONDER: list(range(5)),
        CALCIFICATION: list(range(2)),
        DISPERSION: [5],
        DISPOSAL: [1],
        DUPLICATION: list(range(3)),
        ANIMISMUS: [4],
        EQUILIBRIUM: [0],
        PROJECTION: list(range(1, 3)),
        PURIFICATION: [3],
        UNIFICATION: [5],
        INPUT: [1],
        OUTPUT_STANDARD: [1],
        OUTPUT_REPEATING: list(range(2)),
        CONDUIT: [1],
    }

    @staticmethod
    def channels(type, armlength):
        if armlength < 1:
            raise ValueError("armlength must be a positive integer")

        if armlength == 1:
            return 1
        elif armlength == 2:
            d = {
                Part.BONDER: 2,
                Part.UNBONDER: 2,
                Part.TRIPLEX: 3,
                Part.MULTIBONDER: 4,
                Part.CALCIFICATION: 1,
                Part.DISPERSION: 4,
                Part.DISPOSAL: 1,
                Part.DUPLICATION: 2,
                Part.ANIMISMUS: 4,
                Part.PROJECTION: 2,
                Part.PURIFICATION: 3,
                Part.UNIFICATION: 3,
                Part.INPUT: 1,
                Part.OUTPUT_STANDARD: 1,
                Part.OUTPUT_REPEATING: 1,
            }
            return d[type]
        elif armlength >= 3:
            d = {
                Part.BONDER: 2,
                Part.UNBONDER: 2,
                Part.TRIPLEX: 3,
                Part.MULTIBONDER: 4,
                Part.CALCIFICATION: 1,
                Part.DISPERSION: 5,
                Part.DISPOSAL: 1,
                Part.DUPLICATION: 2,
                Part.ANIMISMUS: 4,
                Part.PROJECTION: 2,
                Part.PURIFICATION: 3,
                Part.UNIFICATION: 5,
                Part.INPUT: 1,
                Part.OUTPUT_STANDARD: 1,
                Part.OUTPUT_REPEATING: 1,
            }
            return d[type]

    # tick part evaluation order:
    # 1. grabs/drops
    # 2. consumers trigger
    # 3. rotates/translates
    # 4. post movement parts trigger (bonders and mutators)
    # TECHNICALLY post movement parts can trigger during step 2 as well, but it
    # never happens in base game with legal solves

    def __init__(
            self, *, decoder=None, name=b'', position=(0, 0), length=0,
            rotation=0, which_reagent_or_product=0, instructions=None,
            track_hexes=None, arm_number=0, conduit_id=0, conduit_hexes=None,
            grabbing=False, grabbed=None
    ):
        self._hash = None
        self.name = name
        self.position = position
        self.length = length
        self.rotation = rotation
        self.which_reagent_or_product = which_reagent_or_product
        self.instructions = instructions or []
        self.track_hexes = track_hexes or []
        self.arm_number = arm_number
        self.grabbing = grabbing
        self.grabbed: list[bool] = grabbed
        self.conduit_id = conduit_id
        self.conduit_hexes = conduit_hexes or []
        if decoder is None:
            return
        self.name = decoder.read_string()
        if decoder.read_struct_format('<B') != (1,):
            raise ValueError('unknown part version number in solution file')
        self.position = decoder.read_struct_format('<ii')
        self.length, self.rotation, self.which_reagent_or_product, ninstrs = decoder.read_struct_format(
            '<IiII'
        )
        for i in range(ninstrs):
            self.instructions.append(
                Instruction(*decoder.read_struct_format('<ic'))
            )
        if self.name == b'track':
            ntrack, = decoder.read_struct_format('<I')
            for i in range(ntrack):
                self.track_hexes.append(decoder.read_struct_format('<ii'))
        self.arm_number, = decoder.read_struct_format('<I')
        if self.name == b'pipe':
            self.conduit_id, npipe = decoder.read_struct_format('<II')
            for i in range(npipe):
                self.conduit_hexes.append(decoder.read_struct_format('<ii'))

    def copy(self):
        return Part(
            name=self.name, position=self.position, length=self.length,
            rotation=self.rotation,
            which_reagent_or_product=self.which_reagent_or_product,
            instructions=self.instructions, track_hexes=self.track_hexes,
            arm_number=self.arm_number,
            conduit_id=self.conduit_id, conduit_hexes=self.conduit_hexes,
            grabbing=self.grabbing, grabbed=self.grabbed
        )

    def encode(self, encoder):
        encoder.write_string(self.name)
        encoder.write_struct_format(
            '<BiiIiII', 1, self.position[0],
            self.position[1], self.length, self.rotation,
            self.which_reagent_or_product, len(self.instructions)
        )
        for instruction in self.instructions:
            encoder.write_struct_format(
                '<ic', instruction.index,
                instruction.instruction
            )
        if self.name == b'track':
            encoder.write_struct_format('<I', len(self.track_hexes))
            for track_hex in self.track_hexes:
                encoder.write_struct_format('<ii', track_hex[0], track_hex[1])
        encoder.write_struct_format('<I', self.arm_number)
        if self.name == b'pipe':
            encoder.write_struct_format(
                '<II', self.conduit_id,
                len(self.conduit_hexes)
            )
            for conduit_hex in self.conduit_hexes:
                encoder.write_struct_format(
                    '<ii', conduit_hex[0],
                    conduit_hex[1]
                )

    def is_programmable(self):
        return self.name in Part.PARTS_PROGRAMMABLE

    def valid_instructions(self, overlaps_track: bool = False):
        if not self.is_programmable():
            return []

        default = [None]
        if overlaps_track:
            default.append(Instruction.GROUP_TRACK)

        if self.name == Part.BERLO:
            return default + Instruction.GROUP_BERLO
        elif self.name == Part.PISTON:
            return default + Instruction.GROUP_PISTON
        else:
            return default + Instruction.GROUP_ARM

    def arm_grab_positions(self) -> list[tuple[int, int]]:
        if self.name not in Part.PARTS_ARMS:
            raise ValueError("This part is not an arm")
        count = 1
        try:
            count = int(chr(self.name[-1]))
        except ValueError:
            pass

        rotations = [(self.rotation + 6 // count * i) % 6 for i in range(count)]
        return [
            hexmath.summate(self.position, hexmath.scale(
                hexmath.ROTATION_VECTORS[rotation % 6], self.length))
            for rotation in rotations
        ]

    def __eq__(self, other):
        if not isinstance(other, Part):
            return False
        return (
                self.name == other.name and
                self.position == other.position and
                self.length == other.length and
                self.rotation == other.rotation and
                self.which_reagent_or_product == other.which_reagent_or_product and
                self.instructions == other.instructions and
                self.track_hexes == other.track_hexes and
                self.arm_number == other.arm_number and
                self.conduit_id == other.conduit_id and
                self.conduit_hexes == other.conduit_hexes and
                self.grabbing == other.grabbing and
                self.grabbed == other.grabbed
        )

    def __str__(self):
        aspects = [
            "pos=%s" % self.position,
            "rot=%d" % self.rotation,
        ]
        if self.name in Part.PARTS_ARMS:
            aspects += "arm_number=%d" % self.arm_number
            aspects += "length=%d" % self.length
            aspects += "grabbing=%s" % (
                False if self.grabbing else self.grabbed)
        elif self.name == Part.INPUT:
            aspects += "reagent=%d" % self.which_reagent_or_product
        elif self.name in {Part.OUTPUT_STANDARD, Part.OUTPUT_REPEATING}:
            aspects += "product=%d" % self.which_reagent_or_product
        elif self.name == Part.TRACK:
            aspects += "index=%d" % self.track_hexes.index(self.position)
            aspects += "track=%s" % self.track_hexes
        elif self.name == Part.CONDUIT:
            aspects += "id=%s" % self.conduit_id
            aspects += "spaces=%s" % self.conduit_hexes
        return "Part(%s, %s)" % (
            self.name, ", ".join(aspects)
        )

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((
                self.name, self.position, self.length, self.rotation,
                self.which_reagent_or_product, tuple(self.instructions),
                tuple(self.track_hexes), self.arm_number, self.conduit_id,
                tuple(self.conduit_hexes), self.grabbing,
                (tuple(self.grabbed) if isinstance(self.grabbed, list)
                 else self.grabbed)
            ))
        return self._hash


class Instruction:
    ROTATE_CW = b'R'
    ROTATE_CCW = b'r'
    EXTEND = b'E'
    RETRACT = b'e'
    GRAB = b'G'
    DROP = b'g'
    PIVOT_CW = b'P'
    PIVOT_CCW = b'p'
    TRACK_PLUS = b'A'
    TRACK_MINUS = b'a'
    REPEAT = b'C'
    RESET = b'X'
    NOOP = b'O'

    OP_NAMES = {
        ROTATE_CW: "R CW",
        ROTATE_CCW: "R CCW",
        EXTEND: "Extend",
        RETRACT: "Retract",
        GRAB: "Grab",
        DROP: "Drop",
        PIVOT_CW: "Pivot CW",
        PIVOT_CCW: "Pivot CCW",
        TRACK_PLUS: "Track Plus",
        TRACK_MINUS: "Track Minus",
        REPEAT: "Repeat",
        RESET: "Reset",
        NOOP: "Noop",
    }

    GROUP_BERLO = [ROTATE_CCW, ROTATE_CW]
    GROUP_ARM = GROUP_BERLO + [GRAB, DROP, PIVOT_CCW, PIVOT_CW]
    GROUP_PISTON = GROUP_ARM + [EXTEND, RETRACT]
    GROUP_TRACK = [TRACK_PLUS, TRACK_MINUS]

    def __init__(self, index, instruction):
        self.index = index
        self.instruction = instruction

    def __repr__(self):
        return "Instruction(%r, %r)" % (self.index, self.instruction)

    def __str__(self):
        return "Instruction(%s, %s)" % (self.index,
        Instruction.OP_NAMES[self.instruction])


class Puzzle:
    # constants for parts_available:
    ARM = 1 << 0
    MULTIARM = 1 << 1
    PISTON = 1 << 2
    TRACK = 1 << 3
    BONDER = 1 << 8
    UNBONDER = 1 << 9
    MULTIBONDER = 1 << 10
    TRIPLEX = 1 << 11
    CALCIFICATION = 1 << 12
    DUPLICATION = 1 << 13
    PROJECTION = 1 << 14
    PURIFICATION = 1 << 15
    ANIMISMUS = 1 << 16
    DISPOSAL = 1 << 17
    QUINTESSENCE = 1 << 18
    GRAB_AND_ROTATE = 1 << 22
    DROP = 1 << 23
    RESET = 1 << 24
    REPEAT = 1 << 25
    PIVOT = 1 << 26
    BERLO = 1 << 28
    DEFAULT_PARTS_AVAILABLE = 0x07C0170F

    FLAG_NAMES = {
        0: "ARM",
        1: "MULTIARM",
        2: "PISTON",
        3: "TRACK",
        8: "BONDER",
        9: "UNBONDER",
        10: "MULTIBONDER",
        11: "TRIPLEX",
        12: "CALCIFICATION",
        13: "DUPLICATION",
        14: "PROJECTION",
        15: "PURIFICATION",
        16: "ANIMISMUS",
        17: "DISPOSAL",
        18: "QUINTESSENCE",
        22: "GRAB_AND_ROTATE",
        23: "DROP",
        24: "RESET",
        25: "REPEAT",
        26: "PIVOT",
        28: "BERLO",
    }

    FLAG_PARTS = {
        0: [Part.ARM1],
        1: [Part.ARM2, Part.ARM3, Part.ARM6],
        2: [Part.PISTON],
        3: [Part.TRACK],
        8: [Part.BONDER],
        9: [Part.UNBONDER],
        10: [Part.MULTIBONDER],
        11: [Part.TRIPLEX],
        12: [Part.CALCIFICATION],
        13: [Part.DUPLICATION],
        14: [Part.PROJECTION],
        15: [Part.PURIFICATION],
        16: [Part.ANIMISMUS],
        17: [Part.DISPOSAL],
        18: [Part.UNIFICATION, Part.DISPERSION],
        28: [Part.BERLO],
    }

    def __init__(
            self, file=None, *, decoder=None, name=b'', creator=0,
            parts_available=DEFAULT_PARTS_AVAILABLE, reagents=None,
            products=None, output_scale=1, production_info=None
    ):
        self._hash = None
        self.name = name
        self.creator = creator
        self.parts_available = parts_available
        self.reagents = reagents or []
        self.products = products or []
        self.output_scale = output_scale
        self.production_info = production_info
        if isinstance(file, str):
            with open(file, 'rb') as f:
                decoder = Decoder(f.read())
        elif isinstance(file, bytes):
            decoder = Decoder(file)
        elif isinstance(file, Puzzle):
            decoder = Decoder(file.to_bytes())
        elif file is not None:
            raise ValueError(
                'Puzzle() may only take a filename, bytes, or another Puzzle as a positional parameter'
            )
        if decoder is None:
            return
        if decoder.read_struct_format('<I') != (3,):
            raise ValueError('unknown version number in puzzle file')
        self.name = decoder.read_string()
        self.creator, self.parts_available, nreagents = decoder.read_struct_format(
            '<QQI'
        )
        for i in range(nreagents):
            self.reagents.append(Molecule(decoder=decoder))
        nproducts, = decoder.read_struct_format('<I')
        for i in range(nproducts):
            self.products.append(Molecule(decoder=decoder))
        self.output_scale, = decoder.read_struct_format('<I')
        if decoder.read_struct_format('<B') != (0,):
            self.production_info = ProductionInfo(decoder=decoder)

    def encode(self, encoder):
        encoder.write_struct_format('<I', 3)
        encoder.write_string(self.name)
        encoder.write_struct_format(
            '<QQI', self.creator, self.parts_available,
            len(self.reagents)
        )
        for m in self.reagents:
            m.encode(encoder)
        encoder.write_struct_format('<I', len(self.products))
        for m in self.products:
            m.encode(encoder)
        encoder.write_struct_format('<I', self.output_scale)
        if self.production_info is not None:
            encoder.write_struct_format('<B', 1)
            self.production_info.encode(encoder)
        else:
            encoder.write_struct_format('<B', 0)

    def to_bytes(self):
        encoder = Encoder()
        self.encode(encoder)
        return encoder.bytes

    def write_to_path(self, filename):
        with open(filename, 'wb') as f:
            f.write(self.to_bytes())

    def are_parts_available(self, parts):
        return self.parts_available & parts == parts

    def are_parts_banned(self, parts):
        return ~self.parts_available & parts == parts

    def flag_names(self):
        return [
            Puzzle.FLAG_NAMES[i]
            if i in Puzzle.FLAG_NAMES
            else str(i)
            for i in self.all_puzzle_flags()
        ]

    def all_puzzle_flags(self):
        return [i for i in range(32) if 1 << i & self.parts_available]

    def nondef_puzzle_flags(self):
        return [i for i in range(32)
            if
            1 << i & self.parts_available & ~Puzzle.DEFAULT_PARTS_AVAILABLE]

    def full_parts_list(self) -> list[bytes]:
        return [
            part
            for flag in self.all_puzzle_flags()
            if flag in Puzzle.FLAG_PARTS
            for part in Puzzle.FLAG_PARTS[flag]
        ]

    def nondef_parts_list(self) -> list[bytes]:
        return [
            part
            for flag in self.nondef_puzzle_flags()
            if flag in Puzzle.FLAG_PARTS
            for part in Puzzle.FLAG_PARTS[flag]
        ]

    def molecules(self):
        return self.reagents + self.products

    def atom_types(self) -> list[int]:
        return sorted(list({
            atom.type
            for mole in self.molecules()
            for atom in mole.atoms
        }))

    def reagent_types(self) -> list[int]:
        return sorted(list({
            atom.type
            for mole in self.reagents
            for atom in mole.atoms
        }))

    def product_types(self) -> list[int]:
        return sorted(list({
            atom.type
            for mole in self.products
            for atom in mole.atoms
        }))

    def bond_types(self) -> list[int]:
        return sorted(list({
            bond.type
            for mole in self.molecules()
            for bond in mole.bonds
        }))

    def __hash__(self):
        if self._hash is None:
            self._hash = hash((
                self.name, self.creator, self.parts_available,
                tuple(self.reagents),
                tuple(self.products), self.output_scale
            ))
        return self._hash

    def __str__(self):
        return "Puzzle(%s)" % self.name


class Sim:
    lv = None

    @classmethod
    def libverify(cls):
        if cls.lv is None:
            try:
                cls.lv = ctypes.cdll.LoadLibrary(
                    str(
                        pathlib.Path.cwd() /
                        'libverify.so'
                    )
                )
            except OSError:
                try:
                    cls.lv = ctypes.cdll.LoadLibrary(
                        str(
                            pathlib.Path.cwd() /
                            'libverify.dll'
                        )
                    )
                except OSError:
                    raise RuntimeError(
                        'unable to find libverify.so or libverify.dll -- to use om.Sim, download <https://github.com/ianh/omsim>, use `make` to build the library, and place it in the search path or working directory'
                    )
            cls.lv.verifier_create_from_bytes.restype = ctypes.c_void_p
            cls.lv.verifier_error.restype = ctypes.c_char_p
            cls.lv.verifier_evaluate_approximate_metric.restype = ctypes.c_double
        return cls.lv

    def __init__(self, puzzle, solution):
        puzzle_bytes = puzzle
        if isinstance(puzzle_bytes, str):
            with open(puzzle_bytes, 'rb') as f:
                puzzle_bytes = f.read()
        elif isinstance(puzzle_bytes, Puzzle):
            encoder = Encoder()
            puzzle_bytes.encode(encoder)
            puzzle_bytes = encoder.bytes
        solution_bytes = solution
        if isinstance(solution_bytes, str):
            with open(solution_bytes, 'rb') as f:
                solution_bytes = f.read()
        elif isinstance(solution_bytes, Solution):
            encoder = Encoder()
            solution_bytes.encode(encoder)
            solution_bytes = encoder.bytes
        self.verifier = ctypes.c_void_p(
            Sim.libverify().verifier_create_from_bytes(
                ctypes.c_char_p(bytes(puzzle_bytes)),
                ctypes.c_int(len(puzzle_bytes)),
                ctypes.c_char_p(bytes(solution_bytes)),
                ctypes.c_int(len(solution_bytes))
            )
        )
        if Sim.libverify().verifier_error(self.verifier):
            raise SimError(Sim.libverify(), self.verifier)

    def metric(self, metric):
        result = Sim.libverify().verifier_evaluate_metric(
            self.verifier,
            ctypes.c_char_p(metric.encode('utf-8'))
        )
        if Sim.libverify().verifier_error(self.verifier):
            raise SimError(Sim.libverify(), self.verifier)
        return result

    def approximate_metric(self, metric):
        result = Sim.libverify().verifier_evaluate_approximate_metric(
            self.verifier, ctypes.c_char_p(metric.encode('utf-8'))
        )
        if Sim.libverify().verifier_error(self.verifier):
            raise SimError(Sim.libverify(), self.verifier)
        return result

    def rate(self):
        cycles = self.metric('per repetition cycles')
        outputs = self.metric('per repetition outputs')
        if outputs == 0:
            return float('inf')
        return Fraction(cycles, outputs)

    def area_at_infinity(self):
        outputs = self.metric('per repetition outputs')
        a2 = self.approximate_metric('per repetition^2 area')
        if a2 != 0:
            if outputs == 0:
                return (float('inf'), float('inf'))
            else:
                return (2, a2 / (outputs * outputs))
        a1 = self.metric('per repetition area')
        if a1 != 0:
            if outputs == 0:
                return (float('inf'), float('inf'))
            else:
                return (1, Fraction(a1, outputs))
        return (0, self.metric('steady state area'))

    def output_intervals(self):
        return OutputIntervals.from_verifier(Sim.libverify(), self.verifier)

    def default_restrictions(self):
        overlap = self.metric("overlap")
        parts_of_type_baron = self.metric("parts of type baron")
        parts_of_type_glyph_disposal = self.metric(
            "parts of type "
            "glyph-disposal"
        )
        duplicate_reagents = self.metric("duplicate reagents")
        duplicate_products = self.metric("duplicate products")
        maximum_track_gap = self.metric("maximum track gap^2")

        return (
                overlap
                + max(0, parts_of_type_baron - 1)
                + max(0, parts_of_type_glyph_disposal - 1)
                + duplicate_reagents
                + duplicate_products
                + max(0, maximum_track_gap - 1)
        )


class OutputIntervals:

    def __init__(self, pattern=''):
        m = re.fullmatch(
            r"(\d+(?: \d+)*)|(?:(\d+(?: \d+)*) )?(?:\[(\d+(?: \d+)*)\])?",
            pattern
        )
        if m is None:
            raise ValueError('invalid syntax for output intervals')
        self.intervals = []
        if m.group(1) is not None:
            self.intervals.extend(map(int, m.group(1).split(' ')))
        if m.group(2) is not None:
            self.intervals.extend(map(int, m.group(2).split(' ')))
        self.repeats_after = len(self.intervals)
        if m.group(3) is not None:
            self.intervals.extend(map(int, m.group(3).split(' ')))
            n = len(self.intervals) - self.repeats_after
            for i in range(1, n // 2 + 1):
                if n % i != 0:
                    continue
                for j in range(n):
                    if self.intervals[self.repeats_after + j] != \
                            self.intervals[self.repeats_after + j % i]:
                        break
                else:
                    self.intervals = self.intervals[:self.repeats_after + i]
                    break
            while len(self.intervals) > 0 and self.repeats_after > 0 and \
                    self.intervals[-1] == self.intervals[
                self.repeats_after - 1]:
                self.repeats_after -= 1
                self.intervals = self.intervals[:-1]

    @classmethod
    def from_verifier(cls, libverify, verifier):
        intervals = cls()
        n = libverify.verifier_number_of_output_intervals(verifier)
        intervals.repeats_after = libverify.verifier_output_intervals_repeat_after(
            verifier
        )
        if intervals.repeats_after < 0:
            intervals.repeats_after = n
        intervals.intervals = []
        for i in range(n):
            intervals.intervals.append(
                libverify.verifier_output_interval(verifier, ctypes.c_int(i))
            )
        if libverify.verifier_error(verifier):
            raise SimError(libverify, verifier)
        return intervals

    def __str__(self):
        initial = ' '.join(map(str, self.intervals[:self.repeats_after]))
        if self.repeats_after == len(self.intervals):
            return initial
        repeating = ' '.join(map(str, self.intervals[self.repeats_after:]))
        if self.repeats_after == 0:
            return f'[{repeating}]'
        return f'{initial} [{repeating}]'

    def __repr__(self):
        return f'OutputIntervals("{str(self)}")'

    def __eq__(self, other):
        return self.intervals == other.intervals and self.repeats_after == other.repeats_after

    def compare_lt(self, other, *, allow_equality=True):
        i, j = 0, 0
        from math import lcm
        repeats = lcm(
            len(self.intervals) - self.repeats_after,
            len(other.intervals) - other.repeats_after
        )
        for iteration in range(
                max(self.repeats_after, other.repeats_after) + repeats
        ):
            if i == len(self.intervals):
                i = self.repeats_after
                if i == len(self.intervals):
                    return False
            if j == len(other.intervals):
                j = other.repeats_after
                if j == len(other.intervals):
                    return True
            if self.intervals[i] < other.intervals[j]:
                return True
            if self.intervals[i] > other.intervals[j]:
                return False
            i, j = i + 1, j + 1
        if self.repeats_after < other.repeats_after:
            return True
        if self.repeats_after > other.repeats_after:
            return False
        return allow_equality

    def __lt__(self, other):
        return self.compare_lt(other, allow_equality=False)

    def __gt__(self, other):
        return other.compare_lt(self, allow_equality=False)

    def __le__(self, other):
        return self.compare_lt(other, allow_equality=True)

    def __ge__(self, other):
        return other.compare_lt(self, allow_equality=True)


class SimError(Exception):

    def __init__(self, libverify, verifier):
        self.message = libverify.verifier_error(verifier).decode('utf-8')
        self.cycle = libverify.verifier_error_cycle(verifier)
        self.location = (libverify.verifier_error_location_u(verifier),
        libverify.verifier_error_location_v(verifier))
        libverify.verifier_error_clear(verifier)
        super().__init__(self.message, self.cycle, self.location)
