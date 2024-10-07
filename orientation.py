import hexmath
import om


class PartOrientation:

    def __init__(self, part, pos, rot, hexes):
        self.part = part
        self.pos = pos
        self.rot = rot
        self.hexes = hexes

    def __hash__(self):
        return hash((self.part, self.pos, self.rot, tuple(self.hexes)))

    def __eq__(self, other):
        if not isinstance(other, PartOrientation):
            return False
        return (self.part == other.part and self.pos == other.pos and self.rot
                == other.rot and self.hexes == other.hexes)

    def __str__(self):
        return "PartOrientation(part=\"%s\", pos=%r, rot=%d, hexes=%s)" % (
            self.part, self.pos, self.rot, self.hexes)


def adjacent_bonder_orientations():
    return {
        PartOrientation(
            "bonder",
            hexmath.pos_from_direction(1, dir),
            (dir + 2) % 6,
            {hexmath.pos_from_direction(1, dir),
                hexmath.pos_from_direction(1, dir + 1)}
        )
        for dir in range(6)
    }


def find_bonder_orientations(mol, bond_index, arm_length=1, arm_pos=(0, 0),
        part_key="b"):
    mol = mol.copy()
    bond = mol.bonds.pop(bond_index)
    mols = mol.find_disconnected()

    grab_positions = [
        hexmath.summate(hexmath.pos_from_direction(arm_length, dir), arm_pos)
        for dir in range(6)
    ]

    bond_atom_is = tuple(
        i
            for i, a in enumerate(mol.atoms)
            if a.position in bond.positions
    )
    assert len(bond_atom_is) == 2

    all_bonder_orientations = set()

    if len(mols) == 1:
        for rotation in range(6):
            rotated_mol = mol.rotate(rotation)
            for grab_pos in grab_positions:
                for atom_grabbed_i in range(len(mol.atoms)):
                    moved_mol = rotated_mol.move_atom_to(grab_pos,
                        atom_grabbed_i)
                    # skip if overlaps arm
                    if all(a.position != arm_pos for a in moved_mol.atoms):
                        # determine where the bonder would be for this
                        positions = tuple(sorted(moved_mol.atoms[i].position
                            for i in bond_atom_is))
                        direction = hexmath.direction_int(
                            hexmath.difference(positions[2], positions[1]))
                        all_bonder_orientations.add(PartOrientation(
                            part=part_key,
                            pos=positions[0], rot=direction,
                            hexes=set(positions)
                        ))
    elif len(mols) == 2:
        # determine the indexes of the two groups:
        mol0_positions = {atom.position for atom in mols[0].atoms}
        mol1_positions = {atom.position for atom in mols[1].atoms}
        mol0_is = [
            i
            for i, a in enumerate(mol.atoms)
            if a.position in mol0_positions
        ]
        mol1_is = [
            i
            for i, a in enumerate(mol.atoms)
            if a.position in mol1_positions
        ]

        smaller_is = mol0_is if len(mol0_is) < len(mol1_is) else mol1_is
        larger_is = mol1_is if len(mol0_is) < len(mol1_is) else mol0_is

        for rotation in range(6):
            rotated_mol = mol.copy().rotate(rotation)
            for grab_pos in grab_positions:
                for atom_grabbed_i in smaller_is:
                    moved_mol = rotated_mol.copy().move_atom_to(grab_pos,
                        atom_grabbed_i)
                    # skip if overlaps arm
                    overlaps_arm = any(a.position == arm_pos for a in
                        moved_mol.atoms)
                    # skip if other side of bonding is not in grabbable space
                    grabbable = any(
                        moved_mol.atoms[i].position in grab_positions
                            for i in larger_is)
                    if not overlaps_arm and grabbable:
                        # determine where the bonder would be for this
                        positions = tuple(sorted(moved_mol.atoms[i].position
                            for i in bond_atom_is))
                        direction = hexmath.direction_int(
                            hexmath.difference(positions[1], positions[0]))
                        all_bonder_orientations.add(PartOrientation(
                            part=part_key,
                            pos=positions[0], rot=direction,
                            hexes=set(positions)
                        ))
    else:
        raise ValueError("Molecule contains other a number of disconnected "
                         "molecules other than 1 or 2 ???")

    return list(all_bonder_orientations)


def find_calcifier_orientations_for_mols(mol, atom_i):
    pass


def find_molecule_orientations(mol: om.Molecule, key, armlength: int = 1):
    orientations = set()
    for arm_rot in range(6):
        grab_pos = hexmath.pos_from_direction(armlength, arm_rot)
        for atom_grab_i in range(len(mol.atoms)):
            # only iterate the rotations if there are multiple atoms
            for part_rot in range(6 if (len(mol.atoms) > 1) else 1):
                rotated_product = mol.copy().move_atom_to(grab_pos, atom_grab_i
                ).rotate(part_rot, grab_pos)
                # skip orientations where they overlap the arm
                if all(atom.position != (0, 0) for atom in
                        rotated_product.atoms):
                    orientations.add(PartOrientation(
                        key,
                        rotated_product.atoms[0].position,
                        part_rot,
                        {atom.position for atom in rotated_product.atoms}
                    ))
    return orientations


class Constraints:

    def __init__(self, orientations):
        self.parts_orientations = {}
        self.hexes_orientations = {}
        # populate the constraint sets with every constraint that the values
        # overlap for
        for o in orientations:
            self.parts_orientations[o.part] = \
                self.parts_orientations.get(o.part, set()) | {o}
            for h in o.hexes:
                self.hexes_orientations[h] = \
                    self.hexes_orientations.get(h, {None}) | {o}
        # check for possibilities that are constrained and enforce them
        for constraint_set in list(self.parts_orientations.values()) + \
                              list(self.hexes_orientations.values()):
            if len(constraint_set) == 1:
                if not self.assign(next(iter(constraint_set))):
                    raise ValueError(
                        "orientation options implied a contradiction")
            elif len(constraint_set) == 0:
                raise ValueError("orientation options implied a contradiction")

    def copy(self):
        pass

    def is_constrained(self, part):
        if part not in self.parts_orientations:
            raise ValueError("part is not listed among constraints")
        orientation_set = self.parts_orientations[part]
        return len(orientation_set) == 1

    def get(self, part):
        if part not in self.parts_orientations:
            raise ValueError("part is not listed among constraints")
        orientation_set = self.parts_orientations[part]
        if len(orientation_set) == 1:
            return next(iter(orientation_set))
        else:
            return None

    def add_part(self, orientation_set):
        parts = {o.part for o in orientation_set}
        if len(parts) > 1:
            raise ValueError("All new orientations must be for the same part")
        part = next(iter(parts))
        if part in self.parts_orientations:
            raise ValueError("All new orientations must be for a NEW part")
        occupied_hexes = {
            h for h, os in self.hexes_orientations.items()
            if len(os) == 1
        }
        non_overlapping_orientations = {
            o for o in orientation_set
            if all(h not in occupied_hexes for h in o.hexes)
        }
        self.parts_orientations[part] = non_overlapping_orientations
        for o in non_overlapping_orientations:
            for h in o.hexes:
                self.hexes_orientations[h].add(o)
        if len(non_overlapping_orientations) == 1:
            self.assign(next(iter(non_overlapping_orientations)))

    def _assign_single(self, group, key, assigned):
        assert key in group
        constraint_set = group[key]
        assert assigned in constraint_set
        to_eliminate = {
            o for o in constraint_set
            if o != assigned and o is not None
        }
        group[key] = {assigned}
        # return orientations to eliminate
        return to_eliminate

    def assign(self, assigned):
        # remove other possible orientations for this part
        to_eliminate = self._assign_single(self.parts_orientations,
            assigned.part, assigned)
        # remove the orientations of other parts that overlap this part
        for hex in assigned.hexes:
            to_eliminate |= self._assign_single(self.hexes_orientations, hex,
                assigned)
        for o in to_eliminate:
            if not self.eliminate(o):
                return False
        return True

    def _elim_single(self, group, key, eliminated):
        assert key in group
        constraint_set = group[key]
        if eliminated in constraint_set:
            constraint_set.remove(eliminated)
            options_left = len(constraint_set)
            if options_left == 0:
                return False, None
            elif options_left == 1:
                last_option = next(iter(constraint_set))
                return True, last_option
        return True, None

    def eliminate(self, eliminated):
        # eliminate the orientation from all constraint sets
        # and propagate when a constraint set is reduced to one option
        # return false when no options left in a set
        to_assign = set()
        valid, new_assign = self._elim_single(self.parts_orientations,
            eliminated.part, eliminated)
        if not valid:
            return False
        to_assign.add(new_assign)
        for h in eliminated.hexes:
            valid, new_assign = self._elim_single(self.hexes_orientations,
                h, eliminated)
            if not valid:
                return False
            to_assign.add(new_assign)
        # iterate over new assertions
        for o in to_assign:
            if o is not None:
                self.assign(o)
        # no previous fail means success
        return True
