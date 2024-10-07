import orientation
import graph
from frame import Frame


class PathFindingNode(graph.BaseNode):

    def __init__(self, frame: Frame):
        self.frame = frame
        self.prev_moves = frozenset()
        self.next_bond = None
        self.positions = dict()
        self.atom_waiting = dict()
        self.part_constraints = orientation.Constraints()

    def get_mol_of_atom_i(self, atom_i):
        atom_pos = self.positions.get(atom_i)
        if atom_pos is not None:
            mols = self.frame.get_molecules()
            return next((
                m, i for m in mols
                for i, a in enumerate(m.atoms)
                if a.position == atom_pos
            ))
        else:
            mol = self.frame.puzzle.reagents[0].copy()
            reagent_position = self.part_constraints.get("r").pos
            return mol.translate(reagent_position)

    def brute_force_path(self, all_moves):
        # step 1: find all available moves,
        # and choose one
        available_moves = [
            m for m in all_moves
            if m not in self.prev_moves and
               self.atom_waiting.get(m[0]) is None and
               (m[1][0] != 'b' or self.next_bond is None)
        ]
        for move in available_moves:
            self._step2(move)

    def _step2(self, move):
        # step 2: choose an orientation for the part involved in the move
        part_indicator = move[1][0]
        target_atom_i = move[0]
        all_atom_is = [int(t) for t in move[1][1:]]
        if part_indicator not in self.part_constraints.parts_orientations:
            # figure out where the atom can reach for each atom in the step
            for atom_i in all_atom_is:
                mol, mol_atom_i = self.get_mol_of_atom_i(atom_i)
                mol_orientations = orientation.find_molecule_orientations(mol,
                    part_indicator)
                atom_positions = [
                    mol.copy().move_atom_to(o.pos).rotate(o.rot).atoms[
                        mol_atom_i].position
                    for o in mol_orientations
                ]


        # step 3: figure out possible destinations to satisfy the move,
        # and choose one

        # step 4: find a movement path that completes the "move" and
        # calculate its instruction sequence

        # step 5: return the next path finding node

    def __hash__(self):
        return hash((self.frame, self.prev_steps,
        self.next_bond, self.positions))

    def __eq__(self, other):
        if not isinstance(other, PathFindingNode):
            return False
        return (
                self.frame == other.frame and
                self.prev_steps == other.prev_steps and
                self.next_bond == other.next_bond and
                self.positions == other.positions
        )
