def moves_from_steps(step_sequence: list):
    moves = []
    positions = {}
    attached = {}
    occupied = {}

    def add_move(target, dest):
        moves.append((target, positions[target], dest))
        positions[target] = dest
        for nb in attached.setdefault(target, []):
            positions[nb] = '-'

    for step in step_sequence:
        type = step[0]
        if type == 'r':
            target0 = int(step[1])
            assert target0 not in positions
            positions[target0] = 'r'
        elif type == 'c':
            target0 = int(step[1])
            assert target0 in positions
            add_move(target0, 'c')
        elif type == 'b':
            target0 = int(step[1])
            target1 = int(step[2])
            assert target0 in positions
            assert target1 in positions
            add_move(target0, 'b')
            add_move(target1, 'b')
            all_att = attached[target0] + [target0] + attached[target1] + [
                target1]
            for a in all_att:
                attached[a] = all_att.copy()
                attached[a].remove(a)
        elif type == 'p':
            add_move(0, 'p')
            attached = {}
            occupied = {}
            positions = {}
    return moves
