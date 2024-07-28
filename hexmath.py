import math

import numpy as np

from sharedtypes import *

rot_cw = np.array([[0, 1], [-1, 1]])


def _rotate_matrix(count):
    if count == 0:
        return np.identity(2).astype(int)
    else:
        return np.matmul(rot_cw, _rotate_matrix(count - 1))


_rotation_matrices = [_rotate_matrix(i) for i in range(6)]
ROTATION_VECTORS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]


def translate(position: Pos2D, translate_pos: Pos2D,
        translate_rotate: Angle) -> Pos2D:
    rotated = rotate(position, translate_rotate)
    return rotated[0] + translate_pos[0], rotated[1] + translate_pos[1]


def rotate(pos: Pos2D | Pos3D, rot_cw_times_60_deg: Angle,
        pivot: Pos2D = (0, 0)) -> Pos2D | Pos3D:
    if not (2 <= len(pos) <= 3):
        raise ValueError("pos must be 2 or 3 dimensional")
    if not (2 <= len(pivot) <= 3):
        raise ValueError("pos_axle must be 2 or 3 dimensional")
    rot_cw_times_60_deg %= 6

    # remember to return matching dims for pos
    return3d = (len(pos) == 3)

    # force args into 3d
    pos3 = cube_pos(pos)
    pivot3 = cube_pos(pivot)

    result = difference(pos3, pivot3)

    rot_m3 = rot_cw_times_60_deg % 3
    if rot_m3 == 2:
        result = (result[2], result[0], result[1])
    elif rot_m3 == 1:
        result = (result[1], result[2], result[0])

    if rot_cw_times_60_deg % 2 == 1:
        result = difference((0, 0, 0), result)

    result = summate(result, pivot3)

    if return3d:
        return result
    else:
        return result[0], result[1]


def cube_pos(pos: Pos2D | Pos3D) -> Pos3D:
    if len(pos) == 2:
        return pos[0], pos[1], -sum(pos)
    elif len(pos) == 3:
        return pos
    raise ValueError("pos must be a tuple of length 2 or 3")


def cab_distance(pos1: Pos2D, pos2: Pos2D = None) -> int:
    diff = cube_pos(pos1)
    if pos2 is not None:
        diff = np.subtract(diff, cube_pos(pos2))
    return np.abs(diff).sum().item() // 2


def direct_distance(pos1: Pos2D, pos2: Pos2D = None) -> float:
    diff = cube_pos(pos1)
    if pos2 is not None:
        diff = np.subtract(diff, cube_pos(pos2))
    return np.sqrt(np.pow(diff, 2).sum() / 2).item()


def difference(a: Pos, b: Pos) -> Pos:
    if len(a) == 2:
        return a[0] - b[0], a[1] - b[1]
    elif len(a) == 3:
        return a[0] - b[0], a[1] - b[1], a[2] - b[2]
    else:
        raise ValueError("a and b must be 2 or 3 dimensional")


def summate(a: Pos, b: Pos) -> Pos:
    if len(a) == 2:
        return a[0] + b[0], a[1] + b[1]
    elif len(a) == 3:
        return a[0] + b[0], a[1] + b[1], a[2] + b[2]
    else:
        raise ValueError("a and b must be 2 or 3 dimensional")


def scale(vec: Pos, scalar: float | int) -> Pos:
    if len(vec) == 2:
        return vec[0] * scalar, vec[1] * scalar
    elif len(vec) == 3:
        return vec[0] * scalar, vec[1] * scalar, vec[2] * scalar
    else:
        raise ValueError("a and b must be 2 or 3 dimensional")


def count_rotates_between(a: Pos2D, b: Pos2D, pivot: Pos2D = (0, 0)) -> int:
    a3 = cube_pos(a)
    b3 = cube_pos(b)

    # equal positions means no rotation, even if axle pos is the same pos
    if a3 == b3:
        return 0

    pivot3 = cube_pos(pivot)
    vec1 = difference(a3, pivot3)
    vec2 = difference(b3, pivot3)
    # confirm that radii are the same
    radius1 = direct_distance(vec1)
    radius2 = direct_distance(vec2)
    if not (radius1 == radius2):
        raise ValueError("pos1 and pos2 are not equal distances from pos_axle")

    # opposite vectors means 180 degrees
    if (np.subtract(vec1, vec2) == 0).all():
        return 3

    # determine left or right polarity for the rotation
    # cross should not be 0 since we tested for opposing and aligned vecs
    # positive means CCW, negative means CW
    cross = np.cross(vec1, vec2).item()

    # if the diff between pos1 and pos2 is the radius, then abs(rotate) is 1
    if direct_distance(a3, b3) == radius1:
        return int(math.copysign(1, cross))

    # the only valid possibility left to test is 2 rotations + or -.
    rotation_test_value = int(math.copysign(2, cross))
    rotated_2 = rotate(a3, rotation_test_value, pivot3)
    if np.equal(rotated_2, b3).all():
        return rotation_test_value

    # no valid rotation exists
    raise ValueError("pos1 cannot be rotated to pos2")
