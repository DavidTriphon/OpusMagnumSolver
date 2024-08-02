import math

import numpy as np

from sharedtypes import *

ROTATION_VECTORS = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]


def pos_from_direction(magnitude: int, direction: int | float,
        offset: int = None) -> Pos:
    if magnitude == 0:
        return 0, 0
    if magnitude < 0:
        raise ValueError("magnitude must be non-negative")
    if offset is None:
        offset = (direction % 1) * magnitude
    # forgive a certain amount of error
    if abs(round(offset) - offset) > 1e-14 * magnitude:
        raise ValueError(
            "the decimal component of rotation must be a valid integer"
            "fraction where the denominator is equal to magnitude")
    offset = round(offset)
    opposite = direction >= 3
    direction = int(direction) % 3
    if direction == 0:
        x, y = magnitude - offset, offset
    elif direction == 1:
        x, y = -offset, magnitude
    else:  # rotation == 2
        x, y = -magnitude, magnitude - offset

    if opposite:
        x, y = -x, -y

    return x, y


def direction_int(vec: Pos) -> int | None:
    # center
    if vec == (0, 0):
        return None

    # horizontal
    if vec[1] == 0:
        if vec[0] > 0:
            return 0
        else:
            return 3
    # vertical
    elif vec[0] == 0:
        if vec[1] > 0:
            return 1
        else:
            return 4
    # sum 0 diagonal
    elif vec[0] == -vec[1]:
        if vec[1] > 0:
            return 2
        else:
            return 5
    # fractional values
    return None


def direction_float(vec: Pos) -> float | None:
    # center
    if vec == (0, 0):
        return None

    x, y = vec

    vec_prod = x * y
    vec_sum = x + y
    vec_diff = x - y

    # no fraction on the dimension lines
    if vec_prod == 0 or vec_sum == 0:
        return direction_int(vec)

    # rotation 0 and 3, upper right and lower left
    if vec_prod > 0:
        base = 0 if y > 0 else 3
        return base + (y / vec_sum)

    # top and bottom diagonal quadrants
    if (vec_sum < 0) != (vec_diff < 0):  # vec_sum * vec_diff < 0
        base = 1 if y > 0 else 4
        return base + (vec_sum / y)

    # left and right diagonal quadrants
    else:  # vec_sum * vec_diff > 0
        base = 2 if y > 0 else 5
        return base + (vec_sum / x)


def angle_directions(rot1: float | int, rot2: float | int) -> float | int:
    return abs((rot1 - rot2 + 3) % 6 - 3)


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
