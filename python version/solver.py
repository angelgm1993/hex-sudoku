import numpy as np
import time
import random

start_time = time.clock()
random.seed(18)  # 18>7>14>22>19>30>31
T = .3112  # .3112 16.6sec

puzzle = np.empty((16, 16), dtype=int)
original_puzzle = np.empty((16, 16), dtype=str)
boxlims = [(0, 4), (4, 8), (8, 12), (12, 16)]
to_solve = []


def read_puzzle(puzzletxt):

    with open(puzzletxt) as f:
        temp = f.read()
        for puzzlerow, row in zip(temp.splitlines(), range(16)):
            for char, col in zip(puzzlerow, range(16)):
                if char in ['*', 'x']:
                    puzzle[row][col] = 16
                    original_puzzle[row][col] = char
                    to_solve.append([row, col])
                else:
                    puzzle[row][col] = int(char, 16)
                    original_puzzle[row][col] = char
    print_puzzle(puzzle, 0)
    print('\n')
    for x1 in boxlims:
        for y1 in boxlims:
            box = puzzle[slice(*x1), slice(*y1)]
            populate_box(box)


def print_puzzle(sudoku, solution):

    print(62 * '=')

    for i in range(16):
        for j in range(16):
            if j % 4 == 0:
                print('|| ', end='')

            print(f'{hex(sudoku[i, j])[2:]: <3}', end="") if solution else \
                print(f'{original_puzzle[i, j]: <3}', end="")

            if j == 15:
                print('||')

        if i % 4 == 3:
            print(62 * '=')


def populate_box(box):
    """
    Populates a Sudoku box with missing numbers.
    TODO safeguard against user feeding a puzzle with dupes in a box, row, or col.
    """
    missing = []
    [missing.append(z) for z in range(16) if z not in np.reshape(box, 16)]
    for i in range(4):
        for j in range(4):
            if box[i, j] == 16:
                box[i, j] = missing.pop()


def calc_energy(arr):
    e = 0
    for i, j in zip(range(16), range(16)):
        e += 16 - len(set(arr[i, :]))
        e += 16 - len(set(arr[:, j]))

    return e


def same_box_pair():
    pair1, pair2 = 0, 0
    check = True
    while check:
        pair1, pair2 = random.sample(to_solve, 2)
        box_i1 = [int(pair1[0] / 4) * 4, int(pair1[0] / 4) * 4 + 4]
        box_j1 = [int(pair1[1] / 4) * 4, int(pair1[1] / 4) * 4 + 4]
        box_i2 = [int(pair2[0] / 4) * 4, int(pair2[0] / 4) * 4 + 4]
        box_j2 = [int(pair2[1] / 4) * 4, int(pair2[1] / 4) * 4 + 4]
        check = not (box_i1 == box_i2 and box_j1 == box_j2)
    return pair1, pair2


def create_neighbor():
    """
    Randomly changes cells within the same box in puzzle if cells were initially blank.

    """
    copy = np.copy(puzzle)
    pair1, pair2 = same_box_pair()

    temp = copy[pair1[0], pair1[1]]
    copy[pair1[0], pair1[1]] = copy[pair2[0], pair2[1]]
    copy[pair2[0], pair2[1]] = temp

    return copy


read_puzzle(str(input("Enter name: ")) + ".txt")
# puzzle a took 10 sec
# puzzle b took 323 sec
# puzzle c took ? sec
# puzzle d took 233 sec
# puzzle e took ? sec

energy = calc_energy(puzzle)
accept_transition = False

while energy != 0:
    puzzle_neighbor = create_neighbor()
    x = calc_energy(puzzle_neighbor)
    deltaE = x - energy

    if deltaE < 0:
        accept_transition = True

    else:
        rand_val = random.random()
        prob = np.exp(-deltaE / T)
        accept_transition = (rand_val < prob)

    if accept_transition:
        energy = x
        puzzle = np.copy(puzzle_neighbor)

print_puzzle(puzzle, 1)

print(f"--- {(time.clock() - start_time): 0.2f} seconds ---")
