import numpy as np
import os
import find_terminal
from structure_mask import StructureMask
from constants import *
import pandas as pd


def correct_point(point):  # correct some points which is zero in annotation
    bias = [-2, -1, 0, 1, 2]
    points = []
    '''
    for x in bias:
        points.append([point[0] + x, point[1], point[2]])
    for y in bias:
        points.append([point[0], point[1] + y, point[2]])
    for z in bias:
        points.append([point[0], point[1], point[2] + z])
    '''
    for x in bias:
        for y in bias:
            for z in bias:
                points.append([point[0] + x, point[1] + y, point[2] + z])
    # calculate distance between point ans its neighbour points
    d = np.zeros([len(points)])
    num = 0
    for n in points:
        chazhi = np.abs(np.array(n) - np.array(point))
        distance = np.sum(chazhi ** 2, axis=0, keepdims=True) ** 0.5
        d[num] = distance
        num = num + 1
    # get region id label for every point
    points_label = [annotation[tuple(i)] for i in points]
    array = np.array(points_label)
    # find the nearest neighbour point and get the region label
    idx_nonzero = np.where(array > 0)[0]  # find points which region label are not zero
    mask_data = array.copy()
    mask_data[idx_nonzero] = 1
    mask_data = mask_data * d
    min_point_d = min(filter(lambda e: e > 0, list(mask_data)))
    idx = np.where(mask_data == min_point_d)
    point_label = array[idx[0][0]]
    return point_label


def find_terminal_belong_to_one_hemisphere(point):
    point = np.round(np.array(point))
    point = np.ceil(np.array(point) / 100)
    if point[VERTICAL_PLANE_IDX] >= VOXEL_SHAPE[VERTICAL_PLANE_IDX] // 2:
        hemisphere = R_HEMISPHERE
    else:
        hemisphere = L_HEMISPHERE
    return hemisphere


def find_terminal_belong_to_one_region(point):
    point = np.round(np.array(point))
    point = np.ceil(np.array(point) / 100)
    point = [int(point[0]), int(point[1]), int(point[2])]
    point = tuple(point)
    point_id = annotation[point]
    if point_id == 0:
        print('000000000')
        # parent_name = 'zero'
        correct_point_id = correct_point(point)
        parent = structure_mask.structure_tree.ancestor_ids([correct_point_id])[0][0]
        parent_name = structure_mask.structure_tree.get_structures_by_id([parent])[0]['acronym']
    else:
        parent = structure_mask.structure_tree.ancestor_ids([point_id])[0][0]
        parent_name = structure_mask.structure_tree.get_structures_by_id([parent])[0]['acronym']

    return parent_name


def get_neuron_projection_label(soma_hemis, terminal_hemis, terminal_numbers):
    terminal_hemis = np.array(terminal_hemis)
    idx_num = len(np.where(terminal_hemis == soma_hemis)[0])
    if idx_num == 0:
        label = 2  # only contra projection
    elif idx_num == terminal_numbers:
        label = 1  # only ipsi projection
    else:
        label = 3  # ipsi and contra projection all exist
    return label


if __name__ == "__main__":
    data_path = '/data/single_neuro_tracing/tracing_to_allen_dendrites_new'
    structure_mask = StructureMask()
    annotation = structure_mask.annotations
    neuron = []
    dendrites_length = []
    branch_number = []

    files = os.listdir(data_path)
    files.sort()
    for file in files:
        print(file)
        file_path = ''.join([data_path, '/', file])
        data = find_terminal.read_swc(file_path)
        data = find_terminal.inspect_data(data)
        data_copy = data
        data_copy = data_copy.values
        point_type = data_copy[:, 1]
        if 3 not in point_type:
            continue

        dendrite_branch, dendrite_max_length = find_terminal.find_branch_number_and_max_length_for_dendrite(data)
        dendrites_length.append(dendrite_max_length)
        branch_number.append(dendrite_branch)
        neuron.append(file[0:-4])

    summary = {'neuron': neuron, 'max_length': dendrites_length, 'branch_number': branch_number}
    df = pd.DataFrame(summary)
    df.to_csv('./IT_neuron_projection_dendrites_branch_and_max_length.csv')
