import numpy as np
import pandas as pd
import re
from get_neuron_terminal_depth_and_branch_degree_another_version import correct_point
from structure_mask import StructureMask
import scipy.io as sio
from matplotlib import pyplot as plt
from scipy.stats import f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def read_isocortex_region():
    region_path = '/data/gene_preprocess/projection_model/AllenMouseProjectionModel/region_information/cortex_region_names_information.csv'
    regions_name = pd.read_csv(region_path, usecols=[1])
    regions_name = regions_name.values
    regions_names = []
    for i in regions_name:
        regions_names.extend(i)

    return regions_names


def remove_some_information(neuron_termials):
    terminals_regions = []
    for k in neuron_termials:  # k is corresponding to neuron
        if k == 'set()':
            terminals_region = []
            terminals_regions.append(terminals_region)
            continue
        terminals = str.split(k, '\'')
        terminals_region = [terminal for terminal in terminals if
                            terminal != ', ' and terminal != '{' and terminal != '}']
        terminals_regions.append(terminals_region)
    return terminals_regions


def remove_layer_information(region):
    if '-' in region:
        rex_compile = re.compile("([A-Z]*)([a-z]*)(-*)([a-z]*)([0-9]*)")
        rex = rex_compile.search(region)

        if rex.group(1) is '':
            rex = region
        else:
            rex = rex.group(1) + rex.group(2) + rex.group(3) + rex.group(4)
    else:
        rex_compile = re.compile("([A-Z]*)([a-z]*)([0-9]*)")
        rex = rex_compile.search(region)

        if rex.group(1) is '':
            rex = region
        else:
            rex = rex.group(1) + rex.group(2)

    return rex


def find_5_type_neuron_idx(ipsi_termials, contra_terminals):
    # get the ratio of neuron project to ipsi, contra, both respectively
    s1 = []
    s2 = []
    s3 = []
    for i in range(len(ipsi_termials)):
        ipsi_terminal = ipsi_termials[i]
        contra_terminal = contra_terminals[i]
        ipsi_num1 = len(ipsi_terminal)
        contra_num2 = len(contra_terminal)

        res = list(set(ipsi_terminal) & set(contra_terminal))
        uion = list(set(ipsi_terminal).union(contra_terminal))
        diff1 = list(set(ipsi_terminal).difference(contra_terminal))
        diff2 = list(set(contra_terminal).difference(ipsi_terminal))
        ratio1 = len(diff1) / len(uion)  # region that the neuron only projects to ipsi
        ratio2 = len(diff2) / len(uion)  # region that the neuron only projects to contra
        ratio3 = len(res) / len(uion)  # region that the neuron both projects to ipsi and contra
        if round(ratio1 + ratio2 + ratio3, 3) != 1:
            raise Exception
        s1.append(ratio1)
        s2.append(ratio2)
        s3.append(ratio3)

    c1 = []  # least one term == 0
    c11 = []  # ipsi != 0, contra == 0, both != 0
    c111 = []
    c112 = []
    c12 = []  # ipsi == 0, contra != 0, both != 0
    c121 = []
    c122 = []
    c13 = []  # ipsi == 0, contra == 0, both != 0
    c14 = []  # ipsi != 0, contra != 0, both == 0
    c141 = []  # ipsi >= contra
    c142 = []  # ipsi < contra
    c2 = []  # only ipsi, only contra, both all != 0
    for j in range(len(ipsi_termials)):
        r1 = s1[j]  # get the ratio of single neuron in three types
        r2 = s2[j]
        r3 = s3[j]
        if r1 != 0 and r2 != 0 and r3 != 0:
            c2.append(j)
        else:
            c1.append(j)
            if r1 != 0 and r2 == 0 and r3 != 0:
                c11.append(j)
                if r1 >= r3:
                    c111.append(j)
                else:
                    c112.append(j)
            elif r1 == 0 and r2 != 0 and r3 != 0:
                c12.append(j)
                if r2 >= r3:
                    c121.append(j)
                else:
                    c122.append(j)
            elif r1 == 0 and r2 == 0 and r3 != 0:
                c13.append(j)
            elif r1 != 0 and r2 != 0 and r3 == 0:
                c14.append(j)
                if r1 >= r2:
                    c141.append(j)
                else:
                    c142.append(j)
    return c1, c11, c111, c112, c12, c121, c122, c13, c14, c141, c142, c2


def get_cortex_regions_from_terminals(neuron_terminals, cortex_region):
    terminals_cortex_regions = []
    for k in neuron_terminals:  # k is corresponding to neuron
        if k == '[]':
            terminals_cortex = []
            terminals_cortex_regions.append(terminals_cortex)
            continue
        terminals = str.split(k, '\'')
        terminals_cortex = [terminal for terminal in terminals if
                            terminal != ', ' and terminal != '[' and terminal != ']' and remove_layer_information(
                                terminal) in cortex_region]
        terminals_cortex_regions.append(terminals_cortex)

    return terminals_cortex_regions


def get_terminals_length_list(neuron_terminals):
    terminals_cortex_regions = []
    for k in neuron_terminals:  # k is corresponding to neuron
        if k == '[]':
            terminals_cortex = []
            terminals_cortex_regions.append(terminals_cortex)
            continue
        terminals = str.split(k[1:-1], ', ')
        terminals_cortex = [float(terminal) for terminal in terminals if
                            terminal != ', ' and terminal != '[' and terminal != ']']
        terminals_cortex_regions.append(terminals_cortex)

    return terminals_cortex_regions


def draw_error_bar2(x):
    from scipy import stats
    x_mean = []
    x_sem = []
    for i in x:
        i_mean = np.mean(i)
        i_sem = float(stats.sem(i))
        x_mean.extend([i_mean])
        x_sem.extend([i_sem])
    return x_mean, x_sem

def draw_error_bar_median(x):
    from scipy import stats
    x_mean = []
    x_sem = []
    for i in x:
        i_mean = np.median(i)
        i_sem = stats.sem(i)
        x_mean.extend([i_mean])
        x_sem.extend([i_sem])
    return x_mean, x_sem


if __name__ == "__main__":
    ## data1_path: terminal length of full terminals:3744 neurons
    data1_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_neuron_projection_full_summary_cortex_ceil_av.csv'
    data1 = pd.read_csv(data1_path, usecols=[1, 4, 5, 6, 7])
    data1 = data1.values
    neuron_id = data1[:, 0]
    ipsi_terminal_length = data1[:, 1]
    ipsi_terminal_max_length = data1[:, 2]
    contra_terminal_length = data1[:, 3]
    contra_terminal_max_length = data1[:, 4]
    ## data2_path: terminal's region name of full terminals:3744 neurons
    data2_path = '/data/single_neuro_tracing/code/preprocess/IT_neuron_projection_summary_ceil_nounique.csv'
    data2 = pd.read_csv(data2_path, usecols=[1, 2, 4, 5])
    data2 = data2.values
    neuron_id1 = data2[:, 0]
    ipsi_terminal_region = data2[:, 2]
    contra_terminal_region = data2[:, 3]
    ## data3_path: cortex neuron projected bilaterally to cortex regions:3004 neurons
    data3_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data3 = pd.read_csv(data3_path, usecols=[0, 1, 2, 3, 4])
    data3 = data3.values
    neuron_id2 = data3[:, 1]
    soma_region = data3[:, 2]
    ipsi_cortex_terminals = data3[:, 3]
    contra_cortex_terminals = data3[:, 4]
    ipsi_cortex_terminals = remove_some_information(ipsi_cortex_terminals)
    contra_cortex_terminals = remove_some_information(contra_cortex_terminals)
    idx = [j for i in range(len(neuron_id2)) for j in range(len(neuron_id)) if neuron_id2[i] == neuron_id[j]]

    ipsi_terminal_region = ipsi_terminal_region[idx]
    ipsi_terminal_length = ipsi_terminal_length[idx]
    contra_terminal_region = contra_terminal_region[idx]
    contra_terminal_length = contra_terminal_length[idx]

    # read isocortex region
    cortex_regions = read_isocortex_region()
    # transform string to list
    ipsi_cortex_terminal_regions = get_cortex_regions_from_terminals(ipsi_terminal_region, cortex_regions)
    contra_cortex_terminal_regions = get_cortex_regions_from_terminals(contra_terminal_region, cortex_regions)
    ipsi_terminal_length_lists = get_terminals_length_list(ipsi_terminal_length)
    contra_terminal_length_lists = get_terminals_length_list(contra_terminal_length)

    # generate max terminal length of four type terminals[i, bi, bc, c] and the number of three regions that neuron
    # projected to [i, b, c]
    count = np.zeros((3004, 3))  # matrix: the number of three type regions that neuron projected to
    max_length = np.zeros((3004, 4))
    total_length = np.zeros((3004, 4))
    #
    # for i in range(len(soma_region)):  # 3004 neurons
    #     ipsi_cortex_terminal_region = ipsi_cortex_terminal_regions[i]
    #     ipsi_terminal_length_list = ipsi_terminal_length_lists[i]
    #     contra_cortex_terminal_region = contra_cortex_terminal_regions[i]
    #     contra_terminal_length_list = contra_terminal_length_lists[i]
    #     ipsi_terminal = ipsi_cortex_terminals[i]
    #     contra_terminal = contra_cortex_terminals[i]
    #     res = list(set(ipsi_terminal) & set(contra_terminal))  # region that the neuron both projects to ipsi and contra
    #     diff1 = list(set(ipsi_terminal).difference(contra_terminal))  # region that the neuron only projects to ipsi
    #     diff2 = list(set(contra_terminal).difference(ipsi_terminal))  # region that the neuron only projects to contra
    #
    #     count[i, 0] = len(diff1)
    #     count[i, 1] = len(res)
    #     count[i, 2] = len(diff2)
    #
    #     only_ipsi_terminal_length_remained = [ipsi_terminal_length_list[j] for j in
    #                                           range(len(ipsi_cortex_terminal_region))
    #                                           if remove_layer_information(ipsi_cortex_terminal_region[j]) in diff1]
    #     only_contra_terminal_length_remained = [contra_terminal_length_list[k] for k in
    #                                             range(len(contra_cortex_terminal_region))
    #                                             if remove_layer_information(contra_cortex_terminal_region[k]) in diff2]
    #     B_ipsi_terminal_length_remained = [ipsi_terminal_length_list[j] for j in
    #                                        range(len(ipsi_cortex_terminal_region))
    #                                        if remove_layer_information(ipsi_cortex_terminal_region[j]) in res]
    #     B_contra_terminal_length_remained = [contra_terminal_length_list[k] for k in
    #                                          range(len(contra_cortex_terminal_region))
    #                                          if remove_layer_information(contra_cortex_terminal_region[k]) in res]
    #     if len(only_ipsi_terminal_length_remained) == 0:
    #         length1 = 0
    #         length1_total = 0
    #     else:
    #         length1 = np.max(np.array(only_ipsi_terminal_length_remained))
    #         length1_total = np.sum(np.array(only_ipsi_terminal_length_remained))
    #     if len(B_ipsi_terminal_length_remained) == 0:
    #         length2 = 0
    #         length2_total = 0
    #     else:
    #         length2 = np.max(np.array(B_ipsi_terminal_length_remained))
    #         length2_total = np.sum(np.array(B_ipsi_terminal_length_remained))
    #     if len(B_contra_terminal_length_remained) == 0:
    #         length3 = 0
    #         length3_total = 0
    #     else:
    #         length3 = np.max(np.array(B_contra_terminal_length_remained))
    #         length3_total = np.sum(np.array(B_contra_terminal_length_remained))
    #     if len(only_contra_terminal_length_remained) == 0:
    #         length4 = 0
    #         length4_total = 0
    #     else:
    #         length4 = np.max(np.array(only_contra_terminal_length_remained))
    #         length4_total = np.sum(np.array(only_contra_terminal_length_remained))
    #     max_length[i, 0] = length1
    #     max_length[i, 1] = length2
    #     max_length[i, 2] = length3
    #     max_length[i, 3] = length4
    #
    #     total_length[i, 0] = length1_total
    #     total_length[i, 1] = length2_total
    #     total_length[i, 2] = length3_total
    #     total_length[i, 3] = length4_total
    #
    # sio.savemat('./new_results/matrix/total_terminal_length_four_type.mat', {'four_type': total_length})
    # sio.savemat('./new_results/matrix/max_terminal_length_four_type.mat', {'four_type': max_length})
    # sio.savemat('./new_results/matrix/the_number_of_regions_projected_to_3004neurons.mat', {'three_count': count})

    max_terminal_length = sio.loadmat('./new_results/matrix/max_terminal_length_four_type.mat')
    numbers_region = sio.loadmat('./new_results/matrix/the_number_of_regions_projected_to_3004neurons.mat')
    max_terminal_length_sorted = max_terminal_length['four_type']
    numbers_region_sorted = numbers_region['three_count']

    [c1, c11, c111, c112, c12, c121, c122, c13, c14, c141, c142, c2] = find_5_type_neuron_idx(ipsi_cortex_terminals,
                                                                                              contra_cortex_terminals)

    # draw terminal length figure
    # max_terminal_length_sorted = max_terminal_length_sorted # only draw figure: if max: delete the words, if not: remained the words
    i_max_terminal_length_sorted_subtype = []
    bi_max_terminal_length_sorted_subtype = []
    bc_max_terminal_length_sorted_subtype = []
    c_max_terminal_length_sorted_subtype = []
    for i in [c11, c12, c13, c14, c2]:
        i_max_terminal_length_sorted_subtype.append(max_terminal_length_sorted[i][:, 0])
        bi_max_terminal_length_sorted_subtype.append(max_terminal_length_sorted[i][:, 1])
        bc_max_terminal_length_sorted_subtype.append(max_terminal_length_sorted[i][:, 2])
        c_max_terminal_length_sorted_subtype.append(max_terminal_length_sorted[i][:, 3])
    save_path = '/data/single_neuro_tracing/code/preprocess/new_results/figures/'
    subtype = ['c11', 'c12', 'c13', 'c14', 'c2']
    i_maxlength_mean, i_maxlength_sem = draw_error_bar2(i_max_terminal_length_sorted_subtype)
    bi_maxlength_mean, bi_maxlength_sem = draw_error_bar2(bi_max_terminal_length_sorted_subtype)
    bc_maxlength_mean, bc_maxlength_sem = draw_error_bar2(bc_max_terminal_length_sorted_subtype)
    c_maxlength_mean, c_maxlength_sem = draw_error_bar2(c_max_terminal_length_sorted_subtype)

    import seaborn
    seaborn.set_style('white')
    plt.figure()
    X = np.arange(0, 5, 1)
    matrix = np.array([i_maxlength_mean, bi_maxlength_mean, bc_maxlength_mean, c_maxlength_mean])
    matrix1 = np.array([i_maxlength_sem, bi_maxlength_sem, bc_maxlength_sem, c_maxlength_sem])
    width = 0
    wid = 0.2
    for i in range(4):
        plt.bar(X + width, matrix[i, :], width=wid)
        width = width + wid
    plt.xticks(X + width / 2, ['IB', 'BC', 'B', 'IC', 'IBC'])
    plt.legend(['si', 'bi', 'bc', 'sc'])
    width = 0
    wid = 0.2
    for i in range(4):
        for j in range(5):
            plt.plot([X[j] + width, X[j] + width], [matrix[i, j] - matrix1[i, j], matrix[i, j] + matrix1[i, j]], color=[0, 0, 0])
        width = width + wid
    plt.ylabel('mm')
    plt.xlabel('types')

    save_fn = save_path + 'max_terminal_length_in_subtype_refine_mean_sem_bar' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()
 
    # terminal_summary = []
    # types = []
    # m = 0
    # for i in bc_max_terminal_length_sorted_subtype:
    #     terminal_summary.extend(i)
    #     types.extend([subtype[m]] * len(i))
    #     m = m + 1
    # a = {'length': terminal_summary, 'type': types}
    # tukey_result = pairwise_tukeyhsd(a['length'], a['type'])
    # print(tukey_result.summary())

    i_terminal_region_count_subtype = []
    b_terminal_region_count_subtype = []
    c_terminal_region_count_subtype = []

    for i in [c11, c12, c13, c14, c2]:
        i_terminal_region_count_subtype.append(numbers_region_sorted[i][:, 0])
        b_terminal_region_count_subtype.append(numbers_region_sorted[i][:, 1])
        c_terminal_region_count_subtype.append(numbers_region_sorted[i][:, 2])

    save_path = '/data/single_neuro_tracing/code/preprocess/new_results/figures/'
    subtype = ['c11', 'c12', 'c13', 'c14', 'c2']
    i_mean, i_sem = draw_error_bar2(i_terminal_region_count_subtype)
    b_mean, b_sem = draw_error_bar2(b_terminal_region_count_subtype)
    c_mean, c_sem = draw_error_bar2(c_terminal_region_count_subtype)

    plt.figure()
    x1 = np.arange(1, 6, 1)
    x2 = x1 - 0.2
    x3 = x1 + 0.2
    plt.bar(x2, i_mean, width=0.2)
    # plt.scatter(x2, i_terminal_region_count_subtype)
    plt.bar(x1, b_mean, width=0.2)
    # plt.scatter(x1, b_terminal_region_count_subtype)
    plt.bar(x3, c_mean, width=0.2)
    # plt.scatter(x3, c_terminal_region_count_subtype)
    plt.xticks(x1, subtype)

    plt.legend(['i', 'b', 'c'])
    save_fn = save_path + 'three_region_counts_in_subtype_mean' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()
    plt.show()


    print(1)
