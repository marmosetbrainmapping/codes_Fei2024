import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt
import get_single_neuron_projection_ratio
import get_region_connection_for_each_layer
import seaborn
from scipy.stats import f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import seaborn


def find_5_type_neuron_idx(ipsi_termials, contra_terminals, soma_region):
    ipsi_termials = get_single_neuron_projection_ratio.remove_some_information(ipsi_termials)
    contra_terminals = get_single_neuron_projection_ratio.remove_some_information(contra_terminals)
    # get the ratio of neuron project to ipsi, contra, both respectively
    s1 = []
    s2 = []
    s3 = []
    for i in range(len(neuron_id)):
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
    for j in range(len(neuron_id)):
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


def get_number_list(a):
    # if a == '[]':
    #     return 0
    tlen = str.split(a[1:-1], ',')
    # print(len(tlen))
    # if len(tlen) == 1:
    #     print(a)
    #     print(tlen)
    terminals_length = [float(terminal) for terminal in tlen]
    return terminals_length


def get_number_list_1(a):
    tlen = str.split(a[1:-1], ' ')
    terminals_length = []
    for terminal in tlen:
        if '\n' in terminal:
            terminals_length.append(float(terminal[0:-1]))
        elif terminal == '':
            continue
        elif terminal == '...':
            continue
        else:
            terminals_length.append(float(terminal))
    return terminals_length


def get_neuron_number_in_different_region(soma_region, region_sets):
    neuron_idx = []
    for i in region_sets:
        c = []
        for j in range(len(soma_region)):
            if get_single_neuron_projection_ratio.remove_layer_information(soma_region[j]) == i:
                c.append(j)
        neuron_idx.append(c)
    return neuron_idx


def get_neuron_number_in_different_layer(soma_region, Layer):
    neuron_idx = []
    for i in Layer:
        c = []
        for j in range(len(soma_region)):
            if get_single_neuron_projection_ratio.determine_layer(soma_region[j]) == i:
                c.append(j)
        neuron_idx.append(c)
    return neuron_idx


def draw_error_bar(x1, x2, x3, x4, x5):
    from scipy import stats
    x1_mean = np.mean(x1)
    x1_sem = stats.sem(x1)
    x2_mean = np.mean(x2)
    x2_sem = stats.sem(x2)
    x3_mean = np.mean(x3)
    x3_sem = stats.sem(x3)
    x4_mean = np.mean(x4)
    x4_sem = stats.sem(x4)
    x5_mean = np.mean(x5)
    x5_sem = stats.sem(x5)
    x_mean = [x1_mean, x2_mean, x3_mean, x4_mean, x5_mean]
    x_sem = [x1_sem, x2_sem, x3_sem, x4_sem, x5_sem]
    return x_mean, x_sem


def draw_error_bar1(x1, x2, x3, x4):
    from scipy import stats
    x1_mean = np.mean(x1)
    x1_sem = stats.sem(x1)
    x2_mean = np.mean(x2)
    x2_sem = stats.sem(x2)
    x3_mean = np.mean(x3)
    x3_sem = stats.sem(x3)
    x4_mean = np.mean(x4)
    x4_sem = stats.sem(x4)

    x_mean = [x1_mean, x2_mean, x3_mean, x4_mean]
    x_sem = [x1_sem, x2_sem, x3_sem, x4_sem]
    return x_mean, x_sem

def draw_error_bar2(x):
    from scipy import stats
    x_mean = []
    x_sem = []
    for i in x:
        i_mean = np.mean(i)
        i_sem = stats.sem(i)
        x_mean.extend([i_mean])
        x_sem.extend([i_sem])
    return x_mean, x_sem

if __name__ == "__main__":
    data_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data = pd.read_csv(data_path, usecols=[0, 1, 2, 3, 4])
    data = data.values
    neuron_id = data[:, 1]
    soma_region = data[:, 2]
    ipsi_termials = data[:, 3]
    contra_terminals = data[:, 4]


    data1_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_neuron_projection_full_summary_cortex_ceil_av.csv'
    data1 = pd.read_csv(data1_path)
    data1 = data1.values
    neuron_full_id = data1[:, 1]
    neuron_full_ipsi_maxdepth = data1[:, 2]
    neuron_full_contra_maxdepth = data1[:, 3]
    neuron_full_ipsi_tlen = data1[:, 4]
    neuron_full_ipsi_maxtlen = data1[:, 5]
    neuron_full_contra_tlen = data1[:, 6]
    neuron_full_contra_maxtlen = data1[:, 7]
    # neuron_full_bn = data1[:, 8]  # branch number
    # neuron_full_bd = data1[:, 9]  # branch degree kind
    # neuron_full_bmd = data1[:, 10]  # branch max degree
    idx = [j for i in range(len(neuron_id)) for j in range(len(neuron_full_id)) if neuron_id[i] == neuron_full_id[j]]
    neuron_ipsi_maxdepth = neuron_full_ipsi_maxdepth[idx]
    neuron_contra_maxdepth = neuron_full_contra_maxdepth[idx]
    neuron_ipsi_tlen = neuron_full_ipsi_tlen[idx]
    neuron_contra_tlen = neuron_full_contra_tlen[idx]
    neuron_ipsi_max_tlen = neuron_full_ipsi_maxtlen[idx]
    neuron_contra_max_tlen = neuron_full_contra_maxtlen[idx]
    # neuron_bn = neuron_full_bn[idx]
    # neuron_bd = neuron_full_bd[idx]
    # neuron_bmd = neuron_full_bmd[idx]
    neuron_ipsi_tlen_list = [get_number_list(i) for i in list(neuron_ipsi_tlen)]  # str to num
    neuron_contra_tlen_list = [get_number_list(i) for i in list(neuron_contra_tlen)]  # str to num

    # branch number of terminals of projected to cortex region  # 3732 neurons
    data2_path = '/data/single_neuro_tracing/code/preprocess/IT_neuron_projection_to_cortex_branch_number_repair.csv'
    data2 = pd.read_csv(data2_path)
    data2 = data2.values
    neuron_id2 = data2[:, 1]
    neuron_cortex_bn = data2[:, 2]
    idx_2 = [j for i in range(len(neuron_id)) for j in range(len(neuron_id2)) if neuron_id[i] == neuron_id2[j]]
    neuron_bn = neuron_cortex_bn[idx_2]

    [c1, c11, c111, c112, c12, c121, c122, c13, c14, c141, c142, c2] = find_5_type_neuron_idx(ipsi_termials,
                                                                                              contra_terminals,
                                                                                              soma_region)
    save_path = '/data/single_neuro_tracing/code/preprocess/new_results/'
    subtype = ['c11', 'c12', 'c13', 'c14', 'c2']
    # ipsi_max_terminals_c11 = neuron_ipsi_max_tlen[c11]
    # contra_max_terminals_c11 = neuron_contra_max_tlen[c11]
    # ipsi_max_terminals_c12 = neuron_ipsi_max_tlen[c12]
    # contra_max_terminals_c12 = neuron_contra_max_tlen[c12]
    # ipsi_max_terminals_c13 = neuron_ipsi_max_tlen[c13]
    # contra_max_terminals_c13 = neuron_contra_max_tlen[c13]
    # ipsi_max_terminals_c14 = neuron_ipsi_max_tlen[c14]
    # contra_max_terminals_c14 = neuron_contra_max_tlen[c14]
    # ipsi_max_terminals_c2 = neuron_ipsi_max_tlen[c2]
    # contra_max_terminals_c2 = neuron_contra_max_tlen[c2]
    # ipsi_ml_mean, ipsi_ml_sem = draw_error_bar(ipsi_max_terminals_c11, ipsi_max_terminals_c12, ipsi_max_terminals_c13,
    #                                            ipsi_max_terminals_c14, ipsi_max_terminals_c2)
    # contra_ml_mean, contra_ml_sem = draw_error_bar(contra_max_terminals_c11, contra_max_terminals_c12,
    #                                                contra_max_terminals_c13, contra_max_terminals_c14,
    #                                                contra_max_terminals_c2)
    # plt.figure()
    # x = np.arange(0, 5, 1)
    # plt.errorbar(x, ipsi_ml_mean, yerr=ipsi_ml_sem, fmt='o', capsize=3, markersize=3)
    # plt.errorbar(x, contra_ml_mean, yerr=contra_ml_sem, fmt='o', capsize=3, markersize=3)
    # plt.xticks(x, subtype)
    # plt.legend(['ipsilateral', 'contralateral'])
    # save_fn = save_path + 'terminal_length_in_subtype_mean_sem' + '.eps'
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.rcParams['ps.fonttype'] = 42
    # plt.savefig(save_fn, format='pdf')
    # plt.close()
    # h_statistic, p_value = kruskal(ipsi_max_terminals_c11, ipsi_max_terminals_c12, ipsi_max_terminals_c13,
    #                                ipsi_max_terminals_c14, ipsi_max_terminals_c2)
    # y = [list(ipsi_max_terminals_c11), list(ipsi_max_terminals_c12), list(ipsi_max_terminals_c13),
    #      list(ipsi_max_terminals_c14), list(ipsi_max_terminals_c2)]
    # y_c = [list(contra_max_terminals_c11), list(contra_max_terminals_c12), list(contra_max_terminals_c13),
    #        list(contra_max_terminals_c14), list(contra_max_terminals_c2)]
    # ipsi_terminal_max_length = []
    # types = []
    # projection = []
    # m = 0
    # for i in y:
    #     leng = len(i)
    #     ipsi_terminal_max_length.extend(i)
    #
    #     types.extend([subtype[m] for k in range(len(i))])
    #     m = m + 1
    # projection.extend(['ipsilateral'] * len(ipsi_terminal_max_length))
    # m = 0
    # contra_terminal_max_length = []
    # for j in y_c:
    #     leng = len(j)
    #     contra_terminal_max_length.extend(j)
    #
    #     types.extend([subtype[m] for k in range(len(j))])
    #     m = m + 1
    # projection.extend(['contralateral'] * len(contra_terminal_max_length))
    # terminal_max_length = []
    # terminal_max_length.extend(ipsi_terminal_max_length)
    # terminal_max_length.extend(contra_terminal_max_length)
    # a = {'terminal_length': terminal_max_length, 'types': types, 'projection':projection}
    # df = pd.DataFrame(a)
    # seaborn.set_style('white')
    # # draw violin plot
    # plt.figure()
    # seaborn.violinplot(y='terminal_length', x='types', hue='projection', data=df)
    # plt.rcParams['ps.fonttype'] = 42
    # plt.rcParams['pdf.fonttype'] = 42
    # save_fn = save_path + 'terminal_max_length_in_subtype_violin' + '.eps'
    # plt.savefig(save_fn, format='pdf')
    # plt.close()
    # tukey_result = pairwise_tukeyhsd(a['terminal_length'], a['types'])
    # print(tukey_result.summary())

    neuron_bn_c11 = neuron_bn[c11]
    neuron_bn_c12 = neuron_bn[c12]
    neuron_bn_c13 = neuron_bn[c13]
    neuron_bn_c14 = neuron_bn[c14]
    neuron_bn_c2 = neuron_bn[c2]
    bn_mean, bn_sem = draw_error_bar(neuron_bn_c11, neuron_bn_c12, neuron_bn_c13,
                                     neuron_bn_c14, neuron_bn_c2)
    plt.figure()
    x = np.arange(0, 5, 1)
    plt.bar(x, bn_mean, color='#E26A6A')
    for i in x:
        plt.plot([i, i], [bn_mean[i] - bn_sem[i], bn_mean[i] + bn_sem[i]], color=[0, 0, 0])
    plt.xticks(x, subtype)
    save_fn = save_path + 'branch_number_in_subtype_mean_sem_bar' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()
    # # plt.errorbar(x, bn_mean, yerr=bn_sem, fmt='o', capsize=3, markersize=3)
    # # plt.xticks(x, subtype)
    # # save_fn = save_path + 'branch_number_in_subtype_mean_sem' + '.eps'
    # # plt.rcParams['pdf.fonttype'] = 42
    # # plt.rcParams['ps.fonttype'] = 42
    # # plt.savefig(save_fn, format='pdf')
    # # plt.close()
    # # h_statistic, p_value = kruskal(neuron_bn_c11, neuron_bn_c12, neuron_bn_c13,
    # #                                neuron_bn_c14, neuron_bn_c2)
    # branch_number = []
    # types = []
    # m = 0
    # for i in [list(neuron_bn_c11), list(neuron_bn_c12), list(neuron_bn_c13),
    #           list(neuron_bn_c14), list(neuron_bn_c2)]:
    #     bl = len(i)
    #     branch_number.extend(i)
    #
    #     types.extend([subtype[m] for k in range(len(i))])
    #     m = m + 1
    #
    # a = {'branch_number': branch_number, 'types': types}
    # df = pd.DataFrame(a)
    # df = df[df['branch_number'] < 500]
    # # df = df[df['branch_number'] != 2095]
    # # draw violin plot
    # plt.figure()
    # seaborn.violinplot(y='branch_number', x='types', data=df)
    # plt.rcParams['ps.fonttype'] = 42
    # plt.rcParams['pdf.fonttype'] = 42
    # save_fn = save_path + 'branch_number_in_subtype_violin_deo' + '.eps'
    # plt.savefig(save_fn, format='pdf')
    # plt.close()
    # tukey_result = pairwise_tukeyhsd(a['branch_number'], a['types'])
    # print(tukey_result.summary())


   