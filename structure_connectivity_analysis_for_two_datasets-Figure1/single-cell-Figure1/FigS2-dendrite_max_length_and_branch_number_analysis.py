import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr
from matplotlib import pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import re
import seaborn


def determine_layer(SOMA):
    ret = re.search('\d+', SOMA)
    if ret is None:
        ret = 'none'
    else:
        ret = ret.group()
    return ret


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


def get_neuron_number_in_different_layer(soma_region, Layer):
    neuron_idx = []
    for i in Layer:
        c = []
        for j in range(len(soma_region)):
            if determine_layer(soma_region[j]) == i:
                c.append(j)
        neuron_idx.append(c)
    return neuron_idx


def get_neuron_number_in_different_region(soma_region, region_sets):
    neuron_idx = []
    for i in region_sets:
        c = []
        for j in range(len(soma_region)):
            if remove_layer_information(soma_region[j]) == i:
                c.append(j)
        neuron_idx.append(c)
    return neuron_idx


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
    ## data_path: cortex neuron projected bilaterally to cortex regions: 3004 neurons
    data_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data = pd.read_csv(data_path, usecols=[0, 1, 2, 3, 4])
    data = data.values
    neuron_id = data[:, 1]
    soma_region = data[:, 2]
    ipsi_cortex_terminals = data[:, 3]
    contra_cortex_terminals = data[:, 4]
    ipsi_cortex_terminals = remove_some_information(ipsi_cortex_terminals)
    contra_cortex_terminals = remove_some_information(contra_cortex_terminals)

    # # dendrite max length and branch number data   # 1919 neurons
    data1_path = '/data/single_neuro_tracing/code/preprocess/IT_neuron_projection_dendrites_branch_and_max_length.csv'
    data1 = pd.read_csv(data1_path)
    data1 = data1.values
    neuron_id1 = data1[:, 1]
    dendrite_length = data1[:, 2]
    dendrite_branch_number = data1[:, 3]
    idx = [j for i in range(len(neuron_id)) for j in range(len(neuron_id1)) if neuron_id[i] == neuron_id1[j]]
    neuron_id1_ = neuron_id1[idx]  # bilateral projection
    dendrite_length = dendrite_length[idx]
    dendrite_branch_number = dendrite_branch_number[idx]

    idx_1 = [i for i in range(len(neuron_id)) if neuron_id[i] in neuron_id1_]
    ipsi_cortex_terminals = np.array(ipsi_cortex_terminals)[idx_1]
    contra_cortex_terminals = np.array(contra_cortex_terminals)[idx_1]

    [c1, c11, c111, c112, c12, c121, c122, c13, c14, c141, c142, c2] = find_5_type_neuron_idx(ipsi_cortex_terminals,
                                                                                              contra_cortex_terminals)
    save_path = '/mnt/Data16Tb/Data/feiyao/Allen_mouse/data/single_neuro_tracing/code/preprocess/new_results/'
    subtype = ['IB', 'BC', 'B', 'IC', 'IBC']
    dendrite_length_subtype = []
    for i in [c11, c12, c13, c14, c2]:
        dendrite_length_subtype.append(dendrite_length[i])

    length_mean, length_sem = draw_error_bar2(dendrite_length_subtype)
    seaborn.set_style('white')
    plt.figure()
    x = np.arange(0, 5, 1)
    plt.bar(x, length_mean, color='#E26A6A')
    for i in x:
        plt.plot([i, i], [length_mean[i] - length_sem[i], length_mean[i] + length_sem[i]], color=[0, 0, 0])
    plt.xticks(x, subtype)
    plt.ylim(0.3, 0.6)
    save_fn = save_path + 'dendrite_max_length_in_subtype_mean_sem_bar' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()

    length = []
    types = []
    m = 0
    for i in [c11, c12, c13, c14, c2]:
        length.extend(dendrite_length[i])
        types.extend([subtype[m] for k in range(len(i))])
        m = m + 1

    a = {'dendrite_length': length, 'types': types}
    df = pd.DataFrame(a)
    plt.figure()
    seaborn.violinplot(x='types', y='dendrite_length', data=df)
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    save_fn = save_path + 'dendrite_max_length_in_subtype_violin' + '.eps'
    plt.savefig(save_fn, format='pdf')
    plt.close()
    tukey_result = pairwise_tukeyhsd(a['dendrite_length'], a['types'])
    print(tukey_result.summary())

    # branch number
    dendrite_bn_subtype = []
    for i in [c11, c12, c13, c14, c2]:
        dendrite_bn_subtype.append(dendrite_branch_number[i])

    bn_mean, bn_sem = draw_error_bar2(dendrite_bn_subtype)
    seaborn.set_style('white')
    plt.figure()
    x = np.arange(0, 5, 1)
    plt.bar(x, bn_mean, color='#E26A6A')
    for i in x:
        plt.plot([i, i], [bn_mean[i] - bn_sem[i], bn_mean[i] + bn_sem[i]], color=[0, 0, 0])
    plt.xticks(x, subtype)
    plt.ylim(30, 50)
    save_fn = save_path + 'dendrite_branch_number_in_subtype_mean_sem_bar' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()

    bn = []
    types = []
    m = 0
    for i in [c11, c12, c13, c14, c2]:
        bn.extend(dendrite_branch_number[i])
        types.extend([subtype[m] for k in range(len(i))])
        m = m + 1

    a = {'dendrite_branch_number': bn, 'types': types}
    df = pd.DataFrame(a)
    plt.figure()
    seaborn.violinplot(x='types', y='dendrite_branch_number', data=df)
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    save_fn = save_path + 'dendrite_branch_number_in_subtype_violin' + '.eps'
    plt.savefig(save_fn, format='pdf')
    plt.close()
    tukey_result = pairwise_tukeyhsd(a['dendrite_branch_number'], a['types'])
    print(tukey_result.summary())
