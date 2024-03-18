import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr, pearsonr
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd


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


def draw_scatter(d1, d2):
    plt.figure()
    plt.scatter(d1, d2)
    slope, intercept, r_value, p_value, std_err = linregress(list(d1), list(d2))
    # print('slope: %f' % slope)
    # print('r value: %f' % r_value)
    print('R2: %f' % np.square(r_value))
    print('P value: %f' % p_value)
    y_fit = slope * d1 + intercept
    R2 = 1 - np.sum(np.square(d2 - y_fit)) / np.sum(np.square(d2 - np.mean(d2)))
    # print('R2: %f' % R2)
    plt.plot(d1, y_fit)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig('./new_results/figures/axon_and_dendrite_new.eps', format='pdf')
    plt.close()



def draw_scatter_for_5_type(d1, d2):
    types = ['IB', 'BC', 'B', 'IC', 'IBC']
    for i in range(5):
        plt.figure()
        dendri = d1[i]
        axon = d2[i]
        slope, intercept, r_value, p_value, std_err = linregress(list(dendri), list(axon))
        y_fit = slope * dendri + intercept


        # print('slope: %f' % slope)
        # print('r value: %f' % r_value)
        print('R2: %f' % np.square(r_value))
        print('P value: %f' % p_value)
        R2 = 1 - np.sum(np.square(axon - y_fit)) / np.sum(np.square(axon - np.mean(axon)))
        # print('R2: %f' % R2)
        plt.scatter(dendri, axon)
        plt.plot(dendri, y_fit)
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        plt.savefig('./new_results/figures/' + types[i] + '_axon_and_dendrite_new.eps', format='pdf')
        plt.close()


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


if __name__ == "__main__":
    ## data_path: cortex neuron projected bilaterally to cortex regions: 3004 neurons
    data_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data = pd.read_csv(data_path, usecols=[0, 1, 2, 3, 4])
    data = data.values
    neuron_id = data[:, 1]
    ipsi_cortex_terminals = data[:, 3]
    contra_cortex_terminals = data[:, 4]
    ipsi_cortex_terminals = remove_some_information(ipsi_cortex_terminals)
    contra_cortex_terminals = remove_some_information(contra_cortex_terminals)

    # # axon length data   # 3732 neurons
    data1_path = '/data/single_neuro_tracing/code/preprocess/IT_neuron_projection_to_cortex_axon_length_repair.csv'
    data1 = pd.read_csv(data1_path)
    data1 = data1.values
    neuron_id1 = data1[:, 1]
    axon_length = data1[:, 2]
    idx = [j for i in range(len(neuron_id)) for j in range(len(neuron_id1)) if neuron_id[i] == neuron_id1[j]]
    neuron_id1_ = neuron_id1[idx]  # bilateral projection
    axon_length = axon_length[idx]

    # # dendrites length data   # 1027 neurons
    data2_path = '/data/single_neuro_tracing/code/preprocess/IT_neuron_projection_dendrites_length_for_new_data.csv'
    data2 = pd.read_csv(data2_path)
    data2 = data2.values
    neuron_id2 = data2[:, 1]  # neuron having dendrites information from downloading data
    dendrites_length = data2[:, 2]

    idx = [j for i in range(len(neuron_id)) for j in range(len(neuron_id2)) if neuron_id[i] == neuron_id2[j]]
    neuron_id2_ = neuron_id2[idx]  # bilateral projection
    dendrites_length = dendrites_length[idx]

    idx_1 = [i for i in range(len(neuron_id1_)) if neuron_id1_[i] in neuron_id2_]
    ipsi_cortex_terminals = np.array(ipsi_cortex_terminals)[idx_1]
    contra_cortex_terminals = np.array(contra_cortex_terminals)[idx_1]
    axon_length = axon_length[idx_1]

    data4 = pd.read_csv(
        '/data/single_neuro_tracing/code/preprocess/IT_cortex_bilateral_neuron_dendrites_type_from_list.csv')
    data4 = data4.values
    neuron_id4 = data4[:, 1]  # neuron having dendrites information from their excel
    types = data4[:, 2]

    # # check difference
    diff1 = list(set(list(neuron_id2_)).difference(set(neuron_id4)))
    diff2 = list(set(list(neuron_id4)).difference(set(neuron_id2_)))

    idx_diff1 = [i for i in range(len(neuron_id2_)) if neuron_id2_[i] not in diff1]
    idx_diff2 = [i for i in range(len(neuron_id4)) if neuron_id4[i] not in diff2]

    ipsi_cortex_terminals = np.array(ipsi_cortex_terminals)[idx_diff1]
    contra_cortex_terminals = np.array(contra_cortex_terminals)[idx_diff1]
    axon_length = axon_length[idx_diff1]
    dendrites_length = dendrites_length[idx_diff1]
    # draw_scatter(dendrites_length, axon_length)
    r_value = []
    p_value = []
    r1, p1 = pearsonr(dendrites_length, axon_length)
    r_value.append(r1)
    p_value.append(p1)

    [c1, c11, c111, c112, c12, c121, c122, c13, c14, c141, c142, c2] = find_5_type_neuron_idx(ipsi_cortex_terminals,
                                                                                              contra_cortex_terminals)

    dendrites_length_subtype = []
    axon_length_subtype = []
    for i in [c11, c12, c13, c14, c2]:
        dendrites_length_subtype.append(dendrites_length[i])
        axon_length_subtype.append(axon_length[i])

    for i in range(5):
        d1 = dendrites_length_subtype[i]
        a1 = axon_length_subtype[i]
        # a = {'dendrites_length': d1, 'axon_length': a1}
        # df = pd.DataFrame(a)
        # df.to_csv('./' + str(i+1) + '.csv')
        r1, p1 = pearsonr(d1, a1)
        r_value.append(r1)
        p_value.append(p1)
    print(r_value)
    print(p_value)

    dendrites_length_mean, dendrites_length_sem = draw_error_bar2(dendrites_length_subtype)
    import seaborn
    seaborn.set_style('white')
    plt.figure()
    x = np.arange(0, 5, 1)
    plt.bar(x, dendrites_length_mean, color='#E26A6A')
    for i in x:
        plt.plot([i, i], [dendrites_length_mean[i] - dendrites_length_sem[i], dendrites_length_mean[i] + dendrites_length_sem[i]], color=[0, 0, 0])
    subtype = ['IB', 'BC', 'B', 'IC', 'IBC']
    plt.xticks(x, subtype)
    plt.ylim(5, 8)
    save_fn = './new_results/figures/total_dendrites_length_in_subtype_mean_sem_bar' + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()
    # subtype = ['IB', 'BC', 'B', 'IC', 'IBC']
    # import seaborn
    # seaborn.set_style('darkgrid')
    # plt.figure()
    # x = np.arange(0, 5, 1)
    # plt.errorbar(x, dendrites_length_mean, yerr=dendrites_length_sem, fmt='o', capsize=3, markersize=3)
    # plt.xticks(x, subtype)
    # save_fn = './new_results/figures/total_dendrites_length_in_subtype_mean_sem_new' + '.eps'
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.rcParams['ps.fonttype'] = 42
    # plt.savefig(save_fn, format='pdf')
    # plt.close()

    length = []
    types = []
    m = 0
    for i in range(len(dendrites_length_subtype)):
        length.extend(list(dendrites_length_subtype[i]))
        types.extend([subtype[m] for i in range(len(dendrites_length_subtype[i]))])
        m = m + 1
    a = {'length':length, 'type': types}
    df = pd.DataFrame(a)
    plt.figure()
    seaborn.violinplot(x='type', y='length', data=df)
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    save_fn = './new_results/figures/total_dendrites_length_in_subtype_violin' + '.eps'
    plt.savefig(save_fn, format='pdf')
    plt.close()
    tukey_result = pairwise_tukeyhsd(a['length'], a['type'])
    print(tukey_result.summary())

    draw_scatter_for_5_type(dendrites_length_subtype, axon_length_subtype)
