import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl

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
    c113 = []
    c12 = []  # ipsi == 0, contra != 0, both != 0
    c121 = []
    c122 = []
    c123 = []
    c13 = []  # ipsi == 0, contra == 0, both != 0
    c14 = []  # ipsi != 0, contra != 0, both == 0
    c141 = []  # ipsi >= contra
    c142 = []  # ipsi < contra
    c143 = []
    c2 = []  # only ipsi, only contra, both all != 0
    c21 = []  # ipsi  is the best
    c211 = []  # i = B = C
    c212 = []  # i = B
    c213 = []  # I = C
    c214 = []  # i is the best
    c22 = []  # B is the best
    c23 = []  # C is the best
    c231 = []  # C = B
    for j in range(len(ipsi_termials)):
        r1 = s1[j]  # get the ratio of single neuron in three types
        r2 = s2[j]
        r3 = s3[j]
        if r1 != 0 and r2 != 0 and r3 != 0:
            c2.append(j)
            if np.argmax(np.array([r1, r2, r3])) == 0:
                if len(np.where(np.array([r1, r2, r3] == np.max(np.array([r1, r2, r3]))))[0]) == 3:
                    c211.append(j)   # I = B = C
                elif len(np.where(np.array([r1, r2, r3] == np.max(np.array([r1, r2, r3]))))[0]) == 2:
                    if np.where(np.array([r1, r2, r3] == np.max(np.array([r1, r2, r3]))))[0][1] == 2:
                        c212.append(j)  # I = B are the best
                    elif np.where(np.array([r1, r2, r3] == np.max(np.array([r1, r2, r3]))))[0][1] == 1:
                        c213.append(j)  # I = c
                else:
                    c214.append(j)
                c21.append(j)
            elif np.argmax(np.array([r1, r2, r3])) == 2:
                c22.append(j)
            elif np.argmax(np.array([r1, r2, r3])) == 1:
                if len(np.where(np.array([r1, r2, r3] == np.max(np.array([r1, r2, r3]))))[0]) > 1:
                    c231.append(j)
                c23.append(j)
        else:
            c1.append(j)
            if r1 != 0 and r2 == 0 and r3 != 0:
                c11.append(j)
                if r1 > r3:
                    c111.append(j)   # i>B
                elif r1 < r3:
                    c112.append(j)  # I < B
                else:
                    c113.append(j)   # I = B
            elif r1 == 0 and r2 != 0 and r3 != 0:
                c12.append(j)
                if r2 > r3:
                    c121.append(j)  # c >B
                elif r2 < r3:
                    c122.append(j)  # C < B
                else:
                    c123.append(j)  # C = B
            elif r1 == 0 and r2 == 0 and r3 != 0:
                c13.append(j)
            elif r1 != 0 and r2 != 0 and r3 == 0:
                c14.append(j)
                if r1 > r2:        # I > C
                    c141.append(j)
                elif r1 < r2:     # I < C
                    c142.append(j)
                else:
                    c143.append(j)  # I = C
    return c1, c11, c111, c112, c113, c12, c121, c122, c123, c13, c14, c141, c142, c143, c2

# def get_ratio_account_for_dendrites_subtype(dendrite):
#     number = len(np.where(dendrite)[0])
#     types = list(set(dendrite))
#     types_sorted = sorted(types)
#     ratio = {}
#     count = {}
#     for i in types_sorted:
#         if i == 0:
#             continue
#         else:
#             key = str(i)
#             number_key = len(np.where(dendrite == i)[0])
#             ratio[key] = number_key / number
#             count[key] = number_key
#     return ratio, count

def get_ratio_account_for_dendrites_subtype(dendrite):
    number = len(np.where(dendrite)[0])
    types = [i for i in range(1, 39)]
    ratio = {}
    count = {}
    for i in types:
        if i == 0:
            continue
        else:
            key = str(i)
            number_key = len(np.where(dendrite == i)[0])
            ratio[key] = number_key / number
            count[key] = number_key
    return ratio, count

def get_ratio_account_for_dendrites_subtype_new(dendrite, dendrite_all):
    # number = len(np.where(dendrite)[0])
    types = [i for i in range(1, 39)]
    ratio = {}
    count = {}
    for i in types:
        if i == 0:
            continue
        else:
            key = str(i)
            number_key = len(np.where(dendrite == i)[0])
            number = len(np.where(dendrite_all == i)[0])
            if number == 0:
                number = 0.000001
            ratio[key] = number_key / number
            count[key] = number_key
    return ratio, count

def draw_stem_figure(dens):
    x = []
    y = []
    for key, value in dens.items():
        x.append(key)
        y.append(value)
    plt.stem(x, y, markerfmt='ro', linefmt='r-')
    plt.xticks(rotation=90)
    save_path = '/data/single_neuro_tracing/code/preprocess/new_results/'
    name = 'dendrites_bilateral_count_stem'
    save_fn = save_path + name + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()

def draw_dot_figure_for_different_subtype(d1, d2, d3, d4, d5):
    matrix = []
    for i in [d1, d2, d3, d4, d5]:
        vec = []
        for k, v in i.items():
            vec.extend([v])
        matrix.append(vec)

    matrix = np.array(matrix)
    matrix_flip = np.flip(matrix, axis=0)
    rows, cols = matrix.shape
    cmap1 = mpl.cm.get_cmap("Set1")
    #colors = cmap1(np.linspace(0, 1, cmap1.N))[0:5]
    colors = ['#9467BD', '#D62728', '#2CA02C', '#FF7F0E', '#1F77B4']
    i = 0
    fig, ax = plt.subplots()
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.grid(True)
    ax.set_axisbelow(True)
    for j in range(rows):
        for k in range(cols):
            size = matrix_flip[j, k] * 2000
            # if 0 < r <= 0.2:
            #     size = 60
            # elif 0.2 < r <= 0.4:
            #     size = 100
            # elif 0.4 < r <= 0.6:
            #     size = 140
            # elif 0.6 < r <= 0.8:
            #     size = 180
            # elif 0.8 < r <= 1:
            #     size = 220
            # elif r == 0:
            #     size = 0
            ax.scatter(k, j, s=size, c=colors[i])
        i = i + 1

    plt.yticks(np.arange(0, 5), ['IBC', 'IC', 'B', 'BC', 'IB'])
    plt.xticks(np.arange(0, 38), np.arange(1, 39), rotation=90)

    save_path = '/data/single_neuro_tracing/code/preprocess/new_results/'
    name = 'dendrites_bilateral_ratio_for_different_dubtype_dot'
    save_fn = save_path + name + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()

if __name__ == "__main__":
    data1_path = '/data/single_neuro_tracing/code/preprocess/new_results/dendrites_classify.csv'
    data1 = pd.read_csv(data1_path)
    data1 = data1.values
    neuron_id1 = data1[:, 1]
    soma_region1 = data1[:, 2]
    dendrites = data1[:, 3]

    # 3704 neurons
    data2_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projection_ceil.csv'
    data2 = pd.read_csv(data2_path)
    data2 = data2.values
    neuron_id2 = data2[:, 1]
    soma_region2 = data2[:, 2]
    ipsi_cortex_terminals_2 = data2[:, 3]
    contra_cortex_terminals_2 = data2[:, 4]
    ipsi_cortex_terminals_2 = remove_some_information(ipsi_cortex_terminals_2)
    contra_cortex_terminals_2 = remove_some_information(contra_cortex_terminals_2)
    idx_2 = [j for i in range(len(neuron_id2)) for j in range(len(neuron_id1)) if neuron_id2[i] == neuron_id1[j]]
    dendrites_2 = dendrites[idx_2]
    idx_oi = [i for i in range(len(ipsi_cortex_terminals_2)) if len(ipsi_cortex_terminals_2[i]) != 0
              and len(contra_cortex_terminals_2[i]) == 0]
    idx_oc = [i for i in range(len(ipsi_cortex_terminals_2)) if len(ipsi_cortex_terminals_2[i]) == 0
              and len(contra_cortex_terminals_2[i]) != 0]
    dendrites_oi = dendrites_2[idx_oi]
    dendrites_oc = dendrites_2[idx_oc]
    dendrites_oi_ratio, dendrites_oi_count = get_ratio_account_for_dendrites_subtype(dendrites_oi)
    dendrites_oc_ratio, dendrites_oc_count = get_ratio_account_for_dendrites_subtype(dendrites_oc)

    #3004 neurons
    data3_path = '/mnt/Data16Tb/Data/feiyao/Allen_mouse/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data3 = pd.read_csv(data3_path)
    data3 = data3.values
    neuron_id3 = data3[:, 1]
    soma_region3 = data3[:, 2]
    ipsi_cortex_terminals_3 = data3[:, 3]
    contra_cortex_terminals_3 = data3[:, 4]
    ipsi_cortex_terminals_3 = remove_some_information(ipsi_cortex_terminals_3)
    contra_cortex_terminals_3 = remove_some_information(contra_cortex_terminals_3)
    idx_3 = [j for i in range(len(neuron_id3)) for j in range(len(neuron_id1)) if neuron_id3[i] == neuron_id1[j]]
    dendrites_3 = dendrites[idx_3]
    dendrites_3_ratio, dendrites_3_count = get_ratio_account_for_dendrites_subtype(dendrites_3)
    draw_stem_figure(dendrites_3_count)

    [c1, c11, c111, c112, c113, c12, c121, c122, c123, c13, c14, c141, c142, c143, c2] = find_5_type_neuron_idx(
        ipsi_cortex_terminals_3, contra_cortex_terminals_3)

    dendrites_c11 = dendrites_3[c11]
    dendrites_c12 = dendrites_3[c12]
    dendrites_c13 = dendrites_3[c13]
    dendrites_c14 = dendrites_3[c14]
    dendrites_c2 = dendrites_3[c2]

    dendrites_c11_ratio, dendrites_c11_count = get_ratio_account_for_dendrites_subtype(dendrites_c11)
    dendrites_c12_ratio, dendrites_c12_count = get_ratio_account_for_dendrites_subtype(dendrites_c12)
    dendrites_c13_ratio, dendrites_c13_count = get_ratio_account_for_dendrites_subtype(dendrites_c13)
    dendrites_c14_ratio, dendrites_c14_count = get_ratio_account_for_dendrites_subtype(dendrites_c14)
    dendrites_c2_ratio, dendrites_c2_count = get_ratio_account_for_dendrites_subtype(dendrites_c2)
    #draw_dot_figure_for_different_subtype(dendrites_c11_count, dendrites_c12_count, dendrites_c13_count,
                                             # dendrites_c14_count, dendrites_c2_count)
    draw_dot_figure_for_different_subtype(dendrites_c11_ratio, dendrites_c12_ratio, dendrites_c13_ratio,
                                          dendrites_c14_ratio, dendrites_c2_ratio)






    print(1)