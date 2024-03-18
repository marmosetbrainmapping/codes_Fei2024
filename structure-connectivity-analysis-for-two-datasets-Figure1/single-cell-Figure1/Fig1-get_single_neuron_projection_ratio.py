import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt
import seaborn
from scipy import signal
from plot_colormap import define_colormap4


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


def determine_layer(SOMA):
    ret = re.search('\d+', SOMA)
    if ret is None:
        ret = 'none'
    else:
        ret = ret.group()
    return ret


def get_number_distribution_of_neuron_in_layer_of_region(region_sets, soma_region):
    X = np.array([1, 2, 3, 4])
    L1 = []
    L2 = []
    L5 = []
    L6 = []
    count_sum = []
    for i in region_sets:
        count = 0
        l1 = 0
        l2 = 0
        l5 = 0
        l6 = 0
        for j in soma_region:
            if remove_layer_information(j) == i:
                count = count + 1
                if determine_layer(j) == '1':
                    l1 = l1 + 1
                elif determine_layer(j) == '2':
                    l2 = l2 + 1
                elif determine_layer(j) == '5':
                    l5 = l5 + 1
                elif determine_layer(j) == '6':
                    l6 = l6 + 1
        if count == 0:
            L1.append(0)
            L2.append(0)
            L5.append(0)
            L6.append(0)
        else:
            L1.append(l1)  # number of neuron in layer 1 of this region
            L2.append(l2)
            L5.append(l5)
            L6.append(l6)
        count_sum.append(count)
        d1 = np.ones(l1) * 1
        d2 = np.ones(l2) * 2
        d5 = np.ones(l5) * 3
        d6 = np.ones(l6) * 4
        data = np.concatenate([d1, d2, d5, d6], axis=0)

        # seaborn.kdeplot(data)
        # plt.show()
        print(count)
    distribution = [L1, L2, L5, L6]
    distribution = np.array(distribution)
    return distribution


def get_number_of_different_type_neuron_in_different_region(region_sets, soma_region):
    A1 = []
    A2 = []
    A3 = []
    A4 = []
    A5 = []
    A6 = []
    A7 = []
    A8 = []
    A9 = []
    A10 = []
    A11 = []

    for j in [c11, c12, c13, c14, c2]:
        regions = soma_region[j]
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        a5 = 0
        a6 = 0
        a7 = 0
        a8 = 0
        a9 = 0
        a10 = 0
        a11 = 0
        for k in regions:
            if remove_layer_information(k) == region_sets[0]:
                a1 = a1 + 1
            elif remove_layer_information(k) == region_sets[1]:
                a2 = a2 + 1
            elif remove_layer_information(k) == region_sets[2]:
                a3 = a3 + 1
            elif remove_layer_information(k) == region_sets[3]:
                a4 = a4 + 1
            elif remove_layer_information(k) == region_sets[4]:
                a5 = a5 + 1
            elif remove_layer_information(k) == region_sets[5]:
                a6 = a6 + 1
            elif remove_layer_information(k) == region_sets[6]:
                a7 = a7 + 1
            elif remove_layer_information(k) == region_sets[7]:
                a8 = a8 + 1
            elif remove_layer_information(k) == region_sets[8]:
                a9 = a9 + 1
            elif remove_layer_information(k) == region_sets[9]:
                a10 = a10 + 1
            elif remove_layer_information(k) == region_sets[10]:
                a11 = a11 + 1
        A1.append(a1)
        A2.append(a2)
        A3.append(a3)
        A4.append(a4)
        A5.append(a5)
        A6.append(a6)
        A7.append(a7)
        A8.append(a8)
        A9.append(a9)
        A10.append(a10)
        A11.append(a11)
    return A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11


def get_distribution_of_neuron_in_layer_of_region(region_sets, soma_region):
    X = np.array([1, 2, 3, 4])
    L1 = []
    L2 = []
    L5 = []
    L6 = []
    count_sum = []
    # plt.figure()
    for i in region_sets:  # each PFC region
        count = 0
        l1 = 0
        l2 = 0
        l5 = 0
        l6 = 0
        for j in soma_region:
            if remove_layer_information(j) == i:
                count = count + 1
                if determine_layer(j) == '1':
                    l1 = l1 + 1
                elif determine_layer(j) == '2':
                    l2 = l2 + 1
                elif determine_layer(j) == '5':
                    l5 = l5 + 1
                elif determine_layer(j) == '6':
                    l6 = l6 + 1
        if count == 0:
            L1.append(0)
            L2.append(0)
            L5.append(0)
            L6.append(0)
        else:
            L1.append(l1 / count)  # probability of neuron in layer 1 of this region
            L2.append(l2 / count)
            L5.append(l5 / count)
            L6.append(l6 / count)
        count_sum.append(count)
        d1 = np.ones(l1) * 1
        d2 = np.ones(l2) * 2
        d5 = np.ones(l5) * 3
        d6 = np.ones(l6) * 4
        data = np.concatenate([d1, d2, d5, d6], axis=0)

        # polyfit data and find peak value

        # kde = seaborn.kdeplot(data=data, shade=True)
        # plt.xlim((0, 5))
        # plt.xticks([1, 2, 3, 4], ['L1', 'L2/3', 'L5', 'L6'])
        # x , y = kde.get_lines()[0].get_data()
        # peak_idx = signal.find_peaks(y)[0]
        # peak_value = x[peak_idx], y[peak_idx]
        # plt.plot(peak_value)

        print(count)
    distribution = [L1, L2, L5, L6]
    distribution = np.array(distribution)
    # # draw distribution
    # plt.figure(figsize=(12, 5))
    # width = 0
    # wid = 0.08
    # for i in range(11):
    #     plt.bar(X + width, distribution[:, i], width=wid, alpha=0.5)
    #     width = width + wid
    # plt.xticks(X + width / 2, ['L1', 'L2/3', 'L5', 'L6'])
    # plt.legend(region_sets)
    # plt.ylabel('probability')
    # plt.xlabel('layer')
    # plt.title('the probability distribution of neurons in all layers of each region ')
    # plt.show()

    return distribution, count_sum
def get_distribution_neuron_in_each_region_of_layer(region_sets, soma_regions):
    L_ratio = []
    L_number = []
    for layer in ['1', '2', '5', '6']:
        count = 0
        a1 = 0
        a2 = 0
        a3 = 0
        a4 = 0
        a5 = 0
        a6 = 0
        a7 = 0
        a8 = 0
        a9 = 0
        a10 = 0
        a11 = 0
        for k in soma_regions:
            if determine_layer(k) == layer:
                count = count + 1
                if remove_layer_information(k) == region_sets[0]:
                    a1 = a1 + 1
                elif remove_layer_information(k) == region_sets[1]:
                    a2 = a2 + 1
                elif remove_layer_information(k) == region_sets[2]:
                    a3 = a3 + 1
                elif remove_layer_information(k) == region_sets[3]:
                    a4 = a4 + 1
                elif remove_layer_information(k) == region_sets[4]:
                    a5 = a5 + 1
                elif remove_layer_information(k) == region_sets[5]:
                    a6 = a6 + 1
                elif remove_layer_information(k) == region_sets[6]:
                    a7 = a7 + 1
                elif remove_layer_information(k) == region_sets[7]:
                    a8 = a8 + 1
                elif remove_layer_information(k) == region_sets[8]:
                    a9 = a9 + 1
                elif remove_layer_information(k) == region_sets[9]:
                    a10 = a10 + 1
                elif remove_layer_information(k) == region_sets[10]:
                    a11 = a11 + 1
        union = np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11])
        if count == 0:
            union_ratio = np.zeros(np.shape(union))
        else:
            union_ratio = union / count
        L_number.append(union)
        L_ratio.append(union_ratio)
    return L_ratio, L_number



def get_distribution_neuron_in_each_layer(soma_regions):
    layer1 = 0
    layer2 = 0
    layer5 = 0
    layer6 = 0
    for i in soma_regions:
        if determine_layer(i) == '1':
            layer1 = layer1 + 1
        elif determine_layer(i) == '2':
            layer2 = layer2 + 1
        elif determine_layer(i) == '5':
            layer5 = layer5 + 1
        elif determine_layer(i) == '6':
            layer6 = layer6 + 1
    return np.array([layer1, layer2, layer5, layer6])


def draw_pie_image(sizes, labels, name, save_path):
    import matplotlib as mpt
    # mpt.rcParams['font.family'] = 'fangsong'
    # labels = 'c11', 'c12', 'c13', 'c14', 'c2'
    # sizes = [26.7, 20.5, 4.3, 2.8, 45.7]
    # labels = 'layer1', 'layer2/3', 'layer5', 'layer6'
    # sizes = [7.23, 34.43, 53.81, 4.53]
    explode = (0,)
    for i in range(len(sizes) - 1):
        explode = explode + (0,)
    plt.figure()
    plt.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90, pctdistance=0.8)
    plt.axis('equal')
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    # plt.savefig(save_path + '/ceil/' + name + '_pie.eps', format='pdf')
    plt.savefig(save_path + name + '_pie.eps', format='pdf')
    plt.close()


def draw_kdeplot(data, xlabels, save_path, name, line_color, shade_color):
    plt.figure()
    kde = seaborn.kdeplot(data=data, color=line_color)
    seaborn.despine(top=True, right=True)
    kde.fill_between(x=kde.get_lines()[0].get_xdata(), y1=kde.get_lines()[0].get_ydata(), color=shade_color)
    plt.xlim((0, 5))
    plt.xticks(range(1, len(xlabels) + 1), xlabels)
    save_fn = save_path + name + '.eps'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_fn, format='pdf')
    plt.close()


if __name__ == "__main__":
    # read neuron file projected to bilateral cortex region
    data_path = '/data/single_neuro_tracing/code/preprocess/new_results/IT_cortex_neuron_projected_to_bilateral_ceil.csv'
    data = pd.read_csv(data_path, usecols=[0, 1, 2, 3, 4])
    data = data.values
    neuron_id = data[:, 1]
    soma_region = data[:, 2]
    ipsi_termials = data[:, 3]
    contra_terminals = data[:, 4]

    ipsi_termials = remove_some_information(ipsi_termials)
    contra_terminals = remove_some_information(contra_terminals)
    soma_region_dl = [remove_layer_information(i) for i in soma_region]
    soma_region_dr = [determine_layer(i) for i in soma_region]
    # region_sets = ['FRP', 'MOs', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl', 'ORBm', 'ORBvl', 'AId', 'AIv']
    region_sets = ['FRP', 'ILA', 'AIv', 'ORBvl', 'ORBm', 'ORBl', 'ACAv', 'MOs', 'AId', 'ACAd', 'PL']  # PFC regions

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


    idx_specific_ipsi = np.where(np.array(s1) == 0)  # neuron whose specific ipsi conn is zero
    idx_specific_contra = np.where(np.array(s2) == 0)  # neuron whose specific contra conn is zero
    idx_specific_both = np.where(np.array(s3) == 0)  # neuron whose bilateral conn is zero

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
    c213 = []   # I = C
    c214 = []  # I is the best
    c22 = []  # B is the best
    c23 = []  # C is the best
    c231 = []  # C = B
    c232 = []  # C>B

    for j in range(len(neuron_id)):
        soma = soma_region[j]
        r1 = s1[j]  # get the ratio of single neuron in three types
        r2 = s2[j]
        r3 = s3[j]
        if r1 != 0 and r2 != 0 and r3 != 0:
            c2.append(j)
            if np.argmax(np.array([r1, r2, r3])) == 0:
                if len(np.where(np.array([r1, r2, r3]) == np.max(np.array([r1, r2, r3])))[0]) == 3:
                    c211.append(j)  # i = b = c
                elif len(np.where(np.array([r1, r2, r3]) == np.max(np.array([r1, r2, r3])))[0]) == 2:
                    if np.where(np.array([r1, r2, r3]) == np.max(np.array([r1, r2, r3])))[0][1] == 2:
                       c212.append(j)  # I =B
                    elif np.where(np.array([r1, r2, r3]) == np.max(np.array([r1, r2, r3])))[0][1] == 1:
                        c213.append(j)  # I = C
                else:
                    c214.append(j) # i is the best
                c21.append(j)
            elif np.argmax(np.array([r1, r2, r3])) == 2:
                c22.append(j)  # B is  the best
            elif np.argmax(np.array([r1, r2, r3])) == 1:
                if len(np.where(np.array([r1, r2, r3]) == np.max(np.array([r1, r2, r3])))[0]) > 1:
                    c231.append(j)  # C = B
                else:
                    c232.append(j)
                c23.append(j)
            '''
            if r1 >= r2:
                c21.append(j)
            else:
                c22.append(j)
            
            if np.argmax(np.array([r1, r2, r3])) == 1:
                c23.append(j)
                '''
        else:
            c1.append(j)
            if r1 != 0 and r2 == 0 and r3 != 0:
                c11.append(j)
                if r1 > r3:
                    c111.append(j)
                elif r1 < r3:
                    c112.append(j)
                else:
                    c113.append(j)
            elif r1 == 0 and r2 != 0 and r3 != 0:
                c12.append(j)
                if r2 > r3:
                    c121.append(j)
                elif r2 < r3:
                    c122.append(j)
                else:
                    c123.append(j)
            elif r1 == 0 and r2 == 0 and r3 != 0:
                c13.append(j)
            elif r1 != 0 and r2 != 0 and r3 == 0:
                c14.append(j)
                if r1 > r2:
                    c141.append(j)
                elif r1 < r2:
                    c142.append(j)
                else:
                    c143.append(j)
    # # draw the ratio of 5 type neuron in each layer and draw pie image [ layer level]
    save_path = './new_results/figures/'
    # L = get_distribution_neuron_in_each_layer(soma_region)
    # L_ = np.concatenate([np.ones(L[0]) * 1, np.ones(L[1]) * 2, np.ones(L[2]) * 3, np.ones(L[3]) * 4])
    # labels = ['L1', 'L2/3', 'L5', 'L6']
    # draw_kdeplot(L_, save_path=save_path, xlabels=labels, name='neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    # L_c2 = get_distribution_neuron_in_each_layer(soma_region[c2])
    # L_c2_ = np.concatenate([np.ones(L_c2[0]) * 1, np.ones(L_c2[1]) * 2, np.ones(L_c2[2]) * 3, np.ones(L_c2[3]) * 4])
    # draw_kdeplot(L_c2_, save_path=save_path, xlabels=labels, name='c2_neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    # L_c11 = get_distribution_neuron_in_each_layer(soma_region[c11])
    # L_c11_ = np.concatenate(
    #     [np.ones(L_c11[0]) * 1, np.ones(L_c11[1]) * 2, np.ones(L_c11[2]) * 3, np.ones(L_c11[3]) * 4])
    # draw_kdeplot(L_c11_, save_path=save_path, xlabels=labels, name='c11_neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    # L_c12 = get_distribution_neuron_in_each_layer(soma_region[c12])
    # L_c12_ = np.concatenate(
    #     [np.ones(L_c12[0]) * 1, np.ones(L_c12[1]) * 2, np.ones(L_c12[2]) * 3, np.ones(L_c12[3]) * 4])
    # draw_kdeplot(L_c12_, save_path=save_path, xlabels=labels, name='c12_neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    # L_c13 = get_distribution_neuron_in_each_layer(soma_region[c13])
    # L_c13_ = np.concatenate(
    #     [np.ones(L_c13[0]) * 1, np.ones(L_c13[1]) * 2, np.ones(L_c13[2]) * 3, np.ones(L_c13[3]) * 4])
    # draw_kdeplot(L_c13_, save_path=save_path, xlabels=labels, name='c13_neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    # L_c14 = get_distribution_neuron_in_each_layer(soma_region[c14])
    # L_c14_ = np.concatenate(
    #     [np.ones(L_c14[0]) * 1, np.ones(L_c14[1]) * 2, np.ones(L_c14[2]) * 3, np.ones(L_c14[3]) * 4])
    # draw_kdeplot(L_c14_, save_path=save_path, xlabels=labels, name='c14_neuron_in_layer', line_color="#000000",
    #              shade_color="#B91F00")
    #
    # summary_layer = np.array([L, L_c11, L_c12, L_c13, L_c14, L_c2]).transpose()
    # df = pd.DataFrame(
    #     {'names': ['all', 'c11', 'c12', 'c13', 'c14', 'c2'], 'L1': summary_layer[0], 'L2_3': summary_layer[1],
    #      'L5': summary_layer[2], 'L6': summary_layer[3]})
    # df.to_csv(save_path + '/layer_distribution.csv')
    # # draw the percentage of subtype in all neurons
    # sizes_all = [np.round(len(c11) / len(neuron_id), decimals=3) * 100,
    #              np.round(len(c12) / len(neuron_id), decimals=3) * 100,
    #              np.round(len(c13) / len(neuron_id), decimals=3) * 100,
    #              np.round(len(c14) / len(neuron_id), decimals=3) * 100,
    #              np.round(len(c2) / len(neuron_id), decimals=3) * 100]
    #
    # draw_pie_image(sizes_all, ['c11', 'c12', 'c13', 'c14', 'c2'], 'L_all', save_path=save_path)
    #
    # # draw the percentage of subtype in each layer neurons
    # Layer = ['layer1', 'layer2-3', 'layer5', 'layer6']
    # for i in range(len(Layer)):
    #     sizes = [np.round(L_c11[i] / L[i], decimals=3) * 100, np.round(L_c12[i] / L[i], decimals=3) * 100,
    #              np.round(L_c13[i] / L[i], decimals=3) * 100, np.round(L_c14[i] / L[i], decimals=3) * 100,
    #              np.round(L_c2[i] / L[i], decimals=3) * 100]
    #     labels = 'c11', 'c12', 'c13', 'c14', 'c2'
    #     draw_pie_image(sizes, labels, Layer[i], save_path)   # layer
    #
    # # neuron number in each layer of region
    # number_d = get_number_distribution_of_neuron_in_layer_of_region(region_sets, soma_region)
    # df = pd.DataFrame(
    #     {'region_name': region_sets, 'L1': number_d[0], 'L2_3': number_d[1], 'L5': number_d[2], 'L6': number_d[3]})
    # df.to_csv(save_path + '/distribution_in_layer_of_region.csv')
    # # subtype neuron number in each region
    number_region = get_number_of_different_type_neuron_in_different_region(region_sets, soma_region)
    number_region = np.array(number_region).transpose()
    # df = pd.DataFrame(
    #     {'region_name': region_sets, 'c11': number_region[0], 'c12': number_region[1], 'c13': number_region[2],
    #      'c14': number_region[3], 'c2': number_region[4]})
    # # df.to_csv(save_path + '/ceil/subtype_in_region.csv')
    # df.to_csv(save_path + '/subtype_in_region.csv')
    # draw the percentage of subtype in each PFC region (matrix)
    number_region = number_region.transpose()
    sum_region = np.sum(number_region, axis=1)
    percentage_subtype_in_region = np.round(number_region / sum_region.reshape((-1,1)), decimals=2)
    percentage_subtype_in_region = percentage_subtype_in_region.transpose()
    # df = pd.DataFrame(
    #     {'region_name': region_sets, 'c11': number_region.transpose()[0], 'c12': number_region.transpose()[1],
    #      'c13': number_region.transpose()[2], 'c14': number_region.transpose()[3],
    #      'c2': number_region.transpose()[4]})
    # # df.to_csv(save_path + '/ceil/percentage_of_subtype_in_region.csv')
    # df.to_csv(save_path + '/percentage_of_subtype_in_region.csv')
    mymap2 = define_colormap4("tab20c", "tab20b")
    fig = plt.gcf()
    fig.set_size_inches((10, 5), forward=False)
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.labeltop'] = True
    seaborn.heatmap(percentage_subtype_in_region, vmin=0, vmax=0.5, cmap=mymap2, annot=percentage_subtype_in_region)  #region
    ax = plt.gca()
    ax.set_yticklabels(labels=['IB', 'BC', 'B', 'IC', 'IBC'], rotation=0)
    ax.set_xticklabels(labels=region_sets, rotation=90)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_path + '/percentage_of_neuron_in_region.eps', format='pdf')
    plt.close()
    #
    # # get the distribution of one neuron type in one region in different layer and the total number of the neuron type in different regions
    # # # draw the percentage of region neurons in all neurons
    [a, count_in_region] = get_distribution_of_neuron_in_layer_of_region(region_sets, soma_region)
    percentage_region = np.round(np.array(count_in_region) / len(neuron_id), decimals=3)
    sorted_idx = np.argsort(percentage_region)
    plt.figure()
    plt.barh(np.array(region_sets)[sorted_idx], percentage_region[sorted_idx], color='#B071BF')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.savefig(save_path + '/neuron_in_region.eps', format='pdf')
    plt.close()

    