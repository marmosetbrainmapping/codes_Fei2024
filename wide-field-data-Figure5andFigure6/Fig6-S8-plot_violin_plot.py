import scipy.io as sio
import seaborn
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd

if __name__ == "__main__":
    data1_path = '/wide-field-data/new_analysis/FC_strength/rem_FC_strength_type5.mat'
    data1 = sio.loadmat(data1_path)['FC']
    FC_3 = ['ipsi', 'hete', 'homo']
    FC_strength = []
    types = []
    m = 0
    for i in range(3):
        s = data1[:, i]
        s = list(s)
        FC_strength.extend(s)
        types.extend([FC_3[m]] * len(s))
        m = m + 1

    a = {'FC strength': FC_strength, 'type': types}
    df = pd.DataFrame(a)

    plt.figure()
    # draw violin plot
    seaborn.violinplot(y='FC strength', x='type', data=df)

    # # # stacked scatter plot
    # seaborn.stripplot(y='FC strength', x='type', data=df, size=4, color="black")
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig('./results/rem_FC_strength_type5.eps', format='pdf')
    plt.close()

    tukey_result = pairwise_tukeyhsd(a['FC strength'], a['type'])
    print(tukey_result)

    # data2_path = '/wide-field-data/new_analysis/condition_varibility/saliences_type5.mat'
    # data2 = sio.loadmat(data2_path)['saliences_FC']
    # FC_3 = ['ipsi', 'hete', 'homo']
    # saliences = []
    # types = []
    # m = 0
    # for i in range(3):
    #     s = data2[:, i]
    #     s = list(s)
    #     saliences.extend(s)
    #     types.extend([FC_3[m]] * len(s))
    #     m = m + 1
    #
    # a = {'salience': saliences, 'type': types}
    # df = pd.DataFrame(a)
    #
    # plt.figure()
    # # draw violin plot
    # seaborn.violinplot(y='salience', x='type', data=df)
    #
    # # # # stacked scatter plot
    # # seaborn.stripplot(y='salience', x='type', data=df, size=4, color="black")
    # plt.rcParams['ps.fonttype'] = 42
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.savefig('./results/saliences_type5.eps', format='pdf')
    # plt.close()
    # tukey_result = pairwise_tukeyhsd(a['salience'], a['type'])
    # print(tukey_result)

    # data3_path = '/wide-field-data/new_analysis/FC_strength/rem_FC_strength_type5_ipsi_hete.mat'
    # data3 = sio.loadmat(data3_path)['FC']
    # FC_2 = ['ipsi', 'hete']
    # FC_strength = []
    # types = []
    # m = 0
    # for i in range(2):
    #     s = data3[:, i]
    #     s = list(s)
    #     FC_strength.extend(s)
    #     types.extend([FC_2[m]] * len(s))
    #     m = m + 1
    #
    # a = {'FC strength': FC_strength, 'type': types}
    # df = pd.DataFrame(a)

    # plt.figure()
    # # draw violin plot
    # seaborn.violinplot(y='FC strength', x='type', data=df)
    #
    # # # # stacked scatter plot
    # # seaborn.stripplot(y='FC strength', x='type', data=df, size=4, color="black")
    # plt.rcParams['ps.fonttype'] = 42
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.savefig('./results/wake_FC_strength_type5_ipsi_hete.eps', format='pdf')
    # plt.close()

    # tukey_result = pairwise_tukeyhsd(a['FC strength'], a['type'])
    # print(tukey_result)

    # data4_path = '/mnt/Data16Tb/Data/feiyao/wield-field-data/new_analysis/condition_varibility/saliences_type5_ipsi_hete.mat'
    # data4 = sio.loadmat(data4_path)['saliences_FC']
    # FC_2 = ['ipsi', 'hete']
    # saliences = []
    # types = []
    # m = 0
    # for i in range(2):
    #     s = data4[:, i]
    #     s = list(s)
    #     saliences.extend(s)
    #     types.extend([FC_2[m]] * len(s))
    #     m = m + 1
    #
    # a = {'salience': saliences, 'type': types}
    # df = pd.DataFrame(a)
    #
    # plt.figure()
    # # draw violin plot
    # seaborn.violinplot(y='salience', x='type', data=df)
    #
    # # # # stacked scatter plot
    # # seaborn.stripplot(y='salience', x='type', data=df, size=4, color="black")
    # plt.rcParams['ps.fonttype'] = 42
    # plt.rcParams['pdf.fonttype'] = 42
    # plt.savefig('./results/saliences_type5_ipsi_hete.eps', format='pdf')
    # plt.close()
    # tukey_result = pairwise_tukeyhsd(a['salience'], a['type'])
    # print(tukey_result)
    #
    

    print(1)