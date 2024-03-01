import scipy.io as sio
import seaborn
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_1samp
from statsmodels.stats import multitest

if __name__ == "__main__":
    # data1_path = '/wield-field-data/new_analysis/fc_strength/rem_fc_strength_type5.mat'
    data1_path = '/wield-field-data/new_analysis/FC_strength/wake-rem_difference_FC_type3.mat'
    data1 = sio.loadmat(data1_path)['FC']
    data1 = np.abs(data1)
    fc_3 = ['ipsi', 'hete', 'homo']
    fc_strength = []
    types = []
    m = 0
    p_values = []
    for i in range(3):
        s = data1[:, i]
        t_stat, p_value = ttest_1samp(s, 0)
        print("t_stst: {0}, p_value: {1}".format(t_stat, p_value))
        p_values.append(p_value)
        s = list(s)
        fc_strength.extend(s)
        types.extend([fc_3[m]] * len(s))
        m = m + 1
    reject, corrected_p_values, _, _ = multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')
    
    a = {'difference fc strength': fc_strength, 'type': types}
    df = pd.DataFrame(a)
    #
    plt.figure()
    # # # draw violin plot
    seaborn.violinplot(y='difference fc strength', x='type', data=df)
    # #
    # # # # # stacked scatter plot
    # # # seaborn.stripplot(y='fc strength', x='type', data=df, size=4, color="black")
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig('/wield-field-data/violin_plot/results/wake-rem_difference_fc_type3_abs.eps', format='pdf')
    plt.close()
    # #
    tukey_result = pairwise_tukeyhsd(a['difference fc strength'], a['type'])
    print(tukey_result)

#    # data1_path =' /wield-field-data/new_analysis/FC_strength/rem_FC_strength_type5.mat'
#    data1_path = '/wield-field-data/new_analysis/FC_strength/rem-nrem_difference_FC_type3_ipsi_hete.mat'
#    data1 = sio.loadmat(data1_path)['FC']
#    data1 = np.abs(data1)
#    FC_2 = ['ipsi', 'hete']
#    FC_strength = []
#    types = []
#    m = 0
#    p_values = []
#    for i in range(2):
#        s = data1[:, i]
#        t_stat, p_value = ttest_1samp(s, 0)
#        print("t_stst: {0}, p_value: {1}".format(t_stat, p_value))
#        p_values.append(p_value)
#        s = list(s)
#        FC_strength.extend(s)
#        types.extend([FC_2[m]] * len(s))
#        m = m + 1
#
#    reject, corrected_p_values, _, _ = multitest.multipletests(p_values, alpha=0.05, method='fdr_bh')
#    a = {'difference FC strength': FC_strength, 'type': types}
#    df = pd.DataFrame(a)
#
#    plt.figure()
#    # draw violin plot
#    seaborn.violinplot(y='difference fc strength', x='type', data=df)
#    #
#    # # # # stacked scatter plot
#    #seaborn.stripplot(y='fc strength', x='type', data=df, size=4, color="black")
#    plt.rcParams['ps.fonttype'] = 42
#    plt.rcParams['pdf.fonttype'] = 42
#    plt.savefig('/wield-field-data/violin_plot/results/rem-nrem_difference_fc_type3_ipsi_hete_abs.eps', format='pdf')
#    plt.close()
#
#    tukey_result = pairwise_tukeyhsd(a['difference FC strength'], a['type'])
#    print(tukey_result)
