import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import seaborn
from plot_colormap import define_colormap, define_colormap2, define_colormap4

c_map = define_colormap4("tab20c", "tab20b")
c_map1 = define_colormap2("tab20b")


def draw_figure(ipsi, contra, save_dir, current, cortex_region_names, PFC):
    fig = plt.gcf()
    fig.set_size_inches((21.5, 3), forward=False)
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.labeltop'] = True

    gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=(0.80, 0.80, 0.02), wspace=0.01)
    heatmap_ax1 = fig.add_subplot(gs[0])
    heatmap_ax2 = fig.add_subplot(gs[1])
    cbar_ax = fig.add_subplot(gs[2])
    ipsi = ipsi + 1e-10
    contra = contra + 1e-10

    # vmin = 10 ** -1.5
    # vmax = 10 ** 1
    vmin = 0
    vmax = 1
    # norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    seaborn.heatmap(ipsi, ax=heatmap_ax1, cbar=False, cmap=c_map, vmin=vmin,
                    vmax=vmax, alpha=1.0, cbar_kws={'label': 'normalized connection strength'})
    seaborn.heatmap(contra, ax=heatmap_ax2, cbar_ax=cbar_ax, yticklabels=False, cmap=c_map,
                     vmin=vmin, vmax=vmax, alpha=1.0, cbar_kws={'label': 'normalized connection strength'})

    heatmap_ax1.set_xticklabels(labels=cortex_region_names, rotation=90, size=5)
    heatmap_ax1.set_yticklabels(labels=PFC, rotation=0, size=5)
    heatmap_ax2.set_xticklabels(labels=cortex_region_names, rotation=90, size=5)
    heatmap_ax1.tick_params(left=False)
    heatmap_ax1.tick_params(top=False)
    heatmap_ax2.tick_params(top=False)
    image_path = save_dir + '/' + current + '/' + 'SC_PFC_reorder_binary_20230424' + ".eps"
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    fig.savefig(image_path, format='pdf')
    plt.close()
    print("Save image to %s" % image_path)


if __name__ == "__main__":

    data_path = '../cre_results'
    #PFC_regions = ['FRP', 'MOs', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl', 'ORBm', 'ORBvl', 'AId', 'AIv']
    PFC_regions = ['FRP', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl', 'ORBm', 'ORBvl', 'AId', 'AIv', 'MOs']
    f = open('../cortex_region_names.txt')
    cortex_regions_name = []
    for i in f.readlines():
        cortex_regions_name.append(i)
    cortex_regions_name = [i.split('\n')[0] for i in cortex_regions_name]

    # # PFC REGION VS cortex region
    # idx = []
    # for j in range(len(PFC_regions)):
    #     region = PFC_regions[j]
    #     for k in range(len(cortex_regions_name)):
    #         if region == cortex_regions_name[k]:
    #             idx.append(k)
    # idx1 = list(set([i for i in range(43)]) - set(idx))
    # idx2 = []
    # idx2.extend(idx)
    # idx2.extend(idx1)
    # cortex_regions_name_order = list(np.array(cortex_regions_name)[idx2])

    idx1 = []
    for j in range(len(PFC_regions)):
        region = PFC_regions[j]
        for k in range(len(cortex_regions_name)):
            if region == cortex_regions_name[k]:
                idx1.append(k)
    # exchange all cortex region order according to thr following file
    cortex_region_name_module_path = '/single_neuro_tracing/code/preprocess/IT_results/cortex_name_module_cell.txt'
    cortex_region_name_module = []
    with open(cortex_region_name_module_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            cortex_region_name_module.append(line)
    idx2 = []
    for i in range(len(cortex_region_name_module)):
        print(i)
        for j in range(len(cortex_regions_name)):
            if cortex_regions_name[j] == cortex_region_name_module[i]:
                idx2.append(int(j))
                break


    wt_path = '../results_new'
    thresold = 10 ** -1.5  
    wt_r2r = np.load(''.join([wt_path, '/', 'R2R_SC_strength.npy']))
    wt_r2l = np.load(''.join([wt_path, '/', 'R2L_SC_strength.npy']))
    wt_r2r[wt_r2r < thresold] = 0
    wt_r2l[wt_r2l < thresold] = 0
    # wt_r2r_pfc = wt_r2r[idx]
    # wt_r2l_pfc = wt_r2l[idx]
    wt_r2r = [np.array(k)[idx2] for k in wt_r2r]
    wt_r2r_pfc = np.array(wt_r2r)[idx1]
    wt_r2l = [np.array(k)[idx2] for k in wt_r2l]
    wt_r2l_pfc = np.array(wt_r2l)[idx1]

    sio.savemat('./R2L_strength_PFC.mat', {'conn': wt_r2l_pfc})
    sio.savemat('./R2R_strength_PFC.mat', {'conn': wt_r2r_pfc})

    wt_P1 = (len(np.where(wt_r2r_pfc)[0]) - 11) / (11 * 43 - 11)  # IPSI SPARSE
    wt_P2 = len(np.where(wt_r2l_pfc)[0]) / (11 * 43)
    wt_r2r_pfc = np.where(wt_r2r_pfc > 0, 1, 0)
    wt_r2l_pfc = np.where(wt_r2l_pfc > 0, 1, 0)
    save_path = '../cre_results'
    draw_figure(wt_r2r_pfc, wt_r2l_pfc, save_path, 'WT', cortex_region_name_module, PFC_regions)

    typelist = ['Emx1', 'layer2-3', 'layer4', 'layer5', 'layer6']
    percentage_ipsi = [wt_P1]
    percentage_contra = [wt_P2]
    for tp in typelist:
        save_dir = ''.join([data_path, '/', tp, '/'])
        file1_path = ''.join([data_path, '/', tp, '/', 'R2R_SC_strength.npy'])
        file2_path = ''.join([data_path, '/', tp, '/', 'R2L_SC_strength.npy'])
        R2R = np.load(file1_path)
        R2L = np.load(file2_path)
        R2R[R2R < thresold] = 0
        R2L[R2L < thresold] = 0
        # R2R_pfc = R2R[idx]
        # R2L_pfc = R2L[idx]
        R2R = [np.array(k)[idx2] for k in R2R]
        R2R_pfc = np.array(R2R)[idx1]
        R2L = [np.array(k)[idx2] for k in R2L]
        R2L_pfc = np.array(R2L)[idx1]
        p1 = (len(np.where(R2R_pfc)[0]) - 11) / (43 * 11 - 11)  # IPSI SPARSE
        p2 = len(np.where(R2L_pfc)[0]) / (43 * 11)  # CONTRA SPARSE
        percentage_ipsi.append(p1)
        percentage_contra.append(p2)
        R2R_pfc = np.where(R2R_pfc > 0, 1, 0)
        R2L_pfc = np.where(R2L_pfc > 0, 1, 0)
        draw_figure(R2R_pfc, R2L_pfc, data_path, tp, cortex_region_name_module, PFC_regions)

    print(1)
