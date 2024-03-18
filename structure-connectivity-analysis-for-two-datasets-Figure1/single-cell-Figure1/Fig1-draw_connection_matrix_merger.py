import scipy.io as sio
from matplotlib import pyplot as plt
import numpy as np
from plot_colormap import define_colormap4
c_map = define_colormap4("tab20c", "tab20b")



def draw_figure(matrix1, matrix2, sour, tar, layer):
    import matplotlib.gridspec as gridspec
    import seaborn
    import matplotlib.colors as colors

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

    matrix1 = matrix1 + 1e-10
    matrix2 = matrix2 + 1e-10
    print(np.min(matrix2[np.where(matrix2 - 1e-10)]))
    vmin = 0
    vmax = 10 ** 0
    # norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    seaborn.heatmap(matrix1, ax=heatmap_ax1, cbar=False, cmap=c_map, vmin=vmin, vmax=vmax)
    seaborn.heatmap(matrix2, ax=heatmap_ax2, cbar_ax=cbar_ax, yticklabels=False, cmap=c_map,
                    vmin=vmin, vmax=vmax, cbar_kws={'label': 'connection strength'})

    heatmap_ax1.set_xticklabels(labels=tar, rotation=90, size=5)
    heatmap_ax1.set_yticklabels(labels=sour, rotation=0, size=5)
    heatmap_ax2.set_xticklabels(labels=tar, rotation=90, size=5)
    heatmap_ax1.tick_params(left=False)
    heatmap_ax1.tick_params(top=False)
    heatmap_ax2.tick_params(top=False)

    # image_path = ''.join(['./IT_results_new/ceil/L', layer, '_conn_reorder_binary_20230424.eps'])
    image_path = ''.join(['./cre_results', layer, '_conn_reorder_merger_binary_20230616.eps'])
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    fig.savefig(image_path, format='pdf')
    plt.close(fig)
    print("Save image to %s" % image_path)
if __name__ == "__main__":
    data_path = '/data/single_neuro_tracing/code/preprocess/new_results/matrix/'
    file1_path = data_path + 'ipsilateral_conn.mat'
    file2_path = data_path + 'contralateral_conn.mat'
    ipsi = sio.loadmat(file1_path)['conn']
    contra = sio.loadmat(file2_path)['conn']

    regions_set = ['FRP', 'ACAd', 'ACAv', 'PL', 'ILA', 'ORBl', 'ORBm', 'ORBvl', 'AId', 'AIv', 'MOs']

    # exchange all cortex region order according to thr following file
    cortex_region_name_module_path = '/data/single_neuro_tracing/code/preprocess/IT_results/cortex_name_module_cell.txt'
    cortex_region_name_module = []
    with open(cortex_region_name_module_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            cortex_region_name_module.append(line)



    idx_ssp = [i for i in range(len(cortex_region_name_module)) if 'SSp' in cortex_region_name_module[i]]

    cortex_region_name_ = cortex_region_name_module[0:idx_ssp[0]]
    cortex_region_name_.extend(['SSp'])
    cortex_region_name_.extend(cortex_region_name_module[idx_ssp[6]+1:43])

    mat1_ssp = ipsi[:, idx_ssp]
    vec1_ssp = np.mean(mat1_ssp, axis=1)
    mat2_ssp = contra[:, idx_ssp]
    vec2_ssp = np.mean(mat2_ssp, axis=1)
    ipsi_ = np.zeros((11, 37))
    contra_ = np.zeros((11, 37))
    ipsi_[:, 0:idx_ssp[0]] = ipsi[:, 0:idx_ssp[0]]
    contra_[:, 0:idx_ssp[0]] = contra[:, 0:idx_ssp[0]]
    ipsi_[:,idx_ssp[0]+1:37] = ipsi[:, idx_ssp[6]+1:43]
    contra_[:, idx_ssp[0]+1:37] = contra[:, idx_ssp[6]+1:43]
    ipsi_[:, idx_ssp[0]] = vec1_ssp
    contra_[:, idx_ssp[0]] = vec2_ssp

    ipsi_ = np.where(ipsi_ > 0, 1, 0)
    contra_ = np.where(contra_ > 0, 1, 0)



    draw_figure(ipsi_, contra_, regions_set, cortex_region_name_, 'all')

