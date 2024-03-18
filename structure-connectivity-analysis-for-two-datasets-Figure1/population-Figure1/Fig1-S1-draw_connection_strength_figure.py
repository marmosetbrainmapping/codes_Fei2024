import os
import time
import seaborn
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from utils.structure_mask import StructureMask
from utils.experiment import Experiment
from utils.constants import *
import allensdk.core.json_utilities as ju
import multiprocessing
import scipy.io as sio
from model.global_model import GlobalModel
import pandas as pd
from cre_results import find_experiments_related_layer
from model.cortex_model import CortexModel
from model.global_model import GlobalModel
import plot_colormap
c_map1 = plot_colormap.define_colormap("tab20c", "tab20b")

def draw_connection_figure(mat, tp, cortex_region_names):
    # for log(x + e)
    epsilon = 1e-10
    mat += epsilon
    ipsilateral_mat = mat[:, 0:43]
    contralateral_mat = mat[:, 43:86]
    fig = plt.gcf()
    fig.set_size_inches((20, 10), forward=False)
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.labeltop'] = True

    gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=(0.80, 0.80, 0.02), wspace=0.01)
    heatmap_ax1 = fig.add_subplot(gs[0])
    heatmap_ax2 = fig.add_subplot(gs[1])
    cbar_ax = fig.add_subplot(gs[2])

    vmin = 10 ** -1.5
    vmax = 10 ** 1

    norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    seaborn.heatmap(ipsilateral_mat, ax=heatmap_ax1, cbar=False, cmap=c_map1, norm=norm, vmin=vmin,
                    vmax=vmax, alpha=1)
    seaborn.heatmap(contralateral_mat, ax=heatmap_ax2, cbar_ax=cbar_ax, yticklabels=False, cmap=c_map1,
                    norm=norm, vmin=vmin, vmax=vmax, alpha=1, cbar_kws={'label': 'normalized connection strength'})
    heatmap_ax1.set_xticklabels(labels=cortex_region_names, rotation=90)
    heatmap_ax1.set_yticklabels(labels=cortex_region_names, rotation=0)
    heatmap_ax2.set_xticklabels(labels=cortex_region_names, rotation=90)
    heatmap_ax1.tick_params(left=False)
    heatmap_ax1.tick_params(top=False)
    heatmap_ax2.tick_params(top=False)

    image_path = './connection_figures/' + tp + "_connection_strength_NEW.svg"
    plt.rcParams['svg.fonttype'] = 'none'
    fig.savefig(image_path)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    fig.savefig('./connection_figures/' + tp + "_connection_strength_NEW.eps", format='pdf')

    # print("Save image to %s" % image_path)
    plt.close()


def draw_connection_difference_figure(mat, tp, cortex_region_names):
    # for log(x + e)
    epsilon = 1e-10
    mat[mat < (10 ** -1.5)] = 0
    mat += epsilon
    ipsilateral_mat = mat[:, 0:43]
    contralateral_mat = mat[:, 43:86]
    difference = ipsilateral_mat - contralateral_mat
    x = difference[np.where(difference<0)]
    print('min difference:')
    print(np.min(difference))
    print('max difference:')
    print(np.max(difference))
    fig = plt.gcf()
    fig.set_size_inches((10, 10), forward=False)
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.labeltop'] = True

    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=(0.80, 0.02), wspace=0.01)
    heatmap_ax1 = fig.add_subplot(gs[0])
    cbar_ax = fig.add_subplot(gs[1])

    vmin = -3
    vmax = 3

    # norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    seaborn.heatmap(difference, ax=heatmap_ax1, cbar_ax=cbar_ax, cmap=c_map,
                     cbar=cbar_ax, vmin=vmin, vmax=vmax, alpha=1, cbar_kws={'label': 'normalized connection strength'})
    heatmap_ax1.set_xticklabels(labels=cortex_region_names, rotation=90)
    heatmap_ax1.set_yticklabels(labels=cortex_region_names, rotation=0)

    heatmap_ax1.tick_params(left=False)
    heatmap_ax1.tick_params(top=False)



    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    fig.savefig('./connection_figures/' + tp + "_connection_strength_difference_thresh.eps", format='pdf')

    plt.close()


if __name__ == "__main__":
    os.chdir('../')
    global_model = GlobalModel()
    cortex_model = CortexModel(global_model)
    cortex_region_name = cortex_model.cortex_region_names

    cortex_region_name_module_path = '../cortex_name_module.txt'
    cortex_region_name_module = []
    with open(cortex_region_name_module_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            cortex_region_name_module.append(line)
    idx = []
    idx1 = []
    idx2 = []
    for i in range(len(cortex_region_name_module)):
        print(i)
        for j in range(len(cortex_region_name)):
            if cortex_region_name[j] == cortex_region_name_module[i]:
                idx.append(int(j))
                idx1.append(int(j + 43))
                break
    idx2.extend(idx)
    idx2.extend(idx1)

    data1_path = '../results_new/mean_mean-202211202030_strength.npy'
    data1 = np.load(data1_path)

    data1 = [np.array(k)[idx2] for k in data1]
    data1 = np.array(data1)[idx]
    print(np.max(data1))

    draw_connection_figure(data1, 'WT', cortex_region_names=cortex_region_name_module)
    # draw_connection_difference_figure(data1, 'WT', cortex_region_names=cortex_region_name_module)
    typelist = ['Emx1', 'layer2-3', 'layer4', 'layer5', 'layer6']
    data_path = '../cre_results/'

    for tp in typelist:
        save_dir = ''.join([data_path, '/', tp, '/'])
        file1_path = ''.join([data_path, '/', tp, '/', 'R2R_SC_strength.npy'])
        file2_path = ''.join([data_path, '/', tp, '/', 'R2L_SC_strength.npy'])
        R2R = np.load(file1_path)
        R2L = np.load(file2_path)
        R2R = [np.array(k)[idx] for k in R2R]
        R2R = np.array(R2R)[idx]
        R2L = [np.array(k)[idx] for k in R2L]
        R2L = np.array(R2L)[idx]
        both = np.concatenate([R2R, R2L], axis=1)
        print(np.max(both))
        draw_connection_figure(both, tp, cortex_region_names=cortex_region_name_module)
        draw_connection_difference_figure(both, tp, cortex_region_names=cortex_region_name_module)
