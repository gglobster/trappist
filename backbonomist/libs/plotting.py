#import matplotlib
#matplotlib.use('Agg')
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import AxesGrid#, make_axes_locatable, Size
import numpy as np

def plot_ctg_stats(ctg_cats, fname, ctg_thresholds):
    """Plot distribution statistics of contigs per genome."""
    small_ctgs, mid_ctgs, big_ctgs = ctg_cats
    kb1, kb2, kb3 = ctg_thresholds
    # make bins
    factor = 10
    small_bins = np.arange(0, kb1, kb1/factor)
    mid_bins = np.arange(kb1, kb2, (kb2-kb1)/factor)
    big_bins = np.arange(kb2, kb3, (kb3-kb2)/factor)
    # plot histograms
    fig = plt.figure(figsize=(16,6))
    if not len(small_ctgs) == 0:
        small_x = fig.add_subplot(131)
        small_x.hist(small_ctgs, small_bins, rwidth=kb1/factor)
        small_x.set_title("length < "+str(kb1)+" kb")
        small_x.set_ylabel("number of contigs")
    if not len(mid_ctgs) == 0:
        mid_x = fig.add_subplot(132)
        mid_x.hist(mid_ctgs, mid_bins, rwidth=(kb2-kb1)/factor)
        mid_x.set_title(str(kb1)+" kb < length < "+str(kb2)+" kb")
        mid_x.set_ylabel("number of contigs")
    if not len(big_ctgs) == 0:
        big_x = fig.add_subplot(133)
        big_x.hist(big_ctgs, big_bins, rwidth=(kb3-kb2)/factor)
        big_x.set_title("length > "+str(kb2)+" kb")
        big_x.set_ylabel("number of contigs")
    plt.savefig(fname)
    plt.clf()

def hits_heatmap_multi(ref_n, segs, g_names, contigs, scores, imgfile):
    """Combine matches heatmaps in subplots."""
    g_names.reverse()
    contigs.reverse()
    scores.reverse()
    ctg_count = 0
    for item in contigs:
        ctg_count += len(item)
    fig = plt.figure(figsize=(len(segs)/2+1, ctg_count/2+3))
    grid = AxesGrid(fig, 111, nrows_ncols=(len(g_names), 1), axes_pad=0.4,
                    cbar_location="top", cbar_mode="single", cbar_size=0.1)
    for i in range(len(g_names)):
        mod_score = [-score for score in scores[i]]
        score_RA = np.array(mod_score)
        hmap = grid[i].pcolor(score_RA, cmap='hot', vmin=-1, vmax=0)
        grid[i].xaxis.set_major_locator(MaxNLocator(len(segs)))
        grid[i].yaxis.set_major_locator(MaxNLocator(len(scores[i])))
        grid.cbar_axes[i].colorbar(hmap)
        grid.cbar_axes[i].set_xticks([-1, -0.5, 0])
        grid.cbar_axes[i].set_xticklabels(['High', 'Medium', 'Low'])
        grid[i].set_xticklabels(segs, size='small')
        grid[i].set_xlabel(ref_n+" reference segments", size='small')
        grid[i].set_title(g_names[i], size='small')
        grid[i].set_yticklabels("")
        labels = grid[i].get_xticklabels()
        for label in labels:
            label.set_rotation(30)
        y_index = 0.5
        for contig in contigs[i]:
            grid[i].text(-0.2, y_index, contig, size='small',
                         horizontalalignment='right',
                         verticalalignment='center',)
            y_index +=1
    plt.savefig(imgfile)
    plt.clf()

#def hits_heatmap(title, segs, contigs, scores, imgfile):
#    """Graph matches as a heatmap.
#
#    Deprecated in favor of the multi-plot version.
#
#    """
#    fig = plt.figure()
#    h = [Size.Fixed(1.5), Size.Fixed(15)]
#    v = [Size.Fixed(2.5), Size.Fixed(5.)]
#    ax = plt.subplot(111)
#    hmap = ax.pcolor(scores, cmap='hot', vmin=0, vmax=1)
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("bottom", size="5%", pad=1)
#    cbar = plt.colorbar(hmap, cax=cax, ticks=[0, 0.5, 1],
#                    orientation='horizontal')
#    cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])
#    ax.xaxis.set_major_locator(MaxNLocator(len(segs)))
#    ax.yaxis.set_major_locator(MaxNLocator(len(contigs)))
#    ax.set_xticklabels(segs, size='small')
#    ax.set_yticklabels("")
#    labels = ax.get_xticklabels()
#    for label in labels:
#        label.set_rotation(30)
#    y_index = 0.5
#    for contig in contigs:
#        ax.text(-0.2, y_index, contig, size='small',
#                horizontalalignment='right',
#                verticalalignment='center',)
#        y_index +=1
#    ax.set_xlabel("Reference segments", size='small')
#    ax.set_title(title)
#    plt.savefig(imgfile)
#    plt.clf()