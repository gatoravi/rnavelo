import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import velocyto as vcy
import numpy as np
import sys

def usage():
    print(" ".join(["python3", sys.argv[0], "loom_file"]))

def main():
    #vlm = vcy.VelocytoLoom("velocyto/cellsorted_wgEncodeCshlLongRnaSeqGm12878CellTotalAlnRep1_withtags_sorted_sampled_385ZR.loom")
    if len(sys.argv) < 2:
        return usage()
    vlm = vcy.VelocytoLoom(sys.argv[1])
    #vlm.set_clusters(vlm.ca["ClusterName"])
    vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
    #vlm.filter_genes(by_detection_levels=True)
    #vlm.normalize("S", size=True, log=False)
    #vlm.normalize("U", size=True,  log=False)
    vlm._normalize_S(relative_size=vlm.S.sum(0),
                     target_size=vlm.S.sum(0).mean())
    vlm._normalize_U(relative_size=vlm.U.sum(0),
                     target_size=vlm.U.sum(0).mean())
    #vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)
    #vlm.fit_gammas()
    #vlm.plot_phase_portraits(["SOX2"])
    #plt.figure(figsize=(10,10))
    vlm.plot_fractions(save2file = sys.argv[1].split(".loom")[0] + ".png")
    # Plot the principal components
    pc_plot(vlm)
    # Calculate the gamma fit?
    calculate_gammas(vlm)
    # Calculate and plot velocity
    calculate_velocity(vlm)

def calculate_velocity(vlm):
    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift()
    vlm.extrapolate_cell_at_t(delta_t=1)
    vlm.estimate_transition_prob(hidim="Sx_sz", embed="pcs", transform="log", psc=1,
                             n_neighbors=1, sampled_fraction=1)
    vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=False)
    vlm.colorandum = "red"
    vlm.calculate_grid_arrows(smooth=0.5, steps=(25, 25, 25, 25), n_neighbors=1)
    plt.figure(None, (9,9))
    #vlm.plot_arrows_embedding(choice=21, scatter_kwargs={"alpha":1, "lw":0.01})
    vlm.plot_grid_arrows(scatter_kwargs_dict={"c": "red", "alpha":0.7, "lw":0.7, "edgecolor":"0.4", "s":70, "rasterized":True}, angles='xy', scale_units = 'xy', headaxislength=2.75, headlength = 5, headwidth=4.8, scale_type="absolute", quiver_scale = 1, plot_random = True, min_magnitude = None, min_mass = 0, plot_dots = False)
    #plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="w", lw=6, zorder=1000000)
    #plt.plot(pc_obj.projections[pc_obj.ixsort,0], pc_obj.projections[pc_obj.ixsort,1], c="k", lw=3, zorder=2000000)
    plt.gca().invert_xaxis()
    #plt.axis("off")
    plt.axis("equal")
    plt.savefig(sys.argv[1].split(".loom")[0] + "_velocity.png")

def calculate_gammas(vlm):
    #vlm.normalize("S", size=True, log=False)
    #vlm.normalize("U", size=True,  log=False)
    vlm.knn_imputation(k = 1)
    vlm.normalize_median()
    vlm.fit_gammas(limit_gamma=False, fit_offset=False)
    plt.figure(None, (17,2.8), dpi=80)
    gs = plt.GridSpec(1,6)
    for i, gn in enumerate(vlm.ra["Gene"][1:5]):
        ax = plt.subplot(gs[i])
        try:
            ix=np.where(vlm.ra["Gene"] == gn)[0][0]
        except:
            continue
        vcy.scatter_viz(vlm.S[ix,:], vlm.U[ix,:], s=5, alpha=0.4, rasterized=True)
        plt.title(gn)
        xnew = np.linspace(0,vlm.S[ix,:].max())
        plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
        plt.ylim(0, np.max(vlm.U[ix,:])*1.02)
        plt.xlim(0, np.max(vlm.S[ix,:])*1.02)
        #minimal_yticks(0, np.max(vlm.U[ix,:])*1.02)
        #minimal_xticks(0, np.max(vlm.S[ix,:])*1.02)
    plt.savefig(sys.argv[1].split(".loom")[0] + "_gn" + "_gamma.png")

def pc_plot(vlm):
    vlm.perform_PCA()
    plt.figure(None, (17,3.5))
    vcy.scatter_viz(vlm.pcs[:,0], vlm.pcs[:,1], s=10)
    plt.xlabel("PC1"); plt.ylabel("PC2")
    plt.savefig(sys.argv[1].split(".loom")[0] + "_pca.png")

if __name__ == "__main__":
    main()
