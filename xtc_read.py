#!usr/bin/env/python 3
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.analysis import psa
import matplotlib.pyplot as plt

from MDAnalysis.analysis import distances
import numpy as np

# Function for ser_tyr_residslculating H-bonds. 
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


import pandas as pd
import seaborn as sns

from MDAnalysis.analysis.distances import distance_array

def main(): 

    # Input protein into MDAnalysis
    print("loading universe")
    nad_pt = mda.Universe("/home/nadzhou/Mehreen/md.tpr", "/home/nadzhou/Mehreen/140-147ns.xtc")

    # # select atoms. In this ser_tyr_residsse only TYR
    # print(nad_pt.select_atoms("resname TYR"))
    # print("ser_tyr_residslculating hydrogens...")
    
    # ser_tyr_residslc_hbonds(nad_pt)

    # Select SER and TYR
    ser_resids = nad_pt.select_atoms('resname SER')
    tyr_resids = nad_pt.select_atoms('resname TYR')


    print("file reading done.")

    print("selecting residues...")

    ser_tyr_dists = distances.self_distance_array(ser_resids.positions)
    # tYr_self_distances = distances.self_distance_array(tyr_resids.positions)
    # # self_distances = np.where(self_distances >= 10, np.max(self_distances), self_distances)


    print(ser_tyr_dists.reshape(-1, 2))
    sns.histplot(ser_tyr_dists.reshape(-1, 2), alpha=0.6)
    plt.show()

    # calc_distance_plot(self_distances, ser_tyr_resids)

    # print(ser_resids.positions)
    # print(tYr_self_distances.shape)



    # calc_pca(nad_pt)

    # sns.lineplot(ser_self_distances, tYr_self_distances)
    # plt.show()



def calc_distance(atomgroup_select_1: "AtomGroup", atom_group_select_2: "AtomGroup"): 
    """Calculate distance between two given atom group"""
    return distance_array(atomgroup_select_1.positions, 
                                   atom_group_select_2.positions
                                   )


def calc_distance_plot(self_distances, residue_select): 
    """ Calculate self distance

    Args: 
        self_distances [np..ndarray]: distances of atoms within the structure
        residue_select [MDAnalysis.Universe]: Select the atoms or residues of interest
    
    """
    n_residue_select = len(residue_select)

    sq_dist_arr = np.zeros((n_residue_select, n_residue_select))
    triu = np.triu_indices_from(sq_dist_arr, k=1)

    sq_dist_arr[triu] = self_distances
    sq_dist_arr.T[triu] = self_distances

    fig, ax = plt.subplots()
    im = ax.pcolor(ser_tyr_resids.resids, residue_select.resids, sq_dist_arr, cmap=plt.cm.get_cmap('Blues_r'))

    # plt.pcolor gives a rectangular grid by default
    # so we need to make our heatmap square
    ax.set_aspect('equal')

    # add figure labels and titles
    plt.ylabel('Residue IDs')
    plt.xlabel('Residue IDs')
    plt.title('Atomic distances between SER or TYR atoms')

    # colorbar
    cbar = fig.colorbar(im)
    plt.savefig("distance_heat_map.png")
    plt.show()


def ser_tyr_residslc_hbonds(universe): 
    """EXPERIMENTAL: 
            calculate hydrogens"""

    nad_bonds = HBA(universe=universe)
    print(nad_bonds.run())
    tau_timeseries, timeseries = nad_bonds.lifetime()

    print(tau_timeseries)



def calc_pca(universe):
    import MDAnalysis.analysis.pca as pca
    
    tyr_pca = pca.PCA(universe, select='backbone')
    tyr_pca.run()

    print(tyr_pca)


if __name__ == '__main__': 
    main()




# EXPERIMENTAL: 
"""
def ser_tyr_residslc_box(navid_pt, haf_pt): 

    ser_tyr_resids1 = nad_pt.select_atoms('name ser_tyr_resids')
    ser_tyr_resids2 = haf_pt.select_atoms('name ser_tyr_resids')

    resids1_box, resids2_box, dist_box = distances.dist(ser_tyr_resids1, ser_tyr_resids2,
                                                    box=[10, 10, 10, 90, 90, 90])



    plt.plot(resids1_box, dist_box)
    plt.ylabel('ser_tyr_resids distance (Angstrom)')
    plt.axvspan(122, 159, zorder=0, alpha=0.2, color='orange', label='LID')
    plt.axvspan(30, 59, zorder=0, alpha=0.2, color='green', label='NMP')
    plt.legend()
    plt.show()



def ser_tyr_residslculate_distance(nad_pt): 
    CORE_sel = 'name ser_tyr_resids and (resid 1:15)'


    labels = ['DCD', 'DCD2', 'XTC', 'NAMD', 'mixed']

    ps = psa.PSAnalysis([nad_pt, nad_pt], 
                        labels=labels,
                        reference=nad_pt, 
                        path_select='name ca')


    ps.generate_paths(align=True, save=False, weights='mass')
    ps.run(metric='hausdorff')

    print(ps.D)
    ps.plot(linkage='ward')
    plt.show()


"""