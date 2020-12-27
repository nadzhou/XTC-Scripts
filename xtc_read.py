#!/home/bin/python

import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.analysis import psa

import numpy as np
from MDAnalysis.analysis import distances
from numpy.linalg import norm

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.analysis.distances import distance_array

import matplotlib.pyplot as plt
import seaborn as sns

import parse_arguments as pa

def main(): 
    parser = pa.parse_arguments()

    print("loading universe...")
    new_universe = mda.Universe(parser.input_tpr, parser.input_xtc)

    if new_universe: 
        analysis_list()
        operation_input = int(input("Select which simulation analysis: "))
        region_selection_check()
        input_resid = input("Input selection parameters: ")

        analysis_routes(operation_input, new_universe, input_resid)


def region_selection_check(): 
    print("\nInitiating atom selection...\n")
    print("If you want to select a range of atoms, use resid\n\t\te.g., resid 1:240\n")
    print("For domain, just write the domain name, \n\t\te.g., SH3\n")
    print("For residue, write resid,\n\t\te.g., resname TER\n")


def analysis_list(): 
    print("\n1. PCA\t2. Hydrogen Bond lifetimes")
    print("3. Angles\t4. RMSF\n")


def analysis_routes(analysis_choice, universe, atom_select): 
    if analysis_choice == 1: 
        calc_pca(universe, atom_select)

    elif analysis_choice == 4: 
        calc_rmsf(universe, atom_select)

    
def calc_rmsf(universe: mda.Universe, atom_select: mda.Universe):
    calphas = universe.select_atoms("name CA")
    rmsf = mda.analysis.rms.RMSF(atom_select, verbose=True).run()
    rmsf2 = mda.analysis.rms.RMSF(calphas, verbose=True).run()

    print(rmsf.rmsf)
    sns.lineplot(atom_select.resnums, rmsf.rmsf, label="SER residues")
    sns.lineplot(calphas.resnums, rmsf2.rmsf, label="C Alpha")

    plt.xlabel("Residue number")
    plt.ylabel("RMSF values")
    plt.title("RMSF")
    plt.savefig("RMSF.png")
    plt.show()


def theta_NMP(u):
    """Calculate the NMP-CORE angle for E. coli AdK in degrees"""
    C = u.select_atoms("resid 115:125 and backbone").center_of_geometry()
    B = u.select_atoms("resid 90:100 and backbone").center_of_geometry()
    A = u.select_atoms("resid 35:55 and backbone").center_of_geometry()
    BA = A - B
    BC = C - B
    theta = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
    return np.rad2deg(theta)



def calc_distances_bn(ser_resids, tyr_resids) -> np.ndarray: 
    """Given two residues, calculate the distance matrix 
    """
    
    print("file reading done.")

    print("selecting residues...")

    ser_tyr_dists = distances.distance_array(ser_resids.positions, tyr_resids.positions)
    # tYr_self_distances = distances.self_distance_array(tyr_resids.positions)
    # # self_distances = np.where(self_distances >= 10, np.max(self_distances), self_distances)

    print(ser_tyr_dists.shape)
    print(ser_tyr_dists)
    sns.heatmap(ser_tyr_dists, alpha=0.6)
    plt.savefig("heatmap.png")
    plt.title("Heat map of two residues")
    plt.xlabel("Serine distances")
    plt.ylabel("Tyrosine distances")
    plt.show()



def ser_tyr_residslc_hbonds(universe): 
    """EXPERIMENTAL: 
            calculate hydrogens"""

    nad_bonds = HBA(universe=universe)
    nad_bonds.run(verbose=True)



def calc_pca(universe, atom_select):
    import MDAnalysis.analysis.pca as pca
    import pandas as pd
    import seaborn as sns
    
    pc = pca.PCA(universe, select=atom_select)
    pc.run(verbose=True)

    # print(tyr_pca.p_components.shape)
    # plt.plot(tyr_pca.mean_atoms.positions)
    # plt.xlabel('Principal component')
    # plt.ylabel('Cumulative variance')
    # plt.show()

    backbone = universe.select_atoms(atom_select)
    transformed = pc.transform(backbone, n_components=5)
    pc_df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(5)])
    pc_df['Time (ps)'] = pc_df.index * universe.trajectory.dt

    print(pc_df.head())
    sns.pairplot(pc_df, hue='Time (ps)')
    plt.savefig("pairgrid_normal2.png")
    # sns.map(plt.scatter, marker='.')
    plt.show()

if __name__ == '__main__': 
    main()

