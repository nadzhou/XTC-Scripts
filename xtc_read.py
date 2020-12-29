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

    analysis_list()
    operation_input = int(input("Select which simulation analysis: "))
    region_selection_check()
    input_resid = input("Input selection parameters: ")

    analysis_routes(operation_input, parser, input_resid)
    plt.show()


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

    elif analysis_choice == 2: 
        calc_hbonds(universe)

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



def calc_hbonds(universe): 
    """EXPERIMENTAL: 
            calculate hydrogens"""

    hbonds = HBA(universe=universe)
    hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
    hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
    hbonds.run(verbose=True)

    tau_timeseries, timeseries = hbonds.lifetime()
    print(tau_timeseries)
    print(timeseries)

    sns.lineplot(tau_timeseries, timeseries)
    plt.show()



def calc_pca(parser, atom_select):
    import MDAnalysis.analysis.pca as pca
    import pandas as pd
    import seaborn as sns

    universe_1 = mda.Universe(parser.input_tpr, parser.input_xtc)
    universe_2 = mda.Universe(parser.input_tpr2, parser.input_xtc2)
    universe_3 = mda.Universe(parser.input_tpr3, parser.input_xtc3)
    universe_4 = mda.Universe(parser.input_tpr4, parser.input_xtc4)
    universe_5 = mda.Universe(parser.input_tpr5, parser.input_xtc5)
    universe_6 = mda.Universe(parser.input_tpr6, parser.input_xtc6)

    pc = pca.PCA(universe_1, select=atom_select)
    pc2 = pca.PCA(universe_2, select=atom_select)
    pc3 = pca.PCA(universe_3, select=atom_select)
    pc4 = pca.PCA(universe_4, select=atom_select)
    pc5 = pca.PCA(universe_5, select=atom_select)
    pc6 = pca.PCA(universe_6, select=atom_select)
    
    pc.run(verbose=True)
    pc2.run(verbose=True)
    pc3.run(verbose=True)
    pc4.run(verbose=True)
    pc5.run(verbose=True)
    pc6.run(verbose=True)

    backbone = universe_1.select_atoms(atom_select)
    backbone2 = universe_2.select_atoms(atom_select)
    backbone3 = universe_3.select_atoms(atom_select)
    backbone4 = universe_4.select_atoms(atom_select)
    backbone5 = universe_5.select_atoms(atom_select)
    backbone6 = universe_6.select_atoms(atom_select)
   
    transformed = pc.transform(backbone, n_components=5)
    transformed2 = pc.transform(backbone2, n_components=5)
    transformed3 = pc.transform(backbone3, n_components=5)
    transformed4 = pc.transform(backbone4, n_components=5)
    transformed5 = pc.transform(backbone5, n_components=5)
    transformed6 = pc.transform(backbone6, n_components=5)

    pc_df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(5)])
    pc_df2 = pd.DataFrame(transformed2,
                  columns=['PC{}'.format(i+1) for i in range(5)])
    pc_df3 = pd.DataFrame(transformed3,
                  columns=['PC{}'.format(i+1) for i in range(5)])

    pc_df4 = pd.DataFrame(transformed4,
                  columns=['PC{}'.format(i+1) for i in range(5)])

    pc_df5 = pd.DataFrame(transformed5,
                  columns=['PC{}'.format(i+1) for i in range(5)])
    pc_df6 = pd.DataFrame(transformed6,
                  columns=['PC{}'.format(i+1) for i in range(5)])      

    pc_df['Time (ps)'] = pc_df.index * universe_1.trajectory.dt
    pc_df2['Time (ps)'] = pc_df2.index * universe_2.trajectory.dt
    pc_df3['Time (ps)'] = pc_df3.index * universe_3.trajectory.dt
    pc_df4['Time (ps)'] = pc_df4.index * universe_4.trajectory.dt
    pc_df5['Time (ps)'] = pc_df5.index * universe_5.trajectory.dt
    pc_df6['Time (ps)'] = pc_df6.index * universe_6.trajectory.dt

    pc.df.toexcel("1.xlsx")
    pc.df2.toexcel("2.xlsx")
    pc.df3.toexcel("3.xlsx")
    pc.df4.toexcel("4.xlsx")
    pc.df5.toexcel("5.xlsx")
    pc.df6.toexcel("6.xlsx")

    print(pc_df.head())
    sns.pairplot(pc_df, color="grey")
    sns.pairplot(pc_df2, color="red")
    sns.pairplot(pc_df3, color="orange")
    sns.pairplot(pc_df4, color="green")
    sns.pairplot(pc_df5, color="blue")
    sns.pairplot(pc_df6, color="purple")

    plt.savefig("pairgrid_normal2.png")
    # sns.map(plt.scatter, marker='.')
    plt.show()

if __name__ == '__main__': 
    main()

