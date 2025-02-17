import MDAnalysis as mda
import pandas as pd
from pathlib import Path
import numpy as np
import os
from scipy.stats import entropy
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
import sys
from libmdanalysis import (
    calculate_quantities,
    get_atoms,
    get_clustering_dict,
    get_cluster_center,
    load_coords,
    write_coords,
    chain_pdb,
    reres_pdb,
    mkensemble_pdb,
)
import argparse
    

if __name__ == "__main__":

    p = argparse.ArgumentParser(description=
    """
    Analyse the diffusion of an ensemble of models.
    """)
    p.add_argument("ensemble_dir", type=str, help="Directory with the ensemble of models")
    p.add_argument("ref_pdb", type=str, help="Reference PDB")
    p.add_argument("--h1_loop_residues", type=str, help="Residues of the H1 loop (format 'start':'end', example '25:32')")
    p.add_argument("--h2_loop_residues", type=str, help="Residues of the H2 loop (format 'start':'end', example '50:57')")
    p.add_argument("--h3_loop_residues", type=str, help="Residues of the H3 loop (format 'start':'end', example '96:114')")
    p.add_argument("--l1_loop_residues", type=str, help="Residues of the L1 loop (format 'start':'end', example '25:32')", default=None, required=False)
    p.add_argument("--l2_loop_residues", type=str, help="Residues of the L2 loop (format 'start':'end', example '49:51')", default=None, required=False)
    p.add_argument("--l3_loop_residues", type=str, help="Residues of the L3 loop (format 'start':'end', example '92:100')", default=None, required=False)
    p.add_argument("--cluster", action="store_true", help="Cluster structures", default=False, required=False)
    args = p.parse_args()

    ensemble_dir = Path(args.ensemble_dir)
    plddt = True if "colabfold" in ensemble_dir.name else False
    ref_pdb = Path(args.ref_pdb)
    print(f"ref pdb {ref_pdb}")
    
    ref_uni = mda.Universe(ref_pdb)
    pdb_files = [f for f in ensemble_dir.glob("*.pdb")]
    print(f"found {len(pdb_files)} pdb_files in {ensemble_dir}")
    # REFERENCE DATA
    if args.l1_loop_residues is None:
        filename = "ensemble_rmsd_vs_ref_heavy.csv"
        plddt_filename = "plddt_vs_ref_heavy.csv"
    else:
        filename = "ensemble_rmsd_vs_ref_full.csv"
        plddt_filename = "plddt_vs_ref_full.csv"
    reference_path = Path(ensemble_dir, filename)
    if reference_path.exists():
        # remove the file
        reference_path.unlink()

    data = []
    plddt_data = []

    for pdb_file in pdb_files:
        pdb_code = pdb_file.stem
        model_uni = mda.Universe(pdb_file)
        rmsd_obj = calculate_quantities(model_uni, ref_uni,
                                        args.h1_loop_residues, args.h2_loop_residues, args.h3_loop_residues,
                                        args.l1_loop_residues, args.l2_loop_residues, args.l3_loop_residues,
                                        h3_only=False, plddt=plddt)
        if args.l1_loop_residues is None:
            data.append([pdb_code, rmsd_obj.rmsd_full, rmsd_obj.rmsd_framework,
                         rmsd_obj.rmsd_h1, rmsd_obj.rmsd_h2, rmsd_obj.rmsd_h3])
        else:
            data.append([pdb_code, rmsd_obj.rmsd_full, rmsd_obj.rmsd_framework,
                         rmsd_obj.rmsd_h1, rmsd_obj.rmsd_h2, rmsd_obj.rmsd_h3,
                         rmsd_obj.rmsd_l1, rmsd_obj.rmsd_l2, rmsd_obj.rmsd_l3])
        if plddt:
            plddt_data.append([pdb_code] + rmsd_obj.plddt_vector)
    # minimal values for each quantity
    if args.l1_loop_residues is None:
        quantities = ["rmsd_full", "rmsd_framework", "rmsd_h1", "rmsd_h2", "rmsd_h3"]
    else:
        quantities = ["rmsd_full", "rmsd_framework", "rmsd_h1", "rmsd_h2", "rmsd_h3", "rmsd_l1", "rmsd_l2", "rmsd_l3"]
    df = pd.DataFrame(data, columns=["model"] + quantities)
    df.to_csv(reference_path, index=False, sep="\t", float_format="%.3f")
    for quantity in quantities:
        min_val = df[quantity].min()
        min_model = df.loc[df[quantity] == min_val, "model"].values[0]
        print(f"min {quantity} {min_val:.2f} model {min_model}")
    if plddt:
        if args.l1_loop_residues is None:
            df_plddt = pd.DataFrame(plddt_data, columns=["model", "plddt_full", "plddt_h1", "plddt_h2", "plddt_h3"])
        else:
            df_plddt = pd.DataFrame(plddt_data, columns=["model", "plddt_full", "plddt_h1", "plddt_h2", "plddt_h3", "plddt_l1", "plddt_l2", "plddt_l3"])
        df_plddt.to_csv(Path(ensemble_dir, plddt_filename), index=False, sep="\t", float_format="%.3f")

    # # now clustering
    if args.cluster is False:
        print("no clustering option specified, exiting...")
        sys.exit(1)
    # the rmsd_exec is in the directory where the python script is contained
    script_dir = Path(__file__).parent
    rmsd_exec = script_dir / "fast-rmsdmatrix/src/fast-rmsdmatrix"
    N_vector = [10, 20, 50, 100, 200, 500, 1000]
    cl_vector = [5, 10, 15, 20, 30, 40]
    cl_method_vector = ["complete", "average", "single"]
    cluster_centers_data = []
    # loops
    h1_start, h1_end = args.h1_loop_residues.split(":")
    h2_start, h2_end = args.h2_loop_residues.split(":")
    h3_start, h3_end = args.h3_loop_residues.split(":")
    loop_residues = list(range(int(h1_start), int(h1_end)+1)) + list(range(int(h2_start), int(h2_end)+1)) + list(range(int(h3_start), int(h3_end)+1))
    print(f"loop residues {loop_residues}")
    h3_residues = list(range(int(h3_start), int(h3_end)+1))
    assert rmsd_obj.h3_residues == h3_residues, f"rmsd_obj.h3_residues {rmsd_obj.h3_residues} != h3_residues {h3_residues}"
    for topn in N_vector:
        print(f"processing top {topn}")
        if topn > len(pdb_files):
            print(f"topn {topn} is larger than the number of pdb files {len(pdb_files)}, skipping")
            continue
        topn_pdb_files = [Path(ensemble_dir, f"model_{n}.pdb") for n in range(topn)]
        ligands = [pdb_file.stem for pdb_file in topn_pdb_files]
        rec_traj_filename = "traj_rec.xyz"
        lig_traj_filename = "traj_lig.xyz"
        n_rec_atoms, n_lig_atoms = write_coords(rec_traj_filename,
                                                lig_traj_filename,
                                                topn_pdb_files,
                                                loop_residues,
                                                h3_residues)
        # prepare command to be executed
        npairs = topn*(topn-1)//2
        cmd = f"'{rmsd_exec}' {rec_traj_filename} 0 {npairs} 0 1 {topn} {n_rec_atoms} {lig_traj_filename} {n_lig_atoms}"
        # execute command with multiprocessing
        import subprocess
        import shlex
        p = subprocess.run(
                shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        # load the rmsd matrix
        rmsd_matrix = np.loadtxt("ilrmsd_0.matrix")
        # extract only the third column
        rmsd_matrix = rmsd_matrix[:,2]
        # iterate over cl_num.
        for cl_num in cl_vector:
            if cl_num > topn:
                # cl_num is larger than topn, skipping
                continue
            for cl_method in cl_method_vector:
                print(f"clustering witn ncl = {cl_num} and linkage {cl_method}")
                clustering = AgglomerativeClustering(n_clusters=cl_num,
                                                     metric="precomputed",
                                                     linkage=cl_method).fit(squareform(rmsd_matrix))
                clusters = clustering.labels_
                probs = np.bincount(clusters)/len(clusters)
                cl_entropy = entropy(probs)
                cl_dict = get_clustering_dict(clusters, ligands[:topn]) 
                assert len(np.unique(clusters)) == cl_num
                # now for each cluster we want to calculate the dispersion of the h3 rmsd vs the reference  
                df_refe = pd.read_csv(reference_path, sep="\t")
                for cl in cl_dict:
                    df_cl = df_refe[df_refe["model"].isin(cl_dict[cl])]
                    # get h3_rmsd
                    h3_rmsd = df_cl["rmsd_h3"]
                    mean = h3_rmsd.mean()
                    std = h3_rmsd.std()
                    max_val = h3_rmsd.max()
                    min_val = h3_rmsd.min()
                    npw = np.where(clusters == cl)[0]
                    # cluster centers
                    cluster_center_idx = get_cluster_center(npw, rmsd_matrix, topn)
                    cluster_center = ligands[cluster_center_idx]
                    cluster_center_h3_rmsd = df_refe.loc[df_refe["model"] == cluster_center][["rmsd_h3"]].values[0][0]
                    cluster_centers_data.append([topn, cl_num, cl_method, cluster_center,
                                                 cluster_center_h3_rmsd, mean, std,
                                                 max_val, min_val, cl_entropy])
        # remove matrix and trajectories
        os.remove("ilrmsd_0.matrix")
        os.remove(rec_traj_filename)
        os.remove(lig_traj_filename)
    df_cluster_centers = pd.DataFrame(cluster_centers_data, columns=["topn", "cl_num", "cl_method",
                                                                     "cluster_center", "cluster_center_h3_rmsd",
                                                                     "mean", "std", "max", "min", "clustering_entropy"])
    df_cluster_centers.to_csv(Path(ensemble_dir, "cluster_centers_h3_rmsd.csv"), index=False, sep="\t", float_format="%.3f")    
    