import MDAnalysis as mda
import pandas as pd
from pathlib import Path
import MDAnalysis.analysis.rms
from MDAnalysis.analysis import align
import numpy as np
from libmdanalysis import (
    calculate_quantities,
    chain_pdb,
    reres_pdb,
    mkensemble_pdb,
)
import os
import shutil
import argparse


from Bio import PDB
import numpy as np

parser = PDB.PDBParser(QUIET=True)
downloader = PDB.PDBList()
# Atomic radii for various atom types. 
# You can comment out the ones you don't care about or add new ones
atom_radii = {
    "C": 1.70, 
    "N": 1.55, 
    "O": 1.52,
    "S": 1.80,
    "F": 1.47, 
    "P": 1.80, 
    "CL": 1.75, 
    "MG": 1.73,
}

def count_backbone_clashes_h3(structure, h3_residues, light_chain_start, clash_cutoff=0.63):
    """
    count clashes between the backbone atoms of the H3 loop and the light chain

    inspired by the implementation of count_clashes at https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/

    Args:
    - structure: Bio.PDB.Structure
    - h3_residues: list of residues of the H3 loop
    - light_chain_start: int, residue number where the light chain starts
    - clash_cutoff: float, cutoff for clashes
    """
    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}
    # Extract atoms for which we have a radii
    backbone_atoms_h3 = [x for x in structure.get_atoms() if x.element in atom_radii and x.parent.id[1] in h3_residues and x.name in ["N", "CA", "C", "O"]]
    light_chain_atoms = [x for x in structure.get_atoms() if x.element in atom_radii and x.parent.id[1] > light_chain_start and x.name in ["N", "CA", "C", "O"]]
    coords = np.array([a.coord for a in light_chain_atoms], dtype="d")
    # Build a KDTree
    kdt = PDB.kdtrees.KDTree(coords)
    # Initialize a list to hold clashes
    clashes = []
    # Iterate through all atoms
    for atom_1 in backbone_atoms_h3:
        # Find atoms that could be clashing
        kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))
        # Get index and distance of potential clashes
        potential_clash = [(a.index, a.radius) for a in kdt_search]
        for ix, atom_distance in potential_clash:
            atom_2 = light_chain_atoms[ix]
            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))
    return len(clashes)
    

if __name__ == "__main__":

    p = argparse.ArgumentParser(description=
    """
    this code creates an ensemble of full antibody models starting from clustered data
    the input is the path to the dataframe containing the clustered data
    """)
    p.add_argument("clustering_data", type=str, help="clustering data")
    p.add_argument("heavy_chain_reference", type=str, help="Reference heavy chain")
    p.add_argument("light_chain_reference", type=str, help="Reference light chain")
    p.add_argument("--h1_loop_residues", type=str, help="Residues of the H1 loop")
    p.add_argument("--h2_loop_residues", type=str, help="Residues of the H2 loop")
    p.add_argument("--h3_loop_residues", type=str, help="Residues of the H3 loop")
    p.add_argument("--pdb_code", type=str, help="pdb code")
    p.add_argument("--topn", type=int, help="topn")
    p.add_argument("--cl_num", type=int, help="cl_num")
    args = p.parse_args()

    clustering_data = Path(args.clustering_data)
    # read the dataframe
    df = pd.read_csv(clustering_data, sep="\t")
    cl_method = "complete"
    filtered_df = df[(df["topn"] == args.topn) & (df["cl_num"] == args.cl_num) & (df["cl_method"] == cl_method)]
    print(f"filtered_df: {filtered_df}")

    print(f"reference heavy chain is {args.heavy_chain_reference}")
    ref_uni = mda.Universe(args.heavy_chain_reference)
    n_residues_heavy = ref_uni.atoms.n_residues
    print(f"number of residues in the heavy chain: {n_residues_heavy}")
    # light chain file
    l_chain_uni = mda.Universe(args.light_chain_reference)

    # create the ensemble directory
    if Path("clustered_alphaflow_ensembles").exists() is False:
        os.mkdir("clustered_alphaflow_ensembles")
    ensemble_dir = Path("clustered_alphaflow_ensembles", f"{args.pdb_code}_ensemble_{args.topn}_{args.cl_num}")
    if ensemble_dir.exists() is False:
        os.mkdir(ensemble_dir)

    models = filtered_df["cluster_center"].unique()
    models_files = [Path(clustering_data.parent / f"{model}.pdb") for model in models]

    center_plddt = []
    center_h3rmsd_af = []
    center_h3rmsd = []
    center_nclashes = []
    center_h3_nclashes = []

    full_models_files = []

    for model in models_files:
        model_uni = mda.Universe(model)
        # calculate the plddt
        plddt = True
        # calculate the quantities
        rmsd_obj = calculate_quantities(model_uni, ref_uni, args.h1_loop_residues, args.h2_loop_residues, args.h3_loop_residues,
                                        l1_loop_residues=None, l2_loop_residues=None, l3_loop_residues=None,
                                        h3_only=False, plddt=plddt)
        # print stuff
        h3_plddt = rmsd_obj.plddt_vector[3]
        # print(f"model: {model} rmsd_full: {rmsd_obj.rmsd_full:.2f} rmsd_framework: {rmsd_obj.rmsd_framework:.2f} rmsd_h1: {rmsd_obj.rmsd_h1:.2f} rmsd_h2: {rmsd_obj.rmsd_h2:.2f} rmsd_h3: {rmsd_obj.rmsd_h3:.2f}")
        # plot correlation
        center_plddt.append(h3_plddt)
        center_h3rmsd_af.append(rmsd_obj.rmsd_h3)
        # let's extract center_h3rmsd from filtered_df
        h3_rmsd = filtered_df[filtered_df["cluster_center"] == model.stem]["cluster_center_h3_rmsd"].values[0]
        center_h3rmsd.append(h3_rmsd)

        model_uni = mda.Universe(model)
        model_uni.atoms.translate(-rmsd_obj.mod_com)
        model_uni.atoms.rotate(rmsd_obj.R)
        model_uni.atoms.translate(rmsd_obj.ref_com)

        #  suppressing warnings
        import warnings
        warnings.simplefilter("ignore")
        rotated_filename = ensemble_dir / f"{model.stem}_heavy_rotated.pdb"
        with MDAnalysis.Writer(rotated_filename) as pdb:
            pdb.write(model_uni)

        # merge with the rotated chains
        merged_filename = Path(ensemble_dir / f"{model.stem}_merged.pdb")
        model_uni = mda.Merge(model_uni.atoms, l_chain_uni.atoms)
        model_uni.atoms.write(merged_filename)

        # chain and reres the merged file
        chain_filename = chain_pdb(merged_filename, "A")
        reres_filename = reres_pdb(chain_filename, 1)

        new_reres_filename = ensemble_dir / f"{model.stem}_merged_reres.pdb"
        shutil.move(reres_filename, new_reres_filename)
        # remove intermediate files. comment for debugging purposes
        os.unlink(chain_filename)
        os.unlink(merged_filename)
        os.unlink(rotated_filename)
        
        full_models_files.append(new_reres_filename)

        # counting clashes
        pdb_parsed_filename = parser.get_structure(new_reres_filename, new_reres_filename)
        clashes_backbone_h3 = count_backbone_clashes_h3(pdb_parsed_filename, rmsd_obj.h3_residues, n_residues_heavy, clash_cutoff=0.5)
        print(f"model: {model} h3_rmsd {h3_rmsd} clashes_backbone_h3: {clashes_backbone_h3}")
        center_h3_nclashes.append(clashes_backbone_h3)


    # save rmsd and clashes data
    data = {"plddt": center_plddt, "h3rmsd_af2": center_h3rmsd_af, "h3rmsd": center_h3rmsd, "h3_nclashes": center_h3_nclashes, "full_models_files": full_models_files}
    df = pd.DataFrame(data)
    df.to_csv(ensemble_dir / f"ensemble_{args.pdb_code}_top{args.topn}_cl{args.cl_num}.csv", index=False, sep="\t", float_format="%.2f")

    # create the ensemble
    mkensemble_pdb(full_models_files)
    shutil.move("ensemble.pdb", ensemble_dir / f"ensemble_{args.pdb_code}_top{args.topn}_cl{args.cl_num}.pdb")
    print(f"ensemble {args.topn} {args.cl_num} created")

    # calculate correlation coefficient
    correlation = np.corrcoef(center_h3rmsd, center_h3_nclashes)[0][1]
    print(f"correlation h3 rmsd and number of h3 clashes: {correlation:.3f}")
    # has the lowest h3_rmsd 0 clashes?
    print(f"lowest h3_rmsd: {min(center_h3rmsd)}")
    print(f"number of clashes for the lowest h3_rmsd: {center_h3_nclashes[center_h3rmsd.index(min(center_h3rmsd))]}")

    # create the ensemble using only the models with at most 1 h3_clashes
    # sort the dataframe by clashes
    df = df.sort_values(by="h3_nclashes")
    noclash_df = df[df["h3_nclashes"] <= 1]
    noclash_df.to_csv(ensemble_dir / f"ensemble_{args.pdb_code}_top{args.topn}_cl{args.cl_num}_noclash.csv", index=False, sep="\t", float_format="%.2f")
    print(f"selecting {noclash_df.shape[0]} models with at most one backbone clash in the H3 loop")

    mkensemble_pdb(noclash_df["full_models_files"])
    shutil.move("ensemble.pdb", ensemble_dir / f"ensemble_{args.pdb_code}_top{args.topn}_cl{args.cl_num}_noclash.pdb")

    # remove all full_models_files. comment for debugging purposes
    for f in full_models_files:
        os.unlink(f)