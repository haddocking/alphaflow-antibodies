import MDAnalysis.analysis.rms
from MDAnalysis.analysis import align
import numpy as np
import MDAnalysis as mda
from pathlib import Path
import os
from pdbtools.pdb_mkensemble import make_ensemble
from pdbtools.pdb_reres import renumber_residues
from pdbtools.pdb_chain import alter_chain


class RMSD_object:
    def __init__(self):
        self.rmsd_full = -1
        self.rmsd_framework = -1
        self.rmsd_h1 = -1
        self.rmsd_h2 = -1
        self.rmsd_h3 = -1
        self.rmsd_l1 = -1
        self.rmsd_l2 = -1
        self.rmsd_l3 = -1
        self.plddt_vector = None
        self.R = None
        self.mod_com = None
        self.ref_com = None


def calculate_quantities(model_uni, ref_uni, h1_loop_residues, h2_loop_residues, h3_loop_residues, l1_loop_residues=None, l2_loop_residues=None, l3_loop_residues=None, h3_only=False, plddt=False):
    loop_dict = {}
    
    h1_start, h1_end = h1_loop_residues.split(":")
    h2_start, h2_end = h2_loop_residues.split(":")
    h3_start, h3_end = h3_loop_residues.split(":")
    h1_start = int(h1_start)
    h1_end = int(h1_end)
    h2_start = int(h2_start)
    h2_end = int(h2_end)
    h3_start = int(h3_start)
    h3_end = int(h3_end)
    
    rmsd_obj = RMSD_object()

    # light chain loop
    if l1_loop_residues is not None:
        l1_start, l1_end = l1_loop_residues.split(":")
        l2_start, l2_end = l2_loop_residues.split(":")
        l3_start, l3_end = l3_loop_residues.split(":")
        l1_start = int(l1_start)
        l1_end = int(l1_end)
        l2_start = int(l2_start)
        l2_end = int(l2_end)
        l3_start = int(l3_start)
        l3_end = int(l3_end)

    # overall rmsd
    sel = "protein and backbone"
    # chains
    model_chains = model_uni.select_atoms("protein").residues.segments.segids
    ref_chains = ref_uni.select_atoms("protein").residues.segments.segids
    # we assume heavy chain is the first
    model_h_chain = model_chains[0]
    ref_h_chain = ref_chains[0]
    
    if l1_loop_residues is None:
        selch = f"{sel} and chainID {model_h_chain}"
        rmsd_full = MDAnalysis.analysis.rms.rmsd(model_uni.select_atoms(selch).positions, ref_uni.select_atoms(selch).positions, superposition=True)
    else:
        selch = f"{sel} and chainID {model_h_chain} or {sel} and chainID {model_chains[1]}"
        rmsd_full = MDAnalysis.analysis.rms.rmsd(model_uni.select_atoms(selch).positions, ref_uni.select_atoms(selch).positions, superposition=True)
    
    rmsd_obj.rmsd_full = rmsd_full

    model_hresidues = model_uni.select_atoms(f"protein and name CA and chainID {model_h_chain}").resids
    ref_hresidues = ref_uni.select_atoms(f"protein and name CA and chainID {ref_h_chain}").resids
    
    loop_residues = list(range(int(h1_start), int(h1_end)+1)) + list(range(int(h2_start), int(h2_end)+1)) + list(range(int(h3_start), int(h3_end)+1))
    mod_framework_hresids = [el for el in model_hresidues if el not in loop_residues]
    ref_framework_hresids = [el for el in ref_hresidues if el not in loop_residues]
    # extracting the framework atoms
    mod_framework_hatoms = model_uni.select_atoms(f"protein and backbone and chainID {model_h_chain} and resid {' '.join([str(el) for el in mod_framework_hresids])}")
    ref_framework_hatoms = ref_uni.select_atoms(f"protein and backbone and chainID {ref_h_chain} and resid {' '.join([str(el) for el in ref_framework_hresids])}")
    
    if l1_loop_residues is not None:
        # light chain exists
        model_l_chain = model_chains[1]
        ref_l_chain = ref_chains[1]

        model_lresidues = model_uni.select_atoms(f"protein and name CA and chainID {model_l_chain}").resids
        ref_lresidues = ref_uni.select_atoms(f"protein and name CA and chainID {ref_l_chain}").resids
        l_loop_residues = list(range(int(l1_start), int(l1_end)+1)) + list(range(int(l2_start), int(l2_end)+1)) + list(range(int(l3_start), int(l3_end)+1))
        mod_framework_lresids = [el for el in model_lresidues if el not in l_loop_residues]
        ref_framework_lresids = [el for el in ref_lresidues if el not in l_loop_residues]
        mod_framework_latoms = model_uni.select_atoms(f"protein and backbone and chainID {model_l_chain} and resid {' '.join([str(el) for el in mod_framework_lresids])}")
        ref_framework_latoms = ref_uni.select_atoms(f"protein and backbone and chainID {ref_l_chain} and resid {' '.join([str(el) for el in ref_framework_lresids])}")
        # extend the framework atoms
        mod_framework_atoms = mda.Merge(mod_framework_hatoms, mod_framework_latoms).atoms
        ref_framework_atoms = mda.Merge(ref_framework_hatoms, ref_framework_latoms).atoms
    else:
        # they are the same
        mod_framework_atoms = mod_framework_hatoms
        ref_framework_atoms = ref_framework_hatoms
    
    mod_com = model_uni.atoms.center_of_mass()
    ref_com = ref_uni.atoms.center_of_mass()
    # print(f"mod fr shape {mod_framework_atoms.positions.shape} ref fr shape {ref_framework_atoms.positions.shape}")
    R, rmsd_framework = align.rotation_matrix(mod_framework_atoms.positions - mod_framework_atoms.center_of_mass(),
                                              ref_framework_atoms.positions - ref_framework_atoms.center_of_mass())
    rmsd_obj.R = R
    rmsd_obj.mod_com = mod_framework_atoms.center_of_mass()
    rmsd_obj.ref_com = ref_framework_atoms.center_of_mass()
    rmsd_obj.rmsd_framework = rmsd_framework
    if not h3_only:
        #H1
        mod_h1_sel = f"not name H* and backbone and chainID {model_h_chain} and ("
        h1_resids = [f" resid {n} " for n in range(h1_start, h1_end+1)]
        mod_h1_sel += " or ".join(h1_resids)
        mod_h1_sel += " )"
        mod_h1_atoms = model_uni.select_atoms(mod_h1_sel)
        h1_seq = mod_h1_atoms.residues.sequence(format="string")
        loop_dict["H1"] = h1_seq
        ref_h1_sel = f"not name H* and backbone and chainID {ref_h_chain} and ( "
        ref_h1_sel += " or ".join(h1_resids)
        ref_h1_sel += " )"
        ref_h1_atoms = ref_uni.select_atoms(ref_h1_sel)
        # h1 rmsd
        mod_h1_com_coords = mod_h1_atoms.select_atoms('not name H* and backbone').center_of_mass()
        ref_h1_com_coords = ref_h1_atoms.select_atoms('not name H* and backbone').center_of_mass()
        mod_h1_atoms.atoms.translate(-rmsd_obj.mod_com)
        mod_h1_atoms.atoms.rotate(R)
        mod_h1_atoms.atoms.translate(rmsd_obj.ref_com)
        rmsd_h1 = mda.analysis.rms.rmsd(mod_h1_atoms.positions, ref_h1_atoms.positions, superposition=False)
        rmsd_obj.rmsd_h1 = rmsd_h1
        ## H2
        mod_h2_sel = f"not name H* and backbone and chainID {model_h_chain} and ( "
        h2_resids = [f" resid {n} " for n in range(h2_start, h2_end+1)]
        mod_h2_sel += " or ".join(h2_resids)
        mod_h2_sel += " )"
        mod_h2_atoms = model_uni.select_atoms(mod_h2_sel)
        h2_seq = mod_h2_atoms.residues.sequence(format="string")
        loop_dict["H2"] = h2_seq
        ref_h2_sel = f"not name H* and backbone and chainID {ref_h_chain} and ( "
        ref_h2_sel += " or ".join(h2_resids)
        ref_h2_sel += " )"
        ref_h2_atoms = ref_uni.select_atoms(ref_h2_sel)
        # h2 rmsd
        mod_h2_com_coords = mod_h2_atoms.select_atoms('not name H* and backbone').center_of_mass()
        ref_h2_com_coords = ref_h2_atoms.select_atoms('not name H* and backbone').center_of_mass()
        mod_h2_atoms.atoms.translate(-rmsd_obj.mod_com)
        mod_h2_atoms.atoms.rotate(R)
        mod_h2_atoms.atoms.translate(rmsd_obj.ref_com)
        rmsd_h2 = mda.analysis.rms.rmsd(mod_h2_atoms.positions, ref_h2_atoms.positions, superposition=False)
        rmsd_obj.rmsd_h2 = rmsd_h2

        if l1_loop_residues is not None:
            # L1
            mod_l1_sel = f"not name H* and backbone and chainID {model_l_chain} and ( "
            l1_resids = [f" resid {n} " for n in range(l1_start, l1_end+1)]
            mod_l1_sel += " or ".join(l1_resids)
            mod_l1_sel += " )"
            mod_l1_atoms = model_uni.select_atoms(mod_l1_sel)
            l1_seq = mod_l1_atoms.residues.sequence(format="string")
            loop_dict["L1"] = l1_seq
            ref_l1_sel = f"not name H* and backbone and chainID {ref_l_chain} and ( "
            ref_l1_sel += " or ".join(l1_resids)
            ref_l1_sel += " )"
            ref_l1_atoms = ref_uni.select_atoms(ref_l1_sel)
            # l1 rmsd
            mod_l1_com_coords = mod_l1_atoms.select_atoms('not name H* and backbone').center_of_mass()
            ref_l1_com_coords = ref_l1_atoms.select_atoms('not name H* and backbone').center_of_mass()
            mod_l1_atoms.atoms.translate(-rmsd_obj.mod_com)
            mod_l1_atoms.atoms.rotate(R)
            mod_l1_atoms.atoms.translate(rmsd_obj.ref_com)
            rmsd_l1 = mda.analysis.rms.rmsd(mod_l1_atoms.positions, ref_l1_atoms.positions, superposition=False)

            # L2
            mod_l2_sel = f"not name H* and backbone and chainID {model_l_chain} and ( "
            l2_resids = [f" resid {n} " for n in range(l2_start, l2_end+1)]
            mod_l2_sel += " or ".join(l2_resids)
            mod_l2_sel += " )"
            mod_l2_atoms = model_uni.select_atoms(mod_l2_sel)
            l2_seq = mod_l2_atoms.residues.sequence(format="string")
            loop_dict["L2"] = l2_seq
            ref_l2_sel = f"not name H* and backbone and chainID {ref_l_chain} and ( "
            ref_l2_sel += " or ".join(l2_resids)
            ref_l2_sel += " )"
            ref_l2_atoms = ref_uni.select_atoms(ref_l2_sel)

            # l2 rmsd
            mod_l2_com_coords = mod_l2_atoms.select_atoms('not name H* and backbone').center_of_mass()
            ref_l2_com_coords = ref_l2_atoms.select_atoms('not name H* and backbone').center_of_mass()
            mod_l2_atoms.atoms.translate(-rmsd_obj.mod_com)
            mod_l2_atoms.atoms.rotate(R)
            mod_l2_atoms.atoms.translate(rmsd_obj.ref_com)
            rmsd_l2 = mda.analysis.rms.rmsd(mod_l2_atoms.positions, ref_l2_atoms.positions, superposition=False)

            # L3
            mod_l3_sel = f"not name H* and backbone and chainID {model_l_chain} and ( "
            l3_resids = [f" resid {n} " for n in range(l3_start, l3_end+1)]
            mod_l3_sel += " or ".join(l3_resids)
            mod_l3_sel += " )"
            mod_l3_atoms = model_uni.select_atoms(mod_l3_sel)
            l3_seq = mod_l3_atoms.residues.sequence(format="string")
            loop_dict["L3"] = l3_seq
            ref_l3_sel = f"not name H* and backbone and chainID {ref_l_chain} and ( "
            ref_l3_sel += " or ".join(l3_resids)
            ref_l3_sel += " )"
            ref_l3_atoms = ref_uni.select_atoms(ref_l3_sel)

            # l3 rmsd
            mod_l3_com_coords = mod_l3_atoms.select_atoms('not name H* and backbone').center_of_mass()
            ref_l3_com_coords = ref_l3_atoms.select_atoms('not name H* and backbone').center_of_mass()
            mod_l3_atoms.atoms.translate(-rmsd_obj.mod_com)
            mod_l3_atoms.atoms.rotate(R)
            mod_l3_atoms.atoms.translate(rmsd_obj.ref_com)
            rmsd_l3 = mda.analysis.rms.rmsd(mod_l3_atoms.positions, ref_l3_atoms.positions, superposition=False)
            # attributes assignment
            rmsd_obj.rmsd_l1 = rmsd_l1
            rmsd_obj.rmsd_l2 = rmsd_l2
            rmsd_obj.rmsd_l3 = rmsd_l3

    # H3
    mod_h3_sel = f"not name H* and backbone and chainID {model_h_chain} and ( "
    h3_resids = [f" resid {n} " for n in range(h3_start, h3_end+1)]
    h3_residues = [n for n in range(h3_start, h3_end+1)]
    rmsd_obj.h3_residues = h3_residues
    mod_h3_sel += " or ".join(h3_resids)
    mod_h3_sel += " )"
    mod_h3_atoms = model_uni.select_atoms(mod_h3_sel)
    # print sequence
    h3_seq = mod_h3_atoms.residues.sequence(format="string")
    loop_dict["H3"] = h3_seq
    ref_h3_sel = f"not name H* and backbone and chainID {ref_h_chain} and ( "
    ref_h3_sel += " or ".join(h3_resids)
    ref_h3_sel += " )"
    ref_h3_atoms = ref_uni.select_atoms(ref_h3_sel)
    
    # h3 rmsd
    mod_h3_com_coords = mod_h3_atoms.select_atoms('not name H* and backbone').center_of_mass()
    ref_h3_com_coords = ref_h3_atoms.select_atoms('not name H* and backbone').center_of_mass()
    mod_h3_atoms.atoms.translate(-rmsd_obj.mod_com)
    mod_h3_atoms.atoms.rotate(R)
    mod_h3_atoms.atoms.translate(rmsd_obj.ref_com)
    rmsd_h3 = mda.analysis.rms.rmsd(mod_h3_atoms.positions, ref_h3_atoms.positions, superposition=False)
    rmsd_obj.rmsd_h3 = rmsd_h3

    plddt_vector = []
    if plddt:
        # gets the unique residues
        plddt_full = np.mean(model_uni.select_atoms('name CA').tempfactors)
        plddt_h1 = np.mean(mod_h1_atoms.select_atoms('name CA').tempfactors)
        plddt_h2 = np.mean(mod_h2_atoms.select_atoms('name CA').tempfactors)
        plddt_h3 = np.mean(mod_h3_atoms.select_atoms('name CA').tempfactors)
        plddt_vector.append(plddt_full)
        plddt_vector.append(plddt_h1)
        plddt_vector.append(plddt_h2)
        plddt_vector.append(plddt_h3)
        
        if l1_loop_residues is not None:
            plddt_l1 = np.mean(mod_l1_atoms.select_atoms('name CA').tempfactors)
            plddt_l2 = np.mean(mod_l2_atoms.select_atoms('name CA').tempfactors)
            plddt_l3 = np.mean(mod_l3_atoms.select_atoms('name CA').tempfactors)
            plddt_vector.append(plddt_l1)
            plddt_vector.append(plddt_l2)
            plddt_vector.append(plddt_l3)
    rmsd_obj.plddt_vector = plddt_vector

    # # last check on loop data
    # import pandas as pd
    # print(f"loop_dict {loop_dict}")
    # df_loop = pd.read_csv("../data/loop_data.csv", sep="\t")
    # # h3_seq must match the sequence in the loop_data.csv
    # for key in loop_dict:
    #     loop_seq = loop_dict[key]
    #     assert (loop_seq == df_loop[key]).any(), f"sequence {loop_seq} not found in loop_data.csv"

    return rmsd_obj

PROT_RES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]

# Backbone
PROT_ATOMS = ["C", "N", "CA", "O"]
# Side chains
PROT_SIDE_CHAINS_DICT = {
    "ALA": ["C", "N", "CA", "O", "CB"],
    "ARG": ["C", "N", "CA", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["C", "N", "CA", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["C", "N", "CA", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["C", "N", "CA", "O", "CB", "SG"],
    "GLN": ["C", "N", "CA", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["C", "N", "CA", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["C", "N", "CA", "O"],
    "HIS": ["C", "N", "CA", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["C", "N", "CA", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["C", "N", "CA", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["C", "N", "CA", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["C", "N", "CA", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["C", "N", "CA", "O", "CB", "CG", "CD"],
    "SER": ["C", "N", "CA", "O", "CB", "OG"],
    "THR": ["C", "N", "CA", "O", "CB", "OG1", "CG2"],
    "TRP": [
        "C",
        "N",
        "CA",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],  # noqa: E501
    "TYR": [
        "C",
        "N",
        "CA",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "CE1",
        "CE2",
        "CZ",
        "OH",
    ],  # noqa: E501
    "VAL": ["C", "N", "CA", "O", "CB", "CG1", "CG2"],
}

def get_atoms(pdb: Path, full: bool = False):
    """
    Identify what is the molecule type of each PDB.

    Parameters
    ----------
    pdb : PosixPath or :py:class:`haddock.libs.libontology.PDBFile`
        PDB file to have its atoms identified

    Returns
    -------
    atom_dic : dict
        dictionary of atoms
    """
    atom_dic: AtomsDict = {}
    atom_dic.update((r, PROT_ATOMS) for r in PROT_RES)
    if full:
        atom_dic.update(PROT_SIDE_CHAINS_DICT)

    with open(pdb) as fh:
        for line in fh.readlines():
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                if (
                    resname not in PROT_RES
                ):
                    if element != "H":
                        if resname not in atom_dic:
                            atom_dic[resname] = []
                        if atom_name not in atom_dic[resname]:
                            atom_dic[resname].append(atom_name)
    return atom_dic

def load_coords(
        pdb_f,
        atoms,
        filter_resdic=None,
        numbering_dic=None,
        model2ref_chain_dict=None,
        add_resname=None,
        ):
    """
    Load coordinates from PDB.

    Parameters
    ----------
    pdb_f : PDBFile

    atoms : dict
        dictionary of atoms

    filter_resdic : dict
        dictionary of residues to be loaded (one list per chain)

    numbering_dic : dict
        dict of numbering dictionaries (one dictionary per chain)

    add_resname : bool
        use the residue name in the identifier

    Returns
    -------
    coord_dic : dict
        dictionary of coordinates (one per chain)

    chain_ranges: dict
        dictionary of chain ranges
    """
    coord_dic: CoordsDict = {}
    chain_dic: ResDict = {}
    idx = 0
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                if model2ref_chain_dict:
                    chain = model2ref_chain_dict[line[21]]
                else:
                    chain = line[21]
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords = np.asarray([x, y, z])
                if numbering_dic and model2ref_chain_dict:
                    try:
                        resnum = numbering_dic[chain][resnum]
                    except KeyError:
                        continue
                if add_resname is True:
                    identifier = (chain, resnum, atom_name, resname)
                else:
                    identifier = (chain, resnum, atom_name)
                if atom_name not in atoms[resname]:
                    continue
                if chain not in chain_dic:
                    chain_dic[chain] = []
                if filter_resdic:
                    # Only retrieve coordinates from the filter_resdic
                    if chain in filter_resdic and resnum in filter_resdic[chain]:
                        coord_dic[identifier] = coords
                        chain_dic[chain].append(idx)
                        idx += 1
                else:
                    # retrieve everything
                    coord_dic[identifier] = coords
                    chain_dic[chain].append(idx)
                    idx += 1
    chain_ranges: ChainsRange = {}
    for chain, indice in chain_dic.items():
        if not indice:
            # this may happen when filter_resdic is defined on a different set
            # of chains
            raise ALIGNError(f"Chain matching error on {pdb_f}, chain {chain}")
        else:
            min_idx = min(indice)
            max_idx = max(indice)
            chain_ranges[chain] = (min_idx, max_idx)
    return coord_dic, chain_ranges

def cond_index(i: int, j: int, n: int) -> float:
    """
    Get the condensed index from two matrix indexes.

    Parameters
    ----------
    i : int
        Index of the first element.
    j : int
        Index of the second element.
    n : int
        Number of observations.
    """
    return n * (n - 1) / 2 - (n - i) * (n - i - 1) / 2 + j - i - 1


def get_cluster_center(npw, rmsd_matrix, n_obs):
    """
    Get the cluster centers.

    Parameters
    ----------
    npw: np.ndarray
        Indexes of the cluster over cluster_list array
    n_obs : int
        Number of overall observations (models).
    rmsd_matrix : np.ndarray
        RMSD matrix.

    Returns
    -------
    cluster_center : int
        Index of cluster center
    """
    intra_cl_distances = {el: 0.0 for el in npw}

    # iterating over the elements of the cluster
    for m_idx in range(len(npw)):
        npws = npw[m_idx + 1:]
        pairs = [int(cond_index(npw[m_idx], npw_el, n_obs)) for npw_el in npws]
        for pair_idx in range(len(pairs)):
            intra_cl_distances[npw[m_idx]] += rmsd_matrix[pairs[pair_idx]]
            intra_cl_distances[npws[pair_idx]] += rmsd_matrix[pairs[pair_idx]]
    cluster_center = min(intra_cl_distances, key=intra_cl_distances.get)  # type: ignore  # noqa : E501
    return cluster_center


def write_coords(rec_traj_filename, lig_traj_filename, topn_pdb_files, loop_residues, h3_residues):
    """
    Write coordinates to a file.

    Parameters
    ----------
    rec_traj_filename : str
        receptor trajectory filename
    lig_traj_filename : str
        ligand trajectory filename
    topn_pdb_files : list
        list of pdb files
    loop_residues : list
        list of loop residues
    h3_residues : list
        list of h3 residues
    """
    prev_keys = []
    with open(rec_traj_filename, "w") as rec_traj_xyz:
        with open(lig_traj_filename, "w") as lig_traj_xyz:
            for mod in topn_pdb_files:
                atoms = {}
                atoms.update(get_atoms(mod))
                ref_coord_dic, _ = load_coords(
                mod, atoms
                )
                if prev_keys != []:
                    if ref_coord_dic.keys() != prev_keys:
                        raise Exception(
                            "The keys of the ref_coord_dic are not the "
                            "same for all the models. Please check the "
                            "input models."
                            )
                # write receptor coords
                framework_coords = {k:ref_coord_dic[k] for k in ref_coord_dic if k[1] not in loop_residues}
                h3_coords = {k:ref_coord_dic[k] for k in ref_coord_dic if k[1] in h3_residues}
                n_rec_atoms = len(framework_coords)
                n_lig_atoms = len(h3_coords)
                rec_traj_xyz.write(f"{n_rec_atoms}{os.linesep}{os.linesep}")
                lig_traj_xyz.write(f"{n_lig_atoms}{os.linesep}{os.linesep}")
                for k, v in framework_coords.items():
                    rec_traj_xyz.write(f"{'-'.join([str(el) for el in k])} {v[0]} {v[1]} {v[2]}{os.linesep}")
                for k, v in h3_coords.items():
                    lig_traj_xyz.write(f"{'-'.join([str(el) for el in k])} {v[0]} {v[1]} {v[2]}{os.linesep}")
            prev_keys = ref_coord_dic.keys()
            return n_rec_atoms, n_lig_atoms
        

def get_clustering_dict(clusters, ligands):
    """
    Gets dictionary of clusters.

    Parameters
    ----------
    clusters : list
        list of cluster IDs
    ligands : list
        names of the ligands

    Returns
    -------
    cl_dict : dict
        dictionary of clusters
    """
    cl_dict = {}
    # loop over clusters
    for cl in range(len(clusters)):
        if clusters[cl] not in cl_dict.keys():
            cl_dict[clusters[cl]] = [ligands[cl]]
        else:
            cl_dict[clusters[cl]].append(ligands[cl])
    return cl_dict


def reres_pdb(inp_pdb_f, resid=1):
    """
    renumber residues in a pdb file
    
    Parameters
    ----------
    inp_pdb_f : Path
        input pdb file
    resid : int
        residue id to start renumbering
    """
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-reres{resid}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in renumber_residues(pdb_fh, resid):
                f.write(line)
    return out_pdb_fname


def chain_pdb(inp_pdb_f, chain_id):
    """
    alter chain in a pdb file

    Parameters
    ----------
    inp_pdb_f : Path
        input pdb file
    chain_id : str
        chain id to alter
    """
    out_pdb_fname = Path(f"{inp_pdb_f.stem}-chain{chain_id}.pdb")
    with open(inp_pdb_f, "r") as pdb_fh:
        with open(out_pdb_fname, "w") as f:
            for line in alter_chain(pdb_fh, chain_id):
                f.write(line)
    return out_pdb_fname


def mkensemble_pdb(inp_pdb_fs):
    """
    make an ensemble pdb file

    Parameters
    ----------
    inp_pdb_fs : list
        list of pdb files
    """
    out_pdb_fname = Path(f"ensemble.pdb")
    ensemble_pdb = make_ensemble(inp_pdb_fs)
    with open(out_pdb_fname, "w") as f:
            for line in ensemble_pdb:
                f.write(line)
    return out_pdb_fname