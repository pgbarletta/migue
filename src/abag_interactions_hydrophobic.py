###
# Supporting script to get the hydrophobic clusters from the AbAg complexes.
###
import numpy as np
import networkx as nx
import random
from typing import List
import MDAnalysis as mda
from itertools import combinations_with_replacement
from MDAnalysis.analysis import distances
from names import TAtom

########################################
# Some useful data and parameters
ANG_AB_ON = 0.8
ANG_AG_ON = -0.1
########################################


def get_shielding_ON(
    ag_C_xyz, ab_C_xyz, ab_ON_all, ang_ab_ON=ANG_AB_ON, ang_ag_ON=ANG_AG_ON
):
    vec_CC = ag_C_xyz - ab_C_xyz
    nvec_CC = vec_CC / np.linalg.norm(vec_CC)

    for ab_ON in ab_ON_all:
        ab_ON_xyz = ab_ON.position

        ab_vec_ON = ab_ON_xyz - ab_C_xyz
        ab_nvec_ON = ab_vec_ON / np.linalg.norm(ab_vec_ON)

        ag_vec_ON = ab_ON_xyz - ag_C_xyz
        ag_nvec_ON = ag_vec_ON / np.linalg.norm(ag_vec_ON)

        # Useful for debugging:
        dot_Cab_Cag_vs_Cab_ON = np.dot(nvec_CC, ab_nvec_ON)
        dot_Cab_Cag_vs_Cag_ON = np.dot(nvec_CC, ag_nvec_ON)
        if (dot_Cab_Cag_vs_Cab_ON > ang_ab_ON) and (dot_Cab_Cag_vs_Cag_ON < ang_ag_ON):
            return True, ab_ON
    return False, None


def get_interactions(
    ag_C, ab_close_carbons, close_polars, ang_ab_ON=ANG_AB_ON, ang_ag_ON=ANG_AG_ON
):
    ag_C_xyz = ag_C.position
    interacting_ab_C = []
    shielding_polars = []
    for ab_c in ab_close_carbons:
        are_shielded, ON_atm = get_shielding_ON(ag_C_xyz, ab_c.position, close_polars)
        if are_shielded:
            shielding_polars.append(ON_atm)
        else:
            interacting_ab_C.append(ab_c)
    return interacting_ab_C, shielding_polars


def construct_graph(pdb_filename, ag_chains, ab_chains, cutoff_carbons: float = 0.45):
    u = mda.Universe(pdb_filename)

    ag_selection = " or segid ".join(ag_chains)
    ab_selection = " or segid ".join(ab_chains)

    ag_ca_atoms = u.select_atoms(f"segid {ag_selection} and type C")
    ab_atoms = u.select_atoms(f"segid {ab_selection}")

    distancias = distances.distance_array(ag_ca_atoms.positions, ab_atoms.positions)
    ag_ca_interface, ab_interface = np.where(distancias < cutoff_carbons)

    G = nx.Graph()
    shielding_atoms: List = []

    ag_ca_set = set()
    for i, (ag_ca_idx, ab_idx) in enumerate(zip(ag_ca_interface, ab_interface)):
        # Only do this once for every ag C atom. They may show up more than once if they
        # are close to many ab atoms.
        if ag_ca_idx in ag_ca_set:
            continue
        ag_ca_set.add(ag_ca_idx)

        # Get the ag C atom.
        ag_c_atm = ag_ca_atoms[ag_ca_idx]

        # Get all the atoms that are close to this ag C atom.
        dist = distances.distance_array(ag_c_atm.position, u.atoms.positions)
        _, close_atms_idx = np.where(dist < cutoff_carbons)
        close_atoms = u.atoms[close_atms_idx]
        # Filter them out.
        ab_close_carbons = close_atoms.select_atoms(f"segid {ab_selection} and type C")
        close_polars = close_atoms.select_atoms(f"type O or type N")

        # Get the list of the non-shielded (interacting) ab C carbons.
        interacting_ab_C_atms, shielding_polar_atms = get_interactions(
            ag_c_atm, ab_close_carbons, close_polars
        )
        shielding_atoms.extend(shielding_polar_atms)

        for ab_c_atm in interacting_ab_C_atms:
            G.add_edge(ag_c_atm, ab_c_atm)

    return G


def get_putative_clusters(G):
    pre_clusteres = []
    for cluster in sorted(nx.connected_components(G), key=len, reverse=True):
        pre_clusteres.append(cluster)

    return pre_clusteres


def merge_clusters(pre_clusteres: List, cutoff_clusters: float = 0.45):
    def clusters_are_close(cluster_1, cluster_2, cutoff_clusters: float = 0.45):
        xyz_cluster_1 = np.array([atm.position for atm in cluster_1])
        xyz_cluster_2 = np.array([atm.position for atm in cluster_2])
        return np.any(
            np.where(
                distances.distance_array(xyz_cluster_1, xyz_cluster_2) < cutoff_clusters
            )
        )

    H = nx.Graph()
    H.add_node(0)
    for i, j in combinations_with_replacement(range(len(pre_clusteres)), 2):
        if i == j:
            continue
        if clusters_are_close(pre_clusteres[i], pre_clusteres[j], cutoff_clusters):
            H.add_edge(i, j)
        else:
            H.add_node(i)
            H.add_node(j)

    clusteres = []
    for connected_clusters in sorted(nx.connected_components(H), key=len, reverse=True):
        new_cluster = []
        for i in connected_clusters:
            new_cluster.extend(pre_clusteres[i])
        clusteres.append(new_cluster)

    # Make sure that the clusters are sorted by size
    idx = np.flip(np.argsort([len(c) for c in clusteres]))
    sorted_clusters = []
    for i in idx:
        sorted_clusters.append(tuple(clusteres[i]))

    return sorted_clusters


def MDAtom_to_TAtom(pdb_filename, clusters: List):
    u = mda.Universe(pdb_filename)
    backbone = u.select_atoms("backbone")
    new_clusters: List = []

    for clu in clusters:
        new_clu: List[TAtom] = []
        for c in clu:
            res = c.residue
            atom = TAtom(
                index=c.index,
                serial=c.id,
                element=c.element,
                is_sidechain=not (c in backbone),
                resSeq=res.resnum,
                resSeq_str=str(res.resnum),
                resname=res.resname,
                chain_ID=c.segid,
                chain_type="",
                CDR=0,
            )
            new_clu.append(atom)
        new_clusters.append(tuple(new_clu))

    return tuple(new_clusters)


########################################


def get_waters(topologia):

    ids_wat = []
    for at in topologia.atoms:
        if at.residue.name == "HOH":
            ids_wat.append(at.index)

    return ids_wat


def draw_clusters(
    pdb_filename, *, interface_atoms_pdb=None, ag_chains, clusters, filename
):
    with open(filename, "w") as fil:
        fil.write(f"from pymol import cmd\n\n")
        fil.write(f'cmd.set("sphere_scale", "0.9")\n\n')
        fil.write(f'cmd.load("{pdb_filename}")\n')
        fil.write(f'cmd.color("salmon", "')
        for chainID in ag_chains:
            if chainID != ".":
                fil.write(f"chain {chainID} or ")
        fil.write(f'chain {ag_chains[-1]}")\n')
        fil.write(f'cmd.color("atomic", "(not elem C)")\n\n')

        if interface_atoms_pdb:
            # Show epitope residues as lines
            epitope_residues = set(
                [atom.resSeq_str for atom in interface_atoms_pdb.antigen.values()]
            )

            for resi in epitope_residues:
                fil.write(
                    f'cmd.show("lines", "resi {resi} and chain {ag_chains[0]}")' + "\n"
                )

            # Show paratope residues as lines
            paratope_residues = set(
                [
                    (atom.resSeq_str, atom.chain_ID)
                    for atom in interface_atoms_pdb.antibody.values()
                ]
            )
            for resi, chainID in paratope_residues:
                fil.write(f'cmd.show("lines", "resi {resi} and chain {chainID}")\n')

        # Finally, show interacting carbons as spheres
        for n, cluster in enumerate(clusters):
            linea = f""
            cluster_id = "cluster_" + str(n + 1)
            fil.write(f'cmd.select("id ')
            for c in cluster:
                linea += f"{c.serial}+"
            fil.write(linea[0:-1])
            fil.write(f'")\n')
            fil.write(f'cmd.set_name("sele", "{cluster_id}")\n')
            fil.write(f'cmd.show("spheres", "{cluster_id}")\n')
            color = "%06x" % random.randint(0, 0xFFFFFF)
            fil.write(f'cmd.color("0x{color}", "{cluster_id}")\n')
            # fil.write(f'cmd.color("Gray", "{cluster_id}")\n')
