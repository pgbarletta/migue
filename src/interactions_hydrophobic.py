###
# Script to get the hydrophobic clusters from the AbAg complexes.
###
from pathlib import Path
import pickle
import logging
from typing import Final

import pandas as pd
from names import *

from abag_interactions_hydrophobic import *

data_dir = Path("/home/pbarletta/labo/22/migue/data")
expdata_dir = Path("/home/pbarletta/labo/22/migue/data/AB-Bind-Database-master")
pdbs_dir = Path("/home/pbarletta/labo/22/migue/run/pdbs")
mutpdbs_dir = Path("/home/pbarletta/labo/22/migue/run/mut_pdbs")
hydro_dir = Path("/home/pbarletta/labo/22/migue/rtdos/hydro")
CUTOFF_CLUSTERS: Final = 4.5
CUTOFF_CARBONS: Final = 4.5

if __name__ == "__main__":
    logging.basicConfig(filename=Path(hydro_dir, "log.log"), level=logging.INFO)

    print("Reading data.")
    experimental_data = Path(expdata_dir, "AB-Bind_experimental_data.csv")
    expdata_df = pd.read_csv(experimental_data, encoding="latin-1")
    pdb_list = tuple(sorted(set(expdata_df["#PDB"])))
    abag_chains = {}
    for partners, pdb_id in zip(expdata_df["Partners(A_B)"], expdata_df["#PDB"]):
        target, binder = partners.split("_")
        abag_chains[pdb_id] = Chains(antibody=binder, antigen=target)

    mut_list = []
    with open(Path(mutpdbs_dir, "mut_list.txt"), "r") as f:
        for mute in f:
            mut_list.append(mute.strip())
    print("Starting.")

    hydrophobic_clusters = {}
    # check_pdb = "HM_3BN9"
    # idx = mut_list.index(check_pdb)
    # for pdb_idcode in [mut_list[idx]]:
    for pdb_idcode in pdb_list:
        logging.info(pdb_idcode)

        pdb_filename = Path(pdbs_dir, pdb_idcode + ".pdb")
        ag_chains = abag_chains[pdb_idcode].antigen
        ab_chains = abag_chains[pdb_idcode].antibody

        try:
            G = construct_graph(pdb_filename, ag_chains, ab_chains, CUTOFF_CARBONS)
            pre_clusteres = get_putative_clusters(G)
            clusters = merge_clusters(pre_clusteres, CUTOFF_CLUSTERS)
            my_clusters = MDAtom_to_TAtom(pdb_filename, clusters)
        except Exception as e:
            logging.error(
                f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during hydrophobic "
                f"interactions calculation. Probably has no hydrophobic interactions. "
                f"Assigning an empty tuple."
            )
            hydrophobic_clusters[pdb_idcode] = ()
            continue
        else:
            hydrophobic_clusters[pdb_idcode] = my_clusters

    with open(Path.joinpath(hydro_dir, "hydrophobic.pkl"), "wb") as file:
        pickle.dump(hydrophobic_clusters, file)

    print(f" --- Done -- ")
