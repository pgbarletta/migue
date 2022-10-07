###
# Script to get the hydrophobic clusters from the AbAg complexes.
###
from asyncio import as_completed
from pathlib import Path
import pickle
import logging
from typing import Final, Dict, Tuple
import concurrent.futures as cf
from collections import defaultdict

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


def get_hydro(
    pdb_idcode: str, pdb_filename: Path, ag_chains, ab_chains
) -> Tuple[str, Tuple]:
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
        return pdb_idcode, tuple()
    else:
        return pdb_idcode, my_clusters


if __name__ == "__main__":
    if not hydro_dir.is_dir():
        hydro_dir.mkdir()

    logging.basicConfig(filename=Path(hydro_dir, "log.log"), level=logging.INFO)

    print("Reading data.")
    experimental_data = Path(expdata_dir, "AB-Bind_experimental_data.csv")
    expdata_df = pd.read_csv(experimental_data, encoding="latin-1")
    pdb_list = list(sorted(set(expdata_df["#PDB"])))
    pdb_set = set(pdb_list)

    abag_chains = {}
    for partners, pdb_id in zip(expdata_df["Partners(A_B)"], expdata_df["#PDB"]):
        target, binder = partners.split("_")
        abag_chains[pdb_id] = Chains(antibody=binder, antigen=target)

    mut_list = []
    with open(Path(mutpdbs_dir, "mut_list.txt"), "r") as f:
        for mute in f:
            mut_list.append(mute.strip())

    pdbs_mut = defaultdict(list)
    for muta in mut_list:
        this_pdb = muta.split("-")[0]
        pdbs_mut[this_pdb].append(muta)

    print("Starting.")

    mut_hydrophobic_clusters: Dict[str, Tuple] = {}
    pdb_hydrophobic_clusters: Dict[str, Tuple] = {}
    # check_pdb = "HM_3BN9"
    # idx = mut_list.index(check_pdb)
    # for pdb_idcode in [mut_list[idx]]:
    futuros = []
    with cf.ProcessPoolExecutor(4) as ex:
        for pdb, mutas in pdbs_mut.items():
            logging.info(f"{pdb} launched.")

            pdb_filename = Path(pdbs_dir, pdb + ".pdb")
            ag_chains = abag_chains[pdb].antigen
            ab_chains = abag_chains[pdb].antibody
            futuros.append(
                ex.submit(
                    get_hydro,
                    pdb,
                    pdb_filename,
                    ag_chains,
                    ab_chains,
                )
            )
            for mut in mutas:
                pdb_filename = Path(mutpdbs_dir, mut + ".pdb")
                futuros.append(
                    ex.submit(get_hydro, mut, pdb_filename, ag_chains, ab_chains)
                )
        for futu in cf.as_completed(futuros):
            if futu.exception():
                logging.error(f"{futu.exception()}.")
                continue
            pdb_idcode, my_clusters = futu.result()
            if pdb_idcode in pdb_set:
                pdb_hydrophobic_clusters[pdb_idcode] = my_clusters
            else:
                mut_hydrophobic_clusters[pdb_idcode] = my_clusters

            logging.info(f"{pdb_idcode} done.")

    with open(Path.joinpath(hydro_dir, "pdb_hydrophobic.pkl"), "wb") as file:
        pickle.dump(pdb_hydrophobic_clusters, file)
    with open(Path.joinpath(hydro_dir, "mut_hydrophobic.pkl"), "wb") as file:
        pickle.dump(mut_hydrophobic_clusters, file)

    print(f" --- Done -- ")
