# Drug-Synergy-KG
Drug Synergy via KG (Knowledge Graph) with intrinsic features as initial embeddings

The critical raw training data for this project include `comb_final.txt`, `drug_init_embs.csv`, `init_copy_number.csv`, `kg_final2.txt`. The key scripts are `dataloader4KGNN.py`, `load_gat.py`, `model.py`, `run.py`, with utilized components located in the `utils` directory. Within this model, data from the cell line `U251MG` is used in place of the cell line `SNB19`, with `U251MG` being the parental cell line of `SNB19`.
