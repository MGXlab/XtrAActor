from Bio import SeqIO
from collections import defaultdict
import sys
from pathlib import Path
import pandas as pd

# def map_positions(cluster_num):
#     records_msa = list(SeqIO.parse(f"RF_results/foldmason/cluster_phyloglm_db_out/cluster_phyloglm_plus_ec6098ref_MSTA_aa.fa","fasta"))
# #cluster_0_1_mixed_plus_ec6098ref_MSTA_aa
#     map_dict_seq = defaultdict(dict)
#     for i in range(len(records_msa)):
#         seq_record = records_msa[i]
#         seq_id = seq_record.id
#         sequence = seq_record.seq
#         counter_seq=0
#         for i in range(len(sequence)):
#             counter_MSA=i
#             map_dict_seq[seq_id][counter_MSA] = counter_seq
#             if sequence[i] == '-':
#                 counter_seq = counter_seq
#             else:
#                 counter_seq = counter_seq + 1
#     return map_dict_seq

# if __name__ == '__main__':
#     cluster_num = sys.argv[1]
#     seq_id_in = sys.argv[2]
#     path_to_rf_results = Path(sys.argv[3])    #path to "RF_results" folder
#     path_to_feat_imp = path_to_rf_results.joinpath("Feature_importance.tsv")
    
#     # run the overall report plots
#     mapping_out_dir = path_to_rf_results.joinpath("Rep_Seq_Position_Mapping/example_phyloglm/zvalue")
#     Path.mkdir(mapping_out_dir,exist_ok=True)

#     # mapping dictionary
#     map_dict_seq = map_seq_to_msa(cluster_num)

#     # get feature importance scores
#     df_feat_ori = pd.read_csv(path_to_feat_imp, sep='\t')
#     print(df_feat_ori)
#     df_feat_ori[["msa_position", "AA"]] = df_feat_ori['index'].str.split("_",expand=True)
#     df_feat_ori_filt = df_feat_ori.loc[df_feat_ori["cluster_num"] == cluster_num]

#     for index, row in df_feat_ori_filt.iterrows():
#         # print(df_feat_ori_filt.loc[index,'position'])
#         df_feat_ori_filt.loc[index,'seq_position'] = int(map_dict_seq[seq_id_in][int(df_feat_ori_filt.loc[index,'msa_position'])])
    
#     df_feat_ori_filt.to_csv(mapping_out_dir.joinpath(f"{cluster_num}_{seq_id_in.split("|")[0]}_feature_importance.tsv"),sep='\t',index=False)

#     print(df_feat_ori_filt)


# def map_msa_to_seq():
#     records = SeqIO.to_dict(
#         SeqIO.parse(
#             "RF_results/foldmason/cluster_phyloglm_db_out/cluster_phyloglm_plus_ec6098ref_MSTA_aa.fa",
#             "fasta"
#         )
#     )
#     msa_to_seq = {}

#     for seq_id, record in records.items():
#         msa_to_seq[seq_id] = {}
#         seq_counter = 0  # 0-based ungapped sequence index

#         for msa_pos, aa in enumerate(record.seq):
#             if aa == "-":
#                 msa_to_seq[seq_id][msa_pos] = None
#             else:
#                 msa_to_seq[seq_id][msa_pos] = seq_counter
#                 seq_counter += 1

#     return msa_to_seq

def map_msa_to_seq():
    records = SeqIO.to_dict(
        SeqIO.parse(
            "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/RF_models_automation/RF_results/foldmason/cluster_0_1_mixed_db_out/cluster_0_1_mixed_MSTA_aa_plus_ec6098ref_2.faa",
            "fasta"
        )
    )
    msa_to_seq = {}
    for seq_id, record in records.items():
        msa_to_seq[seq_id] = {}
        seq_counter = 1  # count residues seen so far

        for msa_pos, aa in enumerate(record.seq):
            if aa == "-":
                msa_to_seq[seq_id][msa_pos] = None
            else:
                msa_to_seq[seq_id][msa_pos] = seq_counter
                seq_counter += 1


    return msa_to_seq

if __name__ == '__main__':

    cluster_num = sys.argv[1]
    seq_id_in = sys.argv[2]
    path_to_rf_results = Path(sys.argv[3])

    path_to_feat_imp = path_to_rf_results / "Feature_importance.tsv"

    mapping_out_dir = path_to_rf_results / "Rep_Seq_Position_Mapping/example_cluster_0_1_mixed/new_model"
    mapping_out_dir.mkdir(parents=True, exist_ok=True)

    # build mapping
    map_dict_seq = map_msa_to_seq()

    # load feature importance
    df_feat_ori = pd.read_csv(path_to_feat_imp, sep='\t')
    df_feat_ori[["msa_position", "AA"]] = df_feat_ori["index"].str.split("_", expand=True)
    df_feat_ori["msa_position"] = df_feat_ori["msa_position"].astype(int)

    df_feat_ori_filt = df_feat_ori.loc[
        df_feat_ori["cluster_num"] == cluster_num
    ].copy()

    # # invert ref_map once
    # msa_to_ref = {msa_col: ref_pos for ref_pos, msa_col in ref_map.items()}

    # df_feat_ori_filt["ec6098_pos"] = (
    #     df_feat_ori_filt["msa_position"]
    #     .map(msa_to_ref)
    # )
    # assign sequence positions safely
    df_feat_ori_filt["seq_position"] = (
        df_feat_ori_filt["msa_position"]
        .map(map_dict_seq[seq_id_in])
    )

    # drop gap-derived rows
    df_feat_ori_filt = df_feat_ori_filt.dropna(subset=["seq_position"])
    df_feat_ori_filt["seq_position"] = df_feat_ori_filt["seq_position"].astype(int)

    # write output
    out_file = mapping_out_dir / f"{cluster_num}_{seq_id_in.split('|')[0]}_feature_importance.tsv"
    df_feat_ori_filt.to_csv(out_file, sep="\t", index=False)

    print(df_feat_ori_filt)

# # --- Build reference map: EC6098 residue → MSA column ---
# def build_ref_map(msa, ref_seq_id):
#     ref_record = next(r for r in msa if r.id == ref_seq_id)
#     ref_map = {}
#     ref_pos = 0  # 1-based
#     for col_idx, aa in enumerate(ref_record.seq):
#         if aa != '-':
#             ref_map[ref_pos] = col_idx
#             ref_pos += 1
#     return ref_map


# # --- Build MSA column → sequence position map for any sequence ---
# def build_msa_to_seq_map(record):
#     msa_to_seq = {}
#     seq_pos = 0  # 1-based residue counting
#     for col_idx, aa in enumerate(record.seq):
#         if aa == "-":
#             msa_to_seq[col_idx] = None
#         else:
#             seq_pos += 1
#             msa_to_seq[col_idx] = seq_pos
#     return msa_to_seq


# if __name__ == '__main__':
#     cluster_num = sys.argv[1]
#     seq_id_in = sys.argv[2]
#     path_to_rf_results = Path(sys.argv[3])

#     path_to_feat_imp = path_to_rf_results / "Feature_importance.tsv"
#     mapping_out_dir = path_to_rf_results / "Rep_Seq_Position_Mapping/example_phyloglm/zvalue"
#     mapping_out_dir.mkdir(parents=True, exist_ok=True)

#     # --- Load MSA ---
#     msa_file = "RF_results/foldmason/cluster_phyloglm_db_out/cluster_phyloglm_plus_ec6098ref_MSTA_aa.fa"
#     msa = list(SeqIO.parse(msa_file, "fasta"))

#     # --- Build reference mapping ---
#     ref_map = build_ref_map(msa, seq_id_in)  # EC6098 residues → MSA columns
#     msa_to_ref = {msa_col: ref_pos for ref_pos, msa_col in ref_map.items()}  # invert

#     # --- Build target sequence MSA → residue mapping ---
#     target_record = next(r for r in msa if r.id == seq_id_in)
#     msa_to_seq = build_msa_to_seq_map(target_record)  # MSA col → target residue

#     # --- Load feature importance ---
#     df_feat_ori = pd.read_csv(path_to_feat_imp, sep="\t")
#     df_feat_ori[["msa_position", "AA"]] = df_feat_ori["index"].str.split("_", expand=True)
#     df_feat_ori["msa_position"] = df_feat_ori["msa_position"].astype(int)

#     # Filter to cluster
#     df_feat_ori_filt = df_feat_ori.loc[df_feat_ori["cluster_num"] == cluster_num].copy()

#     # --- Map MSA columns to reference and target positions ---
#     # EC6098 (reference) positions
#     df_feat_ori_filt["ec6098_pos"] = df_feat_ori_filt["msa_position"].map(msa_to_ref)

#     # Target sequence positions
#     df_feat_ori_filt["seq_position"] = df_feat_ori_filt["msa_position"].map(msa_to_seq)

#     # Drop positions that correspond to gaps in target sequence
#     df_feat_ori_filt = df_feat_ori_filt.dropna(subset=["seq_position"])
#     df_feat_ori_filt["seq_position"] = df_feat_ori_filt["seq_position"].astype(int)

#     # --- Write output ---
#     out_file = mapping_out_dir / f"{cluster_num}_{seq_id_in.split('|')[0]}_feature_importance.tsv"
#     df_feat_ori_filt.to_csv(out_file, sep="\t", index=False)

#     print(df_feat_ori_filt)


