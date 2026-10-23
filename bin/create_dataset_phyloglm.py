import sys
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Create dataset for PhyloGLM")
parser.add_argument("feature_df", help="Path to the feature.tsv file")
parser.add_argument("biome_file", help="Path to the data_df.tsv file (RF model output: gene, features..., biome)")
parser.add_argument("-o", "--output", default="output_dataset.tsv", help="Output file path")
parser.add_argument("-r", "--reference_biome", default="Lake",
                     help="Last ';'-segment of the biome string to treat as the baseline/reference "
                          "class (e.g. 'Lake'); whichever of the two biomes present does NOT end in "
                          "this becomes the phyloglm target/response (default: Lake)")
args = parser.parse_args()

feat_imp_df = pd.read_csv(args.feature_df, sep='\t', usecols=[0,1],names=['col_name','feature_importance_vals'], dtype={'col_name':'str','feature_importance_vals':np.float64},header=0) 
#print(feat_imp_df.shape)
x = len(feat_imp_df)
top_50_ft_df = feat_imp_df.head(n=x)
#top_50_ft_df

# feat_imp_df.loc[feat_imp_df['feature_importance_vals'] >= 0.00001]

df_data = pd.read_csv(args.biome_file, sep='\t')
#print(df_data)
df_data.drop_duplicates('gene',inplace=True)

df_biome = df_data[["gene","biome"]]
# df_biome.to_csv("comb_cluster_biome.tsv", sep='\t', index=False)
df_biome_num = pd.get_dummies(df_biome.set_index("gene"), dtype=int)
#print(df_biome_num)

# get all values of biome (full ";"-delimited strings), choose the one whose
# LAST segment is not the reference biome (e.g. not "Lake"). biome is now the
# full hierarchy path, so we compare against its last segment rather than the
# whole string — this expects exactly 2 biome classes per group, since
# phyloglm needs a binary response.
unique_biomes = df_biome['biome'].dropna().unique()
if len(unique_biomes) != 2:
    sys.exit(
        f"Expected exactly 2 biome classes for a pairwise phyloglm comparison, "
        f"found {len(unique_biomes)}: {list(unique_biomes)}"
    )

def last_segment(biome_str):
    return biome_str.split(';')[-1]

non_reference_biomes = [b for b in unique_biomes if last_segment(b) != args.reference_biome]
if len(non_reference_biomes) != 1:
    sys.exit(
        f"Could not find exactly one non-reference biome (reference='{args.reference_biome}'). "
        f"Biomes present: {list(unique_biomes)}"
    )
target_biome = non_reference_biomes[0]
target_col = f"biome_{target_biome}"

df_data_filtered = df_data[["gene"] + [col for col in top_50_ft_df['col_name'] if col in df_data.columns]]
df_data_filtered_merged = df_biome_num[[target_col]].join(df_data_filtered.set_index("gene"),how="left")
df_data_filtered_merged.reset_index(inplace=True)
#print(df_data_filtered_merged)

# save look up table with feature and gini
gini_lookup = top_50_ft_df.rename(columns={'col_name': 'feature', 'feature_importance_vals': 'gini'})
gini_output = args.output.replace('_phyloglm_dataset.tsv', '_gini.tsv')
gini_lookup.to_csv(gini_output, sep='\t', index=False)

df_data_filtered_merged.to_csv(args.output, sep='\t', index=False)