import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Create iTOL colorstrip dataset for a cluster tree")
parser.add_argument("data_df", help="Path to the data_df.tsv file")
parser.add_argument("-o", "--output", default="dataset_color_strip.txt", help="Output iTOL dataset file")
parser.add_argument("--label", default="Biome", help="Dataset label shown in iTOL legend")
args = parser.parse_args()

# Feste Farbzuordnung pro Biom-Kategorie -- hier anpassen/erweitern je nach euren tatsächlichen Werten
color_map = {
    "Lake": "#1f77b4",
    "Oceanic": "#d62728",
    "Hot (42-90C)": "#ff7f0e",
    # weitere Kategorien nach Bedarf ergänzen
}

df = pd.read_csv(args.data_df, sep='\t', low_memory=False)
df = df.drop_duplicates('protein')

header = f"""DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL {args.label}
COLOR #ff0000
COLOR_BRANCHES 0

LEGEND_TITLE {args.label}
LEGEND_SHAPES {' '.join(['1'] * len(color_map))}
LEGEND_COLORS {' '.join(color_map.values())}
LEGEND_LABELS {' '.join(color_map.keys())}

DATA
"""

lines = []
for _, row in df.iterrows():
    biome = row['ecosystem_subtype']
    color = color_map.get(biome, "#808080")  # grau als Fallback für unbekannte Kategorien
    lines.append(f"{row['protein']} {color} {biome}")

with open(args.output, "w") as f:
    f.write(header)
    f.write("\n".join(lines))
    f.write("\n")