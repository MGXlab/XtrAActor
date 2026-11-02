#!usr/bin/python3

"""
updates on this script:
- removing zero variance features
- added a 5-fold cross validation from which the average accuracy and average gini scores are calculated and reported
- removed shap calculations - do we even use these?
- runs per-group again, called from the nextflow RF_MODEL process:
    python3 rf_model.py <group_XX_msa.faa> <group_XX_metadata.tsv> <group_XX>
  fasta and metadata are the per-group MSA/metadata files nextflow stages in
  (fasta_dir/group_XX.faa -> mafft/trimal -> metadata_dir/group_XX_metadata.tsv);
  the third arg is the group name, used for output file naming
- classification target is the FULL biome string (e.g.
  "Environmental;Aquatic;Freshwater;Lake"), not just the ecosystem_subtype
- gene name is the full fasta header up to the first whitespace, and must
  match the metadata's gene column exactly (no partial/fallback matching)
- RF training is still parallelized (n_jobs=-1) in both the hyperparameter
  search and the per-fold models
"""


# libraries
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn import metrics
from sklearn.metrics import confusion_matrix


def remove_zero_variance_features(data_df_numerical, feature_names):
    """
    Drop one-hot feature columns where every sequence has the same value
    (variance == 0) — these can never help split the data and just add
    noise/slowdown to the model.
    """
    variances = data_df_numerical[feature_names].var()
    zero_var_features = variances[variances == 0].index.tolist()

    if zero_var_features:
        print(f"Removing {len(zero_var_features)} zero-variance feature(s) "
              f"out of {len(feature_names)}")
        data_df_numerical = data_df_numerical.drop(columns=zero_var_features)

    remaining_feature_names = pd.Index([f for f in feature_names if f not in zero_var_features])
    return data_df_numerical, remaining_feature_names


MAX_RATIO = 3.0   # majority:minority ratio after downsampling (1.0 = fully balanced, 2.0 = at most 2:1)

def downsample_majority(data_df_merged, label_col='biome', random_state=42):
    """
    Randomly downsample larger biome classes so no class has more than
    MAX_RATIO x the members of the smallest class. Used for RF training only;
    the full data_df is still written out for phyloglm.
    """
    counts = data_df_merged[label_col].value_counts()
    n_max = int(counts.min() * MAX_RATIO)
    if counts.max() <= n_max:
        return data_df_merged

    balanced = (
        data_df_merged.groupby(label_col, group_keys=False)
        .apply(lambda g: g.sample(n=min(len(g), n_max), random_state=random_state))
    )
    print(f"Downsampling for RF: {counts.to_dict()} -> "
          f"{balanced[label_col].value_counts().to_dict()}")
    return balanced

def parse_seqs_and_metadata(fasta_path, metadata_path):
    """
    fasta_path: per-group MSA fasta file (e.g. group_XX_mafft_aa.fa), one
        record per gene. Gene name = full header text up to the first
        whitespace (record.id).
    metadata_path: per-group tsv with two columns, no header: gene, biome.
        The metadata's 'gene' column must match the fasta gene name exactly.
        The FULL biome string is used as the classification label.
    """
    seq_dict = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = record.seq  # record.id = everything up to the first whitespace

    # for the structural MSA
    data_df = pd.DataFrame.from_dict(seq_dict, orient='index')

    data_df_numerical = pd.get_dummies(data_df, dtype=int)
    data_df_numerical = data_df_numerical.loc[:, ~data_df_numerical.columns.str.contains('-')]  # remove position-gap columns

    data_df_numerical.index.name = 'gene'

    all_feature_names = data_df_numerical.columns  # every column is a feature; index (gene) isn't a column
    data_df_numerical, feature_names = remove_zero_variance_features(data_df_numerical, all_feature_names)
    num_features = len(feature_names)

    biomes_metadata_df = pd.read_csv(metadata_path, sep='\t', usecols=['gene', 'biome'])
    biomes_metadata_df['gene'] = biomes_metadata_df['gene'].str.strip()
    # normalize biome delimiter — some metadata files use ":" instead of ";";
    # without this, "A:B:C" and "A;B;C" would be treated as different classes.
    # (.str.strip()/.str.replace() leave real NaNs untouched, unlike astype(str),
    # so missing biomes still get caught by the dropna(subset=['biome']) below)
    biomes_metadata_df['biome'] = biomes_metadata_df['biome'].str.strip().str.replace(":", ";")
    # guard against duplicate gene rows in the metadata inflating the merge
    biomes_metadata_df.drop_duplicates(subset='gene', keep='first', inplace=True)

    # --- sanity check: the fasta gene names and the metadata 'gene' column
    # are expected to match exactly. Fail loudly with examples from both
    # sides if they don't, rather than downstream with an opaque sklearn
    # error or a silently mismatched merge. ---
    fasta_genes = set(data_df_numerical.index)
    metadata_genes = set(biomes_metadata_df['gene'])
    overlap = fasta_genes & metadata_genes
    print(f"Fasta genes: {len(fasta_genes)} | Metadata genes: {len(metadata_genes)} | Overlap: {len(overlap)}")

    if overlap != fasta_genes:
        missing = sorted(fasta_genes - metadata_genes)
        print(f"{len(missing)} fasta gene(s) have no matching metadata entry. Examples: {missing[:5]}")
        print("Example fasta gene ids:", list(fasta_genes)[:5])
        print("Example metadata gene ids:", list(metadata_genes)[:5])
        sys.exit(
            "Fasta gene ids and metadata gene ids don't match exactly. "
            "Compare the example ids printed above (case, whitespace, any "
            "prefix/suffix) and fix the mismatch at the source."
        )

    # merge the dataframes
    data_df_merged = data_df_numerical.join(biomes_metadata_df.set_index('gene'), how='left')

    # drop rows with no biome label now, so X and y stay aligned downstream
    data_df_merged = data_df_merged.dropna(subset=['biome'])
    data_df_merged = data_df_merged[data_df_merged['biome'] != 'unknown']
    # --- drop biome classes too small to stratify ---
    # StratifiedKFold(n_splits=5) needs >=5 members per class, and the initial
    # 80/20 train_test_split needs >=2. Rather than crashing the whole task
    # (which kills the group's results entirely once retries are exhausted),
    # drop under-represented classes and keep going with what's left.
    MIN_CLASS_SIZE = 5
    biome_counts = data_df_merged['biome'].value_counts()
    small_biomes = biome_counts[biome_counts < MIN_CLASS_SIZE].index.tolist()
    if small_biomes:
        print(f"Dropping {len(small_biomes)} biome(s) with < {MIN_CLASS_SIZE} samples "
              f"(can't be stratified across 5 folds): {small_biomes}")
        data_df_merged = data_df_merged[~data_df_merged['biome'].isin(small_biomes)]

    remaining_biomes = data_df_merged['biome'].value_counts()
    if len(remaining_biomes) < 2:
        sys.exit(
            f"After dropping small biome classes, only {len(remaining_biomes)} "
            f"class(es) remain — need at least 2 to train a classifier. "
            f"Remaining counts:\n{remaining_biomes}"
        )

    print(f"Number of samples in each biome:\n{remaining_biomes}")

    return data_df_merged, num_features, feature_names


def run_the_classifier(data_df_merged, num_features, feature_names):
    """
    In order to use as much data as possible from each biome, 
    class balancing is done INSIDE each split. The training part is
    downsampled (see downsample_majority / MAX_RATIO); every test part keeps
    the real class mix (not downsampled - ratio of the real data stays). 
    Performance is reported as balanced accuracy.
    """

    X = data_df_merged.iloc[:, 0:num_features]     # features
    y_labels = data_df_merged['biome']              # full biome string: used for fitting AND stratification
    classes = sorted(y_labels.unique())             # fixed class order for report / confusion matrix

    # --- HYPERPARAMETER SEARCH ---
    # done once, on a single 80/20 split, purely to pick hyperparameters;
    # the 80 % training part is balanced, the search scores by balanced accuracy
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_labels, train_size=0.8, test_size=0.2, stratify=y_labels, random_state=42
    )
    train_bal = downsample_majority(data_df_merged.loc[X_train.index])
    X_train_bal, y_train_bal = X.loc[train_bal.index], y_labels.loc[train_bal.index]

    param_grid = {
        'max_depth': [5, 10, None],
        'max_features': ['sqrt', 'log2'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [1000, 2000, 5000],
        'class_weight': ['balanced', None]}

    # inner CV of the search needs >= n_splits members per class
    inner_cv = max(2, min(5, int(y_train_bal.value_counts().min())))

    # n_jobs=-1 here parallelizes the search across hyperparameter
    # candidates/CV folds
    random_search = RandomizedSearchCV(
        RandomForestClassifier(random_state=42), param_grid,
        scoring='balanced_accuracy', cv=inner_cv, random_state=42, n_jobs=-1
    )
    random_search.fit(X_train_bal, y_train_bal)
    best_params = random_search.best_estimator_
    print(f"Hyperparameter tuning... best parameters found: {best_params}")

    # --- 5-FOLD CV: balanced accuracy, report, conf_matrix, AND gini importance ---
    # each fold: training part downsampled (different random subset per fold),
    # test part = all held-out sequences with the real class mix.
    # report/conf_matrix come from pooled OUT-OF-FOLD predictions, so every
    # sequence is scored exactly once by a model that did not train on it.
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    fold_importances = []   # one Series per fold
    fold_accuracies = []    # one float per fold (balanced accuracy)
    oof_true = []           # out-of-fold true labels, pooled across folds
    oof_pred = []           # out-of-fold predictions, pooled across folds

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y_labels)):
        # balance the training part of this fold only
        fold_train_bal = downsample_majority(data_df_merged.iloc[train_idx],
                                             random_state=42 + fold_idx)
        X_fold_train = X.loc[fold_train_bal.index]
        y_fold_train = y_labels.loc[fold_train_bal.index]

        # test part keeps the real class mix
        X_fold_test = X.iloc[test_idx]
        y_fold_test = y_labels.iloc[test_idx]

        fold_clf = RandomForestClassifier(**best_params.get_params())
        # n_jobs=-1 here parallelizes tree building for the final per-fold model
        fold_clf.set_params(random_state=42 + fold_idx, n_jobs=-1)
        fold_clf.fit(X_fold_train, y_fold_train)

        fold_pred = fold_clf.predict(X_fold_test)
        fold_acc = metrics.balanced_accuracy_score(y_fold_test, fold_pred)
        fold_accuracies.append(fold_acc)
        print(f"  fold {fold_idx} balanced accuracy: {fold_acc:.4f} "
              f"(trained on {len(fold_train_bal)}, tested on {len(test_idx)})")

        oof_true.append(y_fold_test.values)
        oof_pred.append(fold_pred)

        fold_imp = pd.Series(fold_clf.feature_importances_, index=feature_names, name=f"fold_{fold_idx}")
        fold_importances.append(fold_imp)

    # --- aggregate Gini importances across folds ---
    gini_folds_df = pd.concat(fold_importances, axis=1)  # columns: fold_0 ... fold_4
    fold_cols = [c for c in gini_folds_df.columns if c.startswith('fold_')]
    gini_folds_df['gini_mean'] = gini_folds_df[fold_cols].mean(axis=1)
    gini_folds_df['gini_std'] = gini_folds_df[fold_cols].std(axis=1)
    gini_folds_df = gini_folds_df.sort_values('gini_mean', ascending=False)

    # --- aggregate accuracy/report/conf_matrix across folds (pooled OOF) ---
    # 'accuracy' is the mean BALANCED accuracy across folds
    accuracy = float(np.mean(fold_accuracies))
    accuracy_std = float(np.std(fold_accuracies))
    print(f"Mean CV balanced accuracy: {accuracy:.4f} (+/- {accuracy_std:.4f})")

    y_test_oof = np.concatenate(oof_true)   # pooled out-of-fold true labels
    y_pred_oof = np.concatenate(oof_pred)   # pooled out-of-fold predictions

    report = metrics.classification_report(
        y_test_oof, y_pred_oof, labels=classes, output_dict=True, zero_division=0
    )
    conf_matrix = confusion_matrix(y_test_oof, y_pred_oof, labels=classes)
    print("Pooled out-of-fold classification report:",
          metrics.classification_report(y_test_oof, y_pred_oof, labels=classes, zero_division=0))
    print(f"Pooled out-of-fold confusion matrix (rows = true, cols = predicted, order {classes}):\n",
          conf_matrix)

    print(f"Number of features: {num_features}")

    return best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df


def create_output_files(data_df_merged, run_name, output_dir,
                         best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df):

    features_out_path = output_dir.joinpath(f"{run_name}_Feature_importance.tsv")
    report_out_path = output_dir.joinpath(f"{run_name}_RandomForest_Report.tsv")
    data_df_out_path = output_dir.joinpath(f"{run_name}_data_df.tsv")

    # --- CREATE AN OUTPUT FILE FOR RF MODEL REPORT ---
    # report's keys are the actual biome strings (target_names was passed to
    # classification_report), so iterate over the real per-class entries
    excluded_keys = {'accuracy', 'macro avg', 'weighted avg'}
    class_names = [k for k in report.keys() if k not in excluded_keys]

    report_dict = {}
    for class_name in class_names:
        precision = report[class_name]['precision']
        recall = report[class_name]['recall']
        f1 = report[class_name]['f1-score']
        support = report[class_name]['support']
        report_dict[class_name] = [
            run_name, best_params, accuracy, accuracy_std,
            precision, recall, f1, support, conf_matrix
        ]

    report_df = pd.DataFrame.from_dict(
        report_dict, orient='index',
        columns=['run_name', 'RF_parameters', 'accuracy', 'accuracy_std',
                 'precision', 'recall', 'f1', 'support', 'conf_matrix']
    ).reset_index(names='biome')
    report_df.to_csv(report_out_path, sep='\t', mode='a', index=False,
                      header=not report_out_path.exists())

    # --- CREATE AN OUTPUT FILE FOR THE IMPORTANT FEATURES ---
    # feature | gini_imp | run_name, where gini_imp is the mean gini
    # importance across the 5 CV folds
    features_df = (
        gini_folds_df[['gini_mean']]
        .rename(columns={'gini_mean': 'gini_imp'})
        .reset_index(names='feature')
    )
    features_df['run_name'] = run_name
    print(features_df)

    features_df.to_csv(features_out_path, sep='\t', mode='a', index=False,
                        header=not features_out_path.exists())

    # --- CREATE AN OUTPUT FILE FOR DATA DF ---
    data_df_merged.to_csv(data_df_out_path, sep='\t')

    return


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(f"Usage: python3 {sys.argv[0]} <fasta_path> <metadata_path> <group_name>")

    fasta_path = Path(sys.argv[1])
    metadata_path = Path(sys.argv[2])
    group_name = sys.argv[3]

    if not fasta_path.exists():
        sys.exit(f"Could not find {fasta_path}")
    if not metadata_path.exists():
        sys.exit(f"Could not find {metadata_path}")

    output_dir = Path("RF_results")
    output_dir.mkdir(parents=True, exist_ok=True)

    data_df_merged, num_features, feature_names = parse_seqs_and_metadata(fasta_path, metadata_path)

    best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df = run_the_classifier(
        data_df_merged, num_features, feature_names
    )

    create_output_files(
        data_df_merged, group_name, output_dir,
        best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df
    )

    print(f"{group_name} finished... Model accuracy: {accuracy:.4f} (+/- {accuracy_std:.4f})")