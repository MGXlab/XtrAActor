#!usr/bin/python3

"""
updates on this script:
- removing zero variance features
- added a 5-fold cross validation from which the average accuracy and average gini scores are calculated and reported
- removed shap calculations - do we even use these?

"""

# libraries
import sys
import re
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


def parse_seqs_and_metadata(fasta_path, metadata_path):
    seq_dict = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = record.seq

    # for the structural MSA
    data_df = pd.DataFrame.from_dict(seq_dict, orient='index')

    data_df_numerical = pd.get_dummies(data_df, dtype=int)
    data_df_numerical = data_df_numerical.loc[:, ~data_df_numerical.columns.str.contains('-')]  # remove position-gap columns

    data_df_numerical.reset_index(names='protein', inplace=True)
    data_df_numerical.loc[data_df_numerical['protein'].str.contains('IMG'), 'img_split'] = (
        data_df_numerical.loc[data_df_numerical['protein'].str.contains('IMG'), 'protein']
        .str.split('|').str[0]
    )
    data_df_numerical['img_split'] = data_df_numerical['img_split'].fillna(data_df_numerical['protein'])
    data_df_numerical.set_index('img_split', inplace=True)

    all_feature_names = data_df_numerical.columns[1:]  # everything except "protein"
    data_df_numerical, feature_names = remove_zero_variance_features(data_df_numerical, all_feature_names)
    num_features = len(feature_names)

    biomes_metadata_df = pd.read_csv(metadata_path, header=0, sep=',')
    biomes_metadata_df.drop(columns=["cluster_rep"], inplace=True)
    # guard against duplicate protein rows in the metadata inflating the merge
    biomes_metadata_df.drop_duplicates(subset='protein', keep='first', inplace=True)

    # merge the dataframes
    data_df_merged = data_df_numerical.join(biomes_metadata_df.set_index('protein'), how='left')
    data_df_merged_1 = data_df_merged.reset_index().drop(columns='img_split').set_index('protein')

    # drop rows with no biome label now, so X and y stay aligned downstream
    # (previously this was only done to y, which desyncs X and y and crashes train_test_split)
    data_df_merged_1 = data_df_merged_1.dropna(subset=['ecosystem_subtype'])

    if 'ecosystem_type' in data_df_merged_1.columns:
        temp = data_df_merged_1.groupby('ecosystem_type').count()  # sample counts per biome
        print(f"Number of samples in each biome:\n{temp['ecosystem']}")

    return data_df_merged_1, num_features, feature_names


def run_the_classifier(data_df_merged, num_features, feature_names):

    X = data_df_merged.iloc[:, 0:num_features]           # features
    y_labels = data_df_merged['ecosystem_subtype']        # categorical labels, for stratification
    y = pd.get_dummies(y_labels, dtype=int)                # one-hot, for fitting

    # --- HYPERPARAMETER SEARCH ---
    # done once, on a single 80/20 split, purely to pick hyperparameters
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, train_size=0.8, test_size=0.2, stratify=y, random_state=42
    )

    param_grid = {
        'max_depth': [5, 10, None],
        'max_features': ['sqrt', 'log2'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [1000, 2000, 5000],
        'class_weight': ['balanced', None]}

    random_search = RandomizedSearchCV(
        RandomForestClassifier(random_state=42), param_grid, random_state=42
    )
    random_search.fit(X_train, y_train)
    best_params = random_search.best_estimator_
    print(f"Hyperparameter tuning... best parameters found: {best_params}")

    # --- 5-FOLD CV: accuracy, report, conf_matrix, AND gini importance ---
    # accuracy/report/conf_matrix are computed from pooled OUT-OF-FOLD predictions:
    # every sample is scored only by a fold that did NOT train on it, so this
    # is a genuine held-out/test-data accuracy, just averaged across 5 folds
    # instead of a single fixed 80/20 split.
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    fold_importances = []   # one Series per fold
    fold_accuracies = []    # one float per fold
    oof_true = []            # out-of-fold true labels (one-hot), pooled across folds
    oof_pred = []            # out-of-fold predictions, pooled across folds

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y_labels)):
        X_fold_train, X_fold_test = X.iloc[train_idx], X.iloc[test_idx]
        y_fold_train, y_fold_test = y.iloc[train_idx], y.iloc[test_idx]

        fold_clf = RandomForestClassifier(**best_params.get_params())
        fold_clf.set_params(random_state=42 + fold_idx)
        fold_clf.fit(X_fold_train, y_fold_train)

        fold_pred = fold_clf.predict(X_fold_test)
        fold_acc = metrics.accuracy_score(y_fold_test, fold_pred)
        fold_accuracies.append(fold_acc)
        print(f"  fold {fold_idx} accuracy: {fold_acc:.4f}")

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
    accuracy = float(np.mean(fold_accuracies))
    accuracy_std = float(np.std(fold_accuracies))
    print(f"Mean CV accuracy: {accuracy:.4f} (+/- {accuracy_std:.4f})")

    y_test_oof = np.vstack(oof_true)   # pooled out-of-fold true labels, one-hot
    y_pred_oof = np.vstack(oof_pred)   # pooled out-of-fold predictions, one-hot

    report = metrics.classification_report(
        y_test_oof, y_pred_oof, output_dict=True, target_names=y.columns.tolist()
    )
    conf_matrix = confusion_matrix(y_test_oof.argmax(axis=1), y_pred_oof.argmax(axis=1))
    print("Pooled out-of-fold classification report:",
          metrics.classification_report(y_test_oof, y_pred_oof, target_names=y.columns.tolist()))
    print("Pooled out-of-fold confusion matrix:\n", conf_matrix)

    print(f"Number of features: {num_features}")

    return best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df


def create_output_files(data_df_merged, fasta_path, output_dir,
                         best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df):

    # cluster name pulled robustly from the fasta filename (cluster_N),
    # rather than a fixed-length slice that breaks on filename-length changes
    match = re.search(r'(cluster_\d+)', fasta_path.stem)
    clust_num = match.group(1) if match else fasta_path.stem

    features_out_path = output_dir.joinpath(f"{clust_num}_Feature_importance.tsv")
    report_out_path = output_dir.joinpath(f"{clust_num}_RandomForest_Report.tsv")
    data_df_out_path = output_dir.joinpath(f"{clust_num}_data_df.tsv")

    # --- CREATE AN OUTPUT FILE FOR RF MODEL REPORT ---
    # report's keys are now the actual biome class names (target_names was
    # passed to classification_report), not "0"/"1" — so iterate over the
    # real per-class entries rather than a hardcoded range(0, 2)
    excluded_keys = {'accuracy', 'macro avg', 'weighted avg'}
    class_names = [k for k in report.keys() if k not in excluded_keys]

    report_dict = {}
    for class_name in class_names:
        precision = report[class_name]['precision']
        recall = report[class_name]['recall']
        f1 = report[class_name]['f1-score']
        support = report[class_name]['support']
        report_dict[class_name] = [
            clust_num, best_params, accuracy, accuracy_std,
            precision, recall, f1, support, conf_matrix
        ]

    report_df = pd.DataFrame.from_dict(
        report_dict, orient='index',
        columns=['clust_num', 'RF_parameters', 'accuracy', 'accuracy_std',
                 'precision', 'recall', 'f1', 'support', 'conf_matrix']
    ).reset_index(names='biome')
    report_df.to_csv(report_out_path, sep='\t', mode='a', index=False,
                      header=not report_out_path.exists())

    # --- CREATE AN OUTPUT FILE FOR THE IMPORTANT FEATURES ---
    # feature | gini_imp | cluster_num, where gini_imp is the mean gini
    # importance across the 5 CV folds
    features_df = (
        gini_folds_df[['gini_mean']]
        .rename(columns={'gini_mean': 'gini_imp'})
        .reset_index(names='feature')
    )
    features_df['cluster_num'] = clust_num
    print(features_df)

    features_df.to_csv(features_out_path, sep='\t', mode='a', index=False,
                        header=not features_out_path.exists())

    # --- CREATE AN OUTPUT FILE FOR DATA DF ---
    data_df_merged.to_csv(data_df_out_path, sep='\t')

    return


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: python3 {sys.argv[0]} <fasta_path> <metadata_path>")

    fasta_path = Path(sys.argv[1])
    metadata_path = Path(sys.argv[2])
    output_dir = Path("RF_results")
    output_dir.mkdir(parents=True, exist_ok=True)

    data_df_merged, num_features, feature_names = parse_seqs_and_metadata(fasta_path, metadata_path)

    best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df = run_the_classifier(
        data_df_merged, num_features, feature_names
    )

    create_output_files(
        data_df_merged, fasta_path, output_dir,
        best_params, accuracy, accuracy_std, report, conf_matrix, gini_folds_df
    )

    print(f"{fasta_path.stem} finished... Model accuracy: {accuracy:.4f} (+/- {accuracy_std:.4f})")