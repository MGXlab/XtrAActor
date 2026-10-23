#!usr/bin/python3

# libraries
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV 
from sklearn import metrics
from sklearn.metrics import confusion_matrix
import shap


def parse_seqs_and_metadata(fasta_path, metadata_path):
    seq_dict = {}
    print(fasta_path)

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_dict[record.id] = record.seq

    # convert to a dataframe (one row per protein, one column per alignment position)
    data_df = pd.DataFrame.from_dict(seq_dict, orient='index')

    data_df_numerical = pd.get_dummies(data_df, dtype=int)
    data_df_numerical = data_df_numerical.loc[:, ~data_df_numerical.columns.str.contains('-')]  # remove position-gap columns

    # capture feature names/count BEFORE adding any non-feature columns
    num_features = data_df_numerical.columns.shape[0]
    feature_names = data_df_numerical.columns

    # fasta headers look like "SRR5209028rbin.34|KIEKOPAJ_00675|group_00009"
    # -> isolate is everything before the first "|"
    data_df_numerical.reset_index(names='protein', inplace=True)
    data_df_numerical['isolate'] = data_df_numerical['protein'].str.split('|').str[0]
    data_df_numerical.set_index('protein', inplace=True)
    # print(data_df_numerical)

    # metadata columns: isolate, Soil.type.host.species, Compartment, Data, Label
    metadata_df = pd.read_csv(metadata_path, header=0)

    # make sure both join keys are the same dtype/format before merging
    metadata_df['isolate'] = metadata_df['isolate'].astype(str)
    data_df_numerical['isolate'] = data_df_numerical['isolate'].astype(str)

    metadata_df.set_index('isolate', inplace=True)

    # map isolate -> Compartment directly
    isolate_to_comp = metadata_df['Compartment']  # Series indexed by isolate
    data_df_numerical['Comp'] = data_df_numerical['isolate'].map(isolate_to_comp)

    data_df_merged = data_df_numerical.drop(columns='isolate')

    temp = data_df_merged.groupby('Comp').size()  # count of samples per class
    print(f"Number of samples in each biome:\n{temp}")

    n_unmatched = data_df_merged['Comp'].isna().sum()
    if n_unmatched:
        print(f"WARNING: {n_unmatched} proteins had no matching isolate in metadata")

    return data_df_merged, num_features, feature_names

def calculate_shap_values(clf, X_test, X_train):
    # SHAP values
    
    # Create the SHAP explainer for the Random Forest model
    explainer = shap.TreeExplainer(clf)

    # Calculate SHAP values for the test set
    shap_values = explainer.shap_values(X_test)

    def shap_values_to_list(shap_values, clf):
        shap_as_list=[]
        for i in range(len(clf.classes_)):
            shap_as_list.append(shap_values[:,:,i])
        return shap_as_list

    shap_as_list = shap_values_to_list(shap_values, clf)
    
    feature_names = X_train.columns
    rf_resultX = pd.DataFrame(shap_as_list[0], columns = feature_names)

    vals = np.abs(rf_resultX.values).mean(0)

    shap_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                    columns=['col_name','feature_importance_vals'])
    shap_importance = shap_importance.sort_values(by=['feature_importance_vals'],
                                ascending=False).set_index('col_name')

    return shap_importance

def run_the_classifier(data_df_merged_balanced, num_features, feature_names):

    X = data_df_merged_balanced.iloc[:, 0:num_features]  # Features
    y = data_df_merged_balanced['Comp']  # Labels
    y = pd.get_dummies(y, dtype=int)

    # Split dataset into training set and test set for accuracy calc
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, train_size=0.8, test_size=0.2, stratify=y
    )

    # --- VARIANCE FILTER: drop near-constant alignment-position columns ---
    # fit on X_train only to avoid leaking test-set info into feature selection
    selector = VarianceThreshold(0.9 * (1 - 0.9))
    X_train_filtered = selector.fit_transform(X_train)
    X_test_filtered = selector.transform(X_test)

    kept_mask = selector.get_support()
    kept_feature_names = feature_names[kept_mask]

    print(f"Variance filter: kept {kept_mask.sum()} of {len(feature_names)} features")

    # rewrap as DataFrames so downstream code (feature importances, column names) still works
    X_train = pd.DataFrame(X_train_filtered, columns=kept_feature_names, index=X_train.index)
    X_test = pd.DataFrame(X_test_filtered, columns=kept_feature_names, index=X_test.index)

    # get the best parameters
    param_grid = {
        'max_depth': [5, 10, None],
        'max_features': ['sqrt', 'log2'],
        'min_samples_leaf': [1, 2, 4],
        'min_samples_split': [2, 5, 10],
        'n_estimators': [1000, 2000, 5000],
        'class_weight': ['balanced', None]
    }

    random_search = RandomizedSearchCV(
        RandomForestClassifier(n_jobs=-1),
        param_grid,
        n_jobs=-1
    )

    random_search.fit(X_train, y_train)
    best_params = random_search.best_estimator_
    print(f"Hyperparamter tuning... best parameters found: {best_params}")

    # --- TRAIN THE MODEL FOR ACCURACY CALCULATION AND FEATURES ---
    clf = best_params
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)

    accuracy = metrics.accuracy_score(y_test, y_pred)
    report = metrics.classification_report(y_test, y_pred, output_dict=True)
    conf_matrix = confusion_matrix(y_test.values.argmax(axis=1), y_pred.argmax(axis=1))
    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
    print("Classification report:", metrics.classification_report(y_test, y_pred))
    print("confusion matrix: \n", conf_matrix)

    # NOTE: use kept_feature_names, not feature_names — clf was trained only on
    # the columns that survived the variance filter, so importances_ is that length
    feature_imp = pd.Series(clf.feature_importances_, index=kept_feature_names).sort_values(ascending=False)
    feature_imp_df = feature_imp.to_frame()

    if num_features <= 5000:  # only calc shap if #features is less than 5000, otherwise takes too long.
        shap_imp_df = calculate_shap_values(clf, X_test, X_train)

    else:
        print(f"Number of features ({num_features}) > 5000, shap values not calculated")
        shap_imp_df = pd.DataFrame(0, index=feature_imp_df.index, columns=["shap_value"])

    return kept_feature_names, clf, accuracy, report, conf_matrix, feature_imp_df, shap_imp_df

def create_output_files(data_df_merged_balanced, fasta_path, output_dir, best_params,accuracy, report, conf_matrix,gini_imp, shap_imp_df):
    # --- CREATE AN OUTPUT FILE FOR THE MODEL REPORT ---
    report_dict={}
    for i in range(0, 2):
        clust_num = fasta_path.stem[0:-8]
        precision = report[f"{i}"]['precision']
        recall = report[f"{i}"]['recall']
        f1 = report[f"{i}"]['f1-score']
        support = report[f"{i}"]['support']
        report_dict[i] = [clust_num,best_params,accuracy, precision, recall, f1, support, conf_matrix]

    report_dict

    features_out_path = output_dir.joinpath(f"{clust_num}_Feature_importance.tsv")
    report_out_path = output_dir.joinpath(f"{clust_num}_RandomForest_Report.tsv")
    data_df_out_path = output_dir.joinpath(f"{clust_num}_data_df.tsv")

    # --- CREATE AN OUTPUT FILE FOR RF MODEL REPORT ---
    report_df = pd.DataFrame.from_dict(report_dict, orient='index', columns=['clust_num','RF_parameters','accuracy', 'precision', 'recall', 'f1', 'support', 'conf_matrix'])
    report_df
    report_df.to_csv(report_out_path, sep='\t', mode='a', index=False, 
                    header=not report_out_path.exists())
    
    # --- CREATE AN OUTPUT FILE FOR THE IMPORTANT FEATURES ---
    features_df = shap_imp_df.join(gini_imp).reset_index()
    # print(features_df)

    features_df.rename(columns={"index": "feature", "shap_value": "shap_imp"}, inplace=True)
    features_df['cluster_num'] = fasta_path.stem[0:-8]
    print(features_df)

    features_df.to_csv(features_out_path, sep='\t', mode='a', index=False, 
                    header=not features_out_path.exists())
    
    # --- CREATE AN OUTPUT FILE FOR DATA DF ---
    data_df_merged_balanced.to_csv(data_df_out_path, sep='\t')

    return

if __name__ == '__main__':
    fasta_path =  Path(sys.argv[1])
    metadata_path = Path(sys.argv[2])
    output_dir = Path(f"RF_results")
    
    Path(f"RF_results").mkdir(parents=True, exist_ok=True)

    data_df_merged_balanced, num_features,feature_names = parse_seqs_and_metadata(fasta_path, metadata_path)

    kept_feature_names, best_params,accuracy, report, conf_matrix, gini_imp, shap_imp_df = run_the_classifier(data_df_merged_balanced, num_features,feature_names)

    create_output_files(data_df_merged_balanced, fasta_path, output_dir, best_params,accuracy, report, conf_matrix, gini_imp, shap_imp_df)

    print(f"{fasta_path.stem[0:-8]} finished... Model accuracy: {accuracy}")