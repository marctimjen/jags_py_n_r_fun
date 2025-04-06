# Initialize NumerAPI - the official Python API client for Numerai
from numerapi import NumerAPI
napi = NumerAPI()

# list the datasets and available versions
all_datasets = napi.list_datasets()
dataset_versions = list(set(d.split('/')[0] for d in all_datasets))
print("Available versions:\n", dataset_versions)

# Set data version to one of the latest datasets
DATA_VERSION = "v5.0"

# Print all files available for download for our version
current_version_files = [f for f in all_datasets if f.startswith(DATA_VERSION)]
print("Available", DATA_VERSION, "files:\n", current_version_files)


import json

# download the feature metadata file
napi.download_dataset(f"{DATA_VERSION}/features.json")

# read the metadata and display
feature_metadata = json.load(open(f"{DATA_VERSION}/features.json"))
for metadata in feature_metadata:
  print(metadata, len(feature_metadata[metadata]))






# training data contains features and targets
# training_data = pd.read_parquet("/workspaces/jags_py_n_r_fun/num/data/num.parquet")
# 
# # Take a subset of 1000 rows from training_data
# subset_data = training_data.loc["n1a9a64f5c6228cd"]
# print(subset_data.shape)  # Print the shape of the subset to verify
# # Save the subset to a new file
# subset_data.to_parquet("/workspaces/jags_py_n_r_fun/num/data/num_subset.parquet")


# subset_data = pd.read_parquet("/workspaces/jags_py_n_r_fun/num/data/num_subset.parquet")
# 
# 
# for i in subset_data.columns:
#     print(i)
# 
# 
# print()  # Added to verify the data loaded correctly
# 
# row = subset_data.loc["n1a9a64f5c6228cd"]
# 
# print(row)  # Print the extracted row to verify