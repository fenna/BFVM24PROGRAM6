import pandas as pd
from sklearn import datasets
from sklearn.model_selection import train_test_split

# Load the Iris dataset
iris = datasets.load_iris()
data = pd.DataFrame(data=iris.data, columns=iris.feature_names)
data['target'] = iris.target

# Split the data 
train_data, test_data = train_test_split(data, test_size=0.3) 

# Save preprocessed data
train_data.to_csv(snakemake.output.train, index=False)
test_data.to_csv(snakemake.output.test, index=False)
