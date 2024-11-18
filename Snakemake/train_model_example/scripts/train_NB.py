import pandas as pd
from sklearn.naive_bayes import GaussianNB
import joblib

# Load preprocessed data
data = pd.read_csv(snakemake.input[0])
X = data.drop("target", axis=1)
y = data["target"]

# Train model
model = GaussianNB()
model.fit(X, y)

# Save the trained model
joblib.dump(model, snakemake.output[0])