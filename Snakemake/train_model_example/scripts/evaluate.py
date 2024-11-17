import pandas as pd
from sklearn.metrics import accuracy_score
import joblib

# Load preprocessed data
data = pd.read_csv(snakemake.input.data)
X = data.drop("target", axis=1)
y = data["target"]

# Load the trained model
model = joblib.load(snakemake.input.model)

# Make predictions and calculate accuracy
predictions = model.predict(X)
accuracy = accuracy_score(y, predictions)

# Save the results
# use output[0] since it uses a temp file
with open(snakemake.output[0], "w") as f:
    f.write(f"Accuracy: {accuracy:.2f}\n")
