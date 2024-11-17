import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load results
results_df = pd.read_csv(snakemake.input[0])

# Create a bar plot
plt.figure(figsize=(10, 6))
sns.barplot(x='model', y='accuracy', data=results_df)
plt.title('Model Accuracy Comparison')
plt.xlabel('Model')
plt.ylabel('Accuracy')
plt.ylim(0, 1)  

# Save the plot
plt.savefig(snakemake.output[0])
