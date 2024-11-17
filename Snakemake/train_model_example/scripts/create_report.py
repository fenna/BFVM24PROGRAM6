import pandas as pd

# Collect results from each model
results_files = snakemake.input
all_results = []

for file in results_files:
    with open(file, 'r') as f:
        model_name = file.split('/')[-1].split('_results')[0]
        accuracy = float(f.read().split(":")[-1].strip())
        all_results.append({'model': model_name, 'accuracy': accuracy})

# Create a DataFrame and save to CSV
results_df = pd.DataFrame(all_results)
results_df.to_csv(snakemake.output[0], index=False)
