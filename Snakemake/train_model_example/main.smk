# Snakefile

# Configurations
models = ["random_forest", "svm"]
output_dir = "output/" 


rule all:
    input:
        output_dir + "accuracy_plot.png"

# Rule to preprocess the Iris dataset
rule preprocess_data:
    output:
        processed_data="data/processed_data.csv"
    script:
        "scripts/preprocess.py"

# Rule to train each model
rule train_model:
    input:
        data="data/processed_data.csv"
    output:
        model= output_dir + "models/{model}_model.pkl"
    params:
        model_name=lambda wildcards: wildcards.model
    benchmark: 
        output_dir + "benchmarks/train_{model}.txt"
    script:
        "scripts/train_{wildcards.model}.py"

# Rule to evaluate each model
rule evaluate_model:
    input:
        model=output_dir + "models/{model}_model.pkl",
        data="data/processed_data.csv"
    output:
        temp(output_dir + "{model}_results.txt")
    params:
        model_name=lambda wildcards: wildcards.model
    script:
        "scripts/evaluate.py"


# Rule to create a final report from the results
rule create_report:
    input:
        expand(output_dir + "{model}_results.txt", model=models)
    output:
        output_dir + "final_report.csv"
    script:
        "scripts/create_report.py"

# Rule to plot the results
rule plot_results:
    input:
        output_dir + "final_report.csv"
    output:
        output_dir + "accuracy_plot.png"
    script:
        "scripts/plot_results.py"

