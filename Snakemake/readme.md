
## Project Structure

This repository contains the code and resources an example project with snakemake. The project focuses on utilizing machine learning models to preprocess and analyze the Iris dataset, evaluate model performance, and visualize the results.

![DAG](train_model_example/dag.png)

```plaintext
.
├── data
│   └── processed_data.csv
├── output
│   ├── accuracy_plot.png
│   ├── final_report.csv
│   ├── benchmarks
│   │   ├── train_random_forest.txt
│   │   ├── train_svm.txt
│   └── models
│       ├── random_forest_model.pkl
│       └── svm_model.pkl
├── scripts
│   ├── preprocess.py
│   ├── train_random_forest.py
│   ├── train_svm.py
│   ├── evaluate.py
│   └── plot_results.py
├── main.smk
```

## Installation

To run this project, you need to have Snakemake installed. You can install snakemake using the following command:

```sh
pip install snakemake 
```

## Usage

1. **Preprocess the Data**: Run the Snakemake workflow to preprocess the Iris dataset.
    ```sh
    snakemake -s main.smk preprocess_data
    ```

2. **Train Models**: Train the specified machine learning models.
    ```sh
    snakemake -s main.smk train_model
    ```

3. **Evaluate Models**: Evaluate the performance of the trained models.
    ```sh
    snakemake -s main.smk evaluate_model
    ```

4. **Create Report**: Generate a final report combining the evaluation results.
    ```sh
    snakemake -s main.smk create_report
    ```

5. **Plot Results**: Create a bar chart comparing the model accuracies.
    ```sh
    snakemake -s main.smk plot_results
    ```

6. **Generate DAG**: Visualize the workflow DAG.
    ```sh
    snakemake -s main.smk --dag | dot -Tpdf > dag.pdf
    
    ```

## Benchmarking

The project includes benchmarking for each model training. The benchmarking results are saved in the `output/benchmarks` directory.

## License

This project is licensed under the GNU License. See the `LICENSE` file for more details.

## Contact

f.feenstra@pl.hanze.nl

