
## Parallelizing Machine Learning Workflows with Snakemake

Snakemake is a Workflow management system focused on bioinformatics and reproducible data science pipelines. It automates and paralelize complex data processing pipelines. It’s highly popular for reproducible research workflows, especially in the bioinformatics research. it uses a rule-based approach to define steps in a pipeline, automatically managing dependencies and rerunning only the necessary parts of a pipeline when files change. It’s Python-based. It can be combined with SLURM. This assignment will introduce you in the benefits of Snakemake. In this assignment you will refactor an existing machine learning script to incorporate Snakemake and leverage parallelization efficiently. Additionally, it offers a foundation for documenting your experiments in a more structured and organized manner. 

You will choose at least one of the following workflows:

  a) Dataprocessing of files (for instance images)
  b) Hyperparameter Tuning
  c) Trying Different Machine Learning Models

You can use the script of your current Omics project or a previous developed script. If you want to work on the dataprocessing case, but you do not have a use case yourself please contact the teacher. 

---

### Assignment Details:

**Preparation**

Read the [paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030138) of this [blog](https://waterdata.usgs.gov/blog/snakemake-for-ml-experiments/)

Snakemake is installed on the linux grid. Read https://snakemake.github.io/ for installation instructions if you want to install it on other systems


**Task Overview:**

  1. **Select a Workflow**: Choose one of the workflows described below.
  2. **Refactor the Script**: Take an existing (machine learning) script and refactor it to use a Snakemake workflow.
  3. **Design the Workflow**: Plan and define the Snakemake workflow, specifying rules for parallel execution, resource allocation, and dependency management.
  4. **Implement the Workflow**: Write Snakemake rules and a Snakefile to handle the parallelization of your chosen task.
  5. **Document the Workflow**: Provide clear and concise documentation, detailing your refactoring approach, workflow design decisions, and how to run the workflow. Use a dag picture to visualize the parallelization. 

---

### Workflow Options:

**1. Data Processing of Image Files:**

- Refactor a data processing script that uses multiple files to use Snakemake for parallel preprocessing tasks like resizing, normalization, and feature extraction.

  - Write Snakemake rules to process the files in parallel.
  - Efficiently handle different preprocessing steps.
  - Specify the location of input files and the output.


**2. Hyperparameter Tuning:**

- Refactor a script that performs hyperparameter tuning to run multiple configurations in parallel using Snakemake.

  - Use Snakemake rules to manage parallel execution of different hyperparameter sets.
  - Choose a machine learning model.
  - Compare the performance metrics of each hyperparameter configuration.
  - Use the outcome to document experiments systematically. 
  - Visualize the results for comparison.
  
**3. Trying Different Models:**

- Refactor a script that tests different machine learning models, parallelizing the model training and evaluation using Snakemake.

  - Implement Snakemake rules to train and compare various models.
  - Collect and analyze performance metrics of each model.
  - Use the outcome to document experiments systematically. 
  - Visualize the results for comparison.

---

### Assignment Requirements:

1. **Refactoring**: Take your existing (machine learning) script and modify it to work within a Snakemake framework.
2. **Snakefile**: Create a Snakefile that defines input, output, rules, and resource management.
3. **Parallelization**: Use Snakemake’s parallelization features to optimize your workflow.
4. **Resource Management**: To prevent system overload, ensure efficient resource allocation. While Assemblix has 80 cores available, you are not permitted to monopolize all of them. Remember, you are not alone—share resources considerately with others!
5. **Error Handling**: Implement error handling and logging to manage task failures.
6. **Reproducibility**: Ensure the workflow is reproducible and easy to run across different environments.

---

### Submission:

- Submit your refactored script(s), the Snakefile, and any supporting files.
- Provide a brief report (1-2 pages) explaining your refactoring process, design choices, and a summary of the results.
- Make sure your code is well-documented, well-organized with comments where appropriate.
- Add a readme with clear instructions.

**Bonus Task**: Extend your workflow to run on a slurm cluster, utilizing Snakemake's cluster features.

---


### Assessment Criteria:

1. **Use of Relevant Literature and Sources (10%)**:  

   - References are cited consistently and scientifically, with a focus on literature that supports workflow design, Snakemake usage, and parallel computing.
   - Demonstrates the ability to research and incorporate insights from high-quality academic sources and relevant online documentation on Snakemake and workflow optimization.

3. **Quality of Workflow and Analysis (35%)**:  

   - Constructs a well-functioning Snakemake workflow that correctly handles dependencies and efficiently parallelizes tasks.
   - Visualizations and output summaries are self-explanatory and accurate.
   - Demonstrates the effectiveness of the workflow by analyzing performance metrics and evaluating computational efficiency.

4. **Code Repository and Snakemake Implementation (30%)**:  

   - Code is clean, well-documented, and adheres to best practices for readability and maintainability.
   - The Snakemake workflow is organized logically, with a clear Snakefile, well-defined rules, and comprehensive instructions provided in a Readme file.
   - The repository includes version control and appropriate software licensing to ensure reproducibility.
   
5. **Supplementary Information and Documentation (25%)**:  

   - Includes a description of how the original script was refactored into a Snakemake workflow.
   - Provides a well-justified and clear explanation of the chosen Snakemake workflow design, including rule dependencies, resource allocation, and parallelization strategies.
   - Critically assesses the efficiency of the implemented workflow and considers improvements based on results and feedback.
   - Reflects on challenges faced during implementation and addresses potential limitations in the workflow.


---

Hints for a repository structure:
```{}
├── .gitignore
├── README.md
├── LICENSE.md
├── images
│   ├── dag.png
├── workflow
│   ├── rules
|   │   ├── module1.smk
|   │   └── module2.smk
│   ├── scripts
|   │   ├── script1.py
|   │   └── script2.R
|   └── main.smk
├── config
│   ├── config.yaml
├── results
└── resources
```

References: 

- https://github.com/fenna/BFVH4DSP1
- https://snakemake.github.io/


J. M. Sadler, A. P. Appling, J. S. Read, S. K. Oliver, X. Jia, J. A. Zwart, V. Kumar
Multi-Task Deep Learning of Daily Streamflow and Water Temperature, Water resources research, Volume58, Issue4 https://doi.org/10.1029/2021WR030138


