### **Portfolio Assignment: Lung Cancer Classification with Dask and XGBoost**

---

#### **Objective**:
This assignment is designed to provide students with hands-on experience in handling large-scale datasets using Dask and training a machine learning model with XGBoost. By the end of this project, students will classify lung cancer tumor subtypes (e.g., Squamous vs. Adenocarcinoma) and analyze the scalability of workflows using Dask.

---

#### **Learning Goals**:
1. Understand the basics of parallel computing and scalability with Dask.
2. Use Dask DataFrame and Delayed for data manipulation and parallelized workflows.
3. Integrate Dask with XGBoost for distributed machine learning.
4. Compare the performance of Dask workflows with traditional pandas/numpy workflows.

---

### **Scenario**
You are tasked with building a scalable machine learning pipeline for a lung cancer research project. The goal is to classify tumor subtypes based on gene expression data and clinical features. You will process a large dataset using Dask, train a predictive model with XGBoost, and evaluate the workflow's computational efficiency.

---

### **Dataset Description**

1. **Clinical Data**:
   - **Format**: Excel file containing 89 patient samples.
   - **Content**: Patient tumor source location, tumor stage, and other clinical variables.
   - **Source**: [Lung3.metadata.xls](https://wiki.cancerimagingarchive.net/download/attachments/16056856/Lung3.metadata.xls?version=1&modificationDate=1404237338168&api=v2).

2. **Gene Expression Data**:
   - **Format**: Matrix file in CSV format.
   - **Content**: Normalized expression levels of thousands of genes.
   - **Source**: [GEO Data](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58661/matrix/).

---

### **Tasks**

#### **1. Getting Started with Dask**
- Load the clinical and gene expression data using **Dask DataFrames**.
- Compare the memory usage and loading times with pandas DataFrames.

---

#### **2. Initial Data Preprocessing**
Write **Dask Delayed** functions to:
- Handle missing values in the clinical data 
- Filter low-quality data in the gene expression dataset 
- Compare the runtime of preprocessing tasks with pandas-based workflows.

---

#### **3. Exploratory Data Analysis (EDA)**
Write **Dask Delayed** functions to:
- Calculate summary statistics (e.g., mean, median, standard deviation) for clinical variables and gene expression data.
- Identify the top 10 most variable genes across patients.
- Group patients by clinical variables and compute the average expression of selected genes.
- Visualize the distribution of the target variable (`TumorSubtype`).
- Compare processing time and outputs with pandas workflows.

---

#### **4. Further Preprocessing**
- Normalize numerical features (e.g., gene expression data) using z-score scaling with Dask.
- Encode the target variable (`TumorSubtype`) as binary (e.g., "1" for Squamous, "0" for others).
- Split the dataset into training and test sets using **Dask**.

---

#### **5. Model Training with XGBoost**
- Train an **XGBoost Classifier** using the Dask-processed training dataset.
- Use default hyperparameters for the initial model.
- Evaluate the model on the test set using metrics such as:
  - Accuracy
  - Precision, Recall, F1-Score
  - AUC-ROC
- Measure the runtime for training and inference using both Dask and pandas workflows.

---

#### **6. Performance Analysis**
- Compare the computational performance (runtime, memory usage) of:
  - Preprocessing tasks.
  - Model training and evaluation.
- Summarize the advantages and trade-offs of using Dask for this workflow.

---

> **Bonus Challenge**:
> Use **Dask Delayed** to implement a custom parallelized workflow to compute correlations between specific clinical variables and gene expression levels.

---

### **Submission**

1. **Scripts and Notebooks**:
   - Submit your final scripts with well-documented code.
   - Ensure that the code is modular, organized, and easy to follow.
   - Provide notebooks as research logs (optional)

2. **ReadMe File**:
   - Include a `ReadMe` with clear instructions on running your code, dependencies, and dataset paths.

3. **Report**:
   - Provide a 1-2 page report summarizing:
     - Design choices.
     - Critical comparison of Dask vs. pandas workflows.
     - Challenges faced and solutions implemented.

---

### **Assessment Criteria**

1. **Use of Relevant Literature and Sources (10%)**:
   - Cites scientific literature or online resources related to Dask, XGBoost, and parallel computing.

2. **Workflow Design and Implementation (35%)**:
   - Constructs a functioning and efficient workflow.
   - Accurately performs preprocessing, EDA, and machine learning tasks using dask.

3. **Code Quality and Repository (30%)**:
   - Code is clean, well-documented, and adheres to best practices.
   - Includes a well-organized repository with version control and licensing.

4. **Supplementary Documentation and Analysis (25%)**:
   - Provides a critical reflection on workflow performance.
   - Demonstrates insight into challenges and possible improvements.

---

### **Prerequisites**
- Proficiency in Python programming.
- Familiarity with pandas and numpy.
- Basic knowledge of machine learning concepts (e.g., classification metrics).
- Understand why XGBoost is suitable for this type of data

---

