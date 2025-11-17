### **Portfolio Assignment: Lung Cancer Classification with Dask and XGBoost**

---

#### **Objective**:
This assignment is designed to provide students with hands-on experience in handling large-scale datasets using Dask and training a machine learning model with XGBoost. 

---

#### **Learning Goals**:
1. Understand the basics of parallel computing and scalability with Dask.
2. Use Dask DataFrame and Delayed for data manipulation and parallelized workflows.
3. Integrate Dask with XGBoost for distributed machine learning.
4. Compare the performance of Dask workflows with traditional pandas/numpy workflows.

---

### **Scenario**
You are tasked with building a scalable machine learning pipeline. You will process a large dataset using Dask, train a predictive model with XGBoost, and evaluate the workflow's computational efficiency.

---

### **Dataset Description**
You can use your own project or the data provided by the teacher. 

---

### **Tasks**

#### **1. Getting Started with Dask**
- Load the data using **Dask DataFrames**.
- Compare the memory usage and loading times with pandas DataFrames.

---

#### **2. Initial Data Preprocessing**
Write **Dask Delayed** functions to:
- Handle missing values in the data 
- Optinal: Filter low-quality data in a gene expression dataset 

---

#### **3. Exploratory Data Analysis (EDA)**
Write **Dask Delayed** functions to:
- Calculate summary statistics (e.g., mean, median, standard deviation) for clinical variables and gene expression data.
- Relevant metrics such as top 10 most variable genes across patients.
- Group patients and make useful overviews
- Visualize the distribution of the target variable.
- Compare processing time and outputs with pandas workflows.

---

#### **4. Further Preprocessing**
- Normalize numerical features using z-score scaling with Dask.
- Encode the target variable 
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

