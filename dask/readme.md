## set up

You can follow the setup in the tutorial of dask. https://docs.dask.org/en/stable/install.html.

If this does not work due to dependency issues you can use the instructions below:

```
  conda create -n dask_env python=3.10
  conda activate dask_env
  conda install -c conda-forge dask=2024.11.2 dask-ml=2024.4.4
  pip install dask-expr==1.1.19
  conda install -c conda-forge ipykernel notebook
  python -m ipykernel install --user --name=dask_env
  conda install -c conda-forge pandas
  conda install -c conda-forge numpy
  
```
## data 

The data used in the notebook `intro_dask.ipynb` is to be found at `assemblix2019:/data/datasets/DS6/cluster_chem/cleaned.zip`

The data used in the ML framework is to be found at `assemblix2019:/data/datasets/PROG6/subset_1000.csv`. This is a subset of the cell type expression data of Zheng[1]. The framework should however work with any tabular data that has features and a label.


## Example framework XGBoost

The example framework for XGBoost can be found at https://github.com/fenna/ml-framework

## References

[1] Zheng, et all (2017). Massively parallel
digital transcriptional profiling of single cells. Nature Communications, 8(1), 14049.
https://doi.org/10.1038/ncomms14049


```

