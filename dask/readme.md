## set up

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

The data used in the notebook `intro_dask.ipynb` is to be found at assemblix2019:/data/datasets/DS6/cluster_chem
copy the file `all_nps.txt` to <your workdirectory>/data 
```
#split bigfile into smaller files
cd data
mkdir all_nps_chunks
split -l 5000  data/all_nps.txt all_nps_chunks/part_
``
