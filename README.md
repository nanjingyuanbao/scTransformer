# scTransformer
scTransformer, a Transformer-based model, can annotate cell types. 

# Requirement
```
R 4.2.2
Python 3.9.12
PyTorch >= 1.5.0
numpy
pandas
scipy
sklearn
Scanpy
random
```
## Input file
```
reference dataset.

cell type label of reference dataset.

query dataset.
```
## Output file
```
After training the scTransformer model, the model will be save at: "log/scTransformer.tar".
The model prediction is saved in the log/y_predicts.npy.
