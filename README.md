# Instructions

This code allows one to work with local views and slack in the minimisation program discussed in the paper at [https://doi.org/10.37236/7743](https://doi.org/10.37236/7743). 

# Prerequisites

- Sage (http://www.sagemath.org/), tested with v8.0.
- Jupyter, included in Sage.
- Additional Python packages: `cytoolz`, `joblib`, and `tqdm`. These are best installed from the command line with `sage -pip install cytoolz joblib tqdm`.

# Use

The code for computing the necessary representative local views and slack is in the files `local_views.pyx` and `properties.pyx`. The Jupyter notebook `d4.ipynb` allows the user to step through the computations that produce the desired results. 

Open Sage with the Jupyter notebook (e.g. `sage -n jupyter` on the command line) and open `d4.ipynb` in the Jupyter interface. Then step through the cells, executing the commands in order.

# Description of included ancillary files

- `d4.ipynb`, a Jupyter notebook for performing computations
- `local_views.pyx`, Cython code for computing sets of representative local views 
- `properties.pyx`, Cython code for computing properties of local views
- `README.md`, this readme file with instructions for using the code

Copyright (C) 2018 Ewan Davies, email: research@ewandavies.org, url: www.ewandavies.org

