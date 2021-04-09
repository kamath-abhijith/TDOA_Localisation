# TIME DIFFERENCE OF ARRIVAL (TDOA) LOCALISATION

Time difference of arrival (TDOA) is an estimation method for target localisation. Repository contains Python3 scripts for:

- Cramer-Rao Lower Bound for unbiased target estimation,
- Maximum-Likelihood estimation using gradient descent,
- Best-Linear-Unbiased estimation.

## Documentation

`docs/instructions.pdf` contains the necessary instructions for the assignment. `docs/solutions.pdf` contains theory and algorithms.

## Installation

Clone this repository and install the requirments using
```shell
git clone https://github.com/kamath-abhijith/TDOA_Localisation
conda create --name <env> --file requirements.txt
```

## Run

- `data/TDOA_data.mat` contains the data. See `data/README_hw1` for the structure.
- Run `tdoa_ml.py` to generate Figure 1(a).
- Run `tdoa_blue.py` to generate Figure 1(b).
- Run `tdoa_montecarlo.py` to generate Figure 2.