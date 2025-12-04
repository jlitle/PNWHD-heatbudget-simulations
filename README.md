# Mussel Heat Budget Simulation for the Pacific Northwest Heat Dome

This repository contains a minimal reproducible example of the heat budget simulation code used in our study on mussel heat stress during the Pacific Northwest Heat Dome. 

Here, we demonstrate the core simulation pipeline. Full datasets and analysis scripts are not included here; see the Data Availability section below for links to archived data and analysis scripts.

---

## Repository Structure

```text
repo-root/
├─ R/
│  └─ heatbudget_functions.R
├─ examples/
│  └─ run_transfer_simulations.R
├─ data/
│  └─ example_inputs.csv
└─ README.md
```

---

## Running the Example

1. Clone the repository:

git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git

2. Open `examples/run_transfer_simulations.R` in R or RStudio.

3. Run the script. Our example uses `data/example_inputs.csv` and demonstrates how the simulation functions work.  

> Note: The outputs from this minimal example do not reproduce results from the study; they illustrate how to run the heat budget simulation.

---

## Dependencies

The example simulation script requires the following R packages:

- tidyverse (includes ggplot2, dplyr, tidyr, etc.)  
- TrenchR  

Install them with:

install.packages("tidyverse")
install.packages("TrenchR")

---

## Data Availability

- This repo includes minimal example input data for the heat budget simulation: `data/example_inputs.csv`  
- We archived full datasets and analysis code supporting the paper in Dryad (currently private for peer review) : [DOI link here]  

---

## Citation

Citation information pending publication
