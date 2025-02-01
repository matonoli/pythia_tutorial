# pythia_tutorial

This repository provides a tutorial for using **Pythia8** with **ROOT** for event generation and analysis. Users can explore event simulation via Jupyter Notebooks, Python scripts, and C++ executables.

## Getting Started

### Run on GitHub Codespaces
Click the badge below to launch the environment instantly:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/matonoli/pythia_tutorial?quickstart=1)

### Run Locally
Clone the repository and set up the Conda environment:

```bash
git clone https://github.com/matonoli/pythia_tutorial.git
cd pythia_tutorial
conda env create -f environment.yml
conda activate pythia8-env
```

## Usage

- **Jupyter Notebooks**: Open and run the `pythia_exercise.ipynb` file and select the (base) Python 3.10.6 kernel
- **Python Scripts**: Run example scripts from the terminal:
  ```bash
  python simulate_jpsi.py
  ```
- **C++ Code**: Compile and execute C++ examples:
  ```bash
  g++ -o simulate_jpsi simulate_jpsi.cpp `root-config --cflags --libs` -I$PYTHIA8/include -L$PYTHIA8/lib -lpythia8  -Wl,-rpath,$PYTHIA8/lib
  ./simulate_jpsi
  ```

