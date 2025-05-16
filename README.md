# Simulation and Performance Analysis of Ramp Metering Control Strategies in Traffic Networks

Musel Tabares Pardo

Filippo La Fauci

Daniel Limmy

## Overview

This project focuses on modeling and comparing traffic flow on a freeway divided into multiple sections. It examines the dynamics of traffic density, flow, and ramp metering by considering both open-loop (no control) and closed-loop (ALINEA control) strategies for managing on-ramp traffic. The primary goal is to assess how different control approaches impact overall traffic conditions, particularly at merging points where vehicles enter the freeway.

## Objectives

- **Modeling Traffic Dynamics:**  
  Develop a discretized mathematical model of a freeway segmented into \( N \) sections. Each section is characterized by parameters such as length, maximum density, free-flow speed, and congested wave speed.

- **On-Ramp Demand and Control:**  
  - Simulate various on-ramp demand profiles (e.g., baseline, high constant, low constant, oscillatory, step change, and random fluctuations).  
  - Implement ALINEA, a closed-loop ramp metering strategy, to adjust on-ramp inflows based on real-time measurements of mainline traffic density.

- **Scenario Analysis:**  
  Evaluate the performance of both open-loop (no control) and closed-loop (ALINEA control) strategies across 10 distinct scenarios. These scenarios vary traffic conditions by altering demand profiles, initial densities, free-flow speeds, and maximum densities.

- **Performance Metrics:**  
  Assess key performance indicators including:
  - Downstream density,
  - Outflow (vehicles exiting the freeway),
  - Overall average density over time.  
  Visual comparisons using time series plots and density heatmaps help in understanding the impact of the control strategies.

## Methodology

1. **Traffic Model Development:**  
   - The freeway is divided into \( N \) sections.
   - A simulation is run over a time horizon \( K \) with a fixed sampling time \( T \).
   - Boundary conditions are imposed at the upstream and downstream ends.

2. **Demand Profile Generation:**  
   - A helper function generates on-ramp demand matrices based on various predefined scenarios.
   - Demand profiles include constant, oscillatory, step change, and random fluctuations to represent realistic traffic demand variations.

3. **Simulation of Control Strategies:**  
   Two simulation approaches are implemented:
   - **No Control (Open-Loop):**  
     The ramp inflows follow the predetermined demand profiles without any adjustments.
   - **ALINEA Control (Closed-Loop):**  
     The ramp metering rate is adjusted in real time based on the difference between the measured mainline density and a target density using the ALINEA control law.

4. **Visualization and Analysis:**  
   For each scenario, detailed plots are generated:
   - **Ramp-by-Ramp Analysis:**  
     For every on-ramp, graphs display:
     - Ramp demand vs. time,
     - Downstream density vs. time,
     - Downstream flow vs. time for both control strategies.
   - **Overall Performance:**  
     A comparative plot of the average density over the entire freeway for both strategies.
   - **Density Heatmaps:**  
     Heatmaps illustrate the spatial and temporal evolution of traffic density across all freeway sections for both strategies, using a common color scale for direct comparison.

## Expected Outcomes

- **Performance Insights:**  
  Identify conditions under which ALINEA control significantly improves freeway performance compared to no control.
  
- **Sensitivity Analysis:**  
  Gain insights into how different on-ramp demand scenarios and traffic parameters affect overall traffic flow.

- **Visual Comparisons:**  
  Produce visual and quantitative comparisons that can help inform traffic management strategies and further research in ramp metering contro
  

## Prerequisites

- **Python 3.7+**
- [**Poetry**](https://python-poetry.org/) installed globally if used for environment
- [**Miniconda / Anaconda**](https://docs.conda.io/en/latest/miniconda.html) if using Conda environment

---

## Installation

### Option 1: Using Poetry

1. **Clone the Repository**

   ```bash
   git clone https://github.com/musel25/actm_control.git
   cd actm_control
   ```

2. **Install Poetry dependencies**

   If you want the virtual environment to live inside the project folder:

   ```bash
   poetry config virtualenvs.in-project true
   ```

   Then install dependencies:

   ```bash
   poetry install
   ```

---

### Option 2: Using Conda

1. **Clone the Repository**

   ```bash
   git clone https://github.com/musel25/actm_control.git
   cd actm_control
   ```

2. **Create and activate the Conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate actm_control
   ```

3. **Use Poetry with Conda Python (optional)**

   Inside the activated Conda environment:

   ```bash
   poetry config virtualenvs.create false
   poetry install
   ```

---

## Running the Project

### From Jupyter Notebook (Browser):

```bash
poetry run jupyter notebook actm_simulations.ipynb
```

### From VS Code:

- Press `Ctrl+Shift+P` → select **Python: Select Interpreter**
- Choose the Python interpreter:

   - If using Poetry venv in project folder:
     ```bash
     .venv/bin/python
     ```
   - If using Conda:
     ```bash
     ~/miniconda3/envs/actm_control/bin/python
     ```

- Open and run `actm_simulations.ipynb`.

---
