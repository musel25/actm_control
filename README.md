# Simulation and Performance Analysis of Ramp Metering Control Strategies in Traffic Networks

Musel Tabares Pardo

Filippo La Fauci

Daniel Limmy Akejelu

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

## Results and Scenario-Specific Discussion

Below is a scenario-by-scenario analysis, outlining the behavior of both No Control and ALINEA strategies for each case:

### 1. Baseline
- **No Control:** High ramp inflows at the start cause a rapid increase in downstream density, leading to congestion and a drop in outflow. The average density remains elevated for a long period, and congestion spreads through the network.
- **ALINEA:** Ramp inflows are regulated, keeping downstream density near the target. Congestion is prevented or quickly resolved, resulting in more stable outflow and lower average density.

### 2. High Constant Demand
- **No Control:** Persistent high demand overwhelms the freeway, causing severe and sustained congestion. Downstream flow drops sharply after an initial peak.
- **ALINEA:** The controller limits ramp inflow to avoid exceeding the target density. While some congestion may occur due to high demand, it is less severe and dissipates faster than in the uncontrolled case.

### 3. Low Constant Demand
- **No Control:** Demand is low enough that congestion does not develop. Both strategies yield similar results, with stable densities and flows.
- **ALINEA:** The controller has little effect, as the system is already operating under optimal conditions.

### 4. Oscillatory Demand
- **No Control:** Density and flow oscillate in response to the demand, with possible congestion during high-demand phases. The system struggles to recover after each surge.
- **ALINEA:** The controller smooths out the oscillations, preventing large spikes in density and maintaining more consistent flow.

### 5. Additional On-ramp
- **No Control:** The extra on-ramp increases total inflow, making the system more prone to congestion. Downstream densities rise quickly, and congestion is widespread.
- **ALINEA:** The controller manages all ramps, distributing inflow to keep densities near the target. Congestion is contained and less severe.

### 6. Lower Free-Flow Speed
- **No Control:** Reduced free-flow speed lowers the system's capacity, so congestion develops more easily and persists longer.
- **ALINEA:** The controller compensates by further limiting ramp inflow, helping to keep densities in check and improving recovery.

### 7. Lower Maximum Density
- **No Control:** The freeway saturates at a lower density, so congestion occurs sooner and is more severe.
- **ALINEA:** The controller adapts to the new capacity, keeping densities below the lower threshold and reducing congestion duration.

### 8. Step Change in Demand
- **No Control:** A sudden increase in demand causes a sharp spike in density and a rapid onset of congestion, with slow recovery.
- **ALINEA:** The controller quickly responds to the step change, limiting inflow and preventing a major breakdown.

### 9. Random Demand Fluctuations
- **No Control:** Random surges in demand lead to unpredictable spikes in density and flow, with frequent congestion events.
- **ALINEA:** The controller dampens the effect of random fluctuations, maintaining more stable traffic conditions and reducing the frequency and severity of congestion.

### 10. Higher Initial Density
- **No Control:** Starting with a higher density means the system is already close to congestion. Even moderate inflows can trigger breakdowns.
- **ALINEA:** The controller is more aggressive in limiting ramp inflow, helping the system recover to optimal density more quickly.

---

This scenario-specific analysis demonstrates that ALINEA consistently improves traffic performance, especially under challenging or variable conditions. The benefits are most pronounced in scenarios with high, variable, or sudden demand, or when the freeway's capacity is reduced.
