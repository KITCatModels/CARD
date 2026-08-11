# CARD – Computer Aided Reactor Design
[![DOI](https://zenodo.org/badge/1330736885.svg)](https://doi.org/10.5281/zenodo.21888110)

A collection of Jupyter notebooks for teaching **Computer Aided Reactor Design (CARD)** to chemical engineering students. The material accompanies a series of training sessions that introduce students to Python-based modeling and simulation of chemical reactors, using open-source scientific software.

## Overview

The training sessions reinforce previously introduced problem-solving strategies through hands-on programming exercises in Python, using open-source software. Each session begins with a brief presentation of the mathematical formulation of the problem, followed by a live demonstration of its implementation in a Jupyter notebook based on a simplified example. Students are then assigned a slightly more complex programming task, completed either during a tutor-guided session or independently at home; solutions are presented and discussed at the beginning of the following session.

Initially, students are introduced to the fundamental concepts of Python programming, including variables, loops, conditional statements, and functions, which serve as the foundation for subsequent implementations. Next, commonly used libraries for solving typical engineering problems are introduced, including NumPy for numerical computations, SciPy and CasADi for root finding and numerical integration, and Matplotlib for data visualization. Building on these concepts, a typical modeling workflow is presented, consisting of four steps: mathematical problem formulation, numerical solution formulation, code implementation, and data analysis.

The complexity of the sessions increases progressively as simplifying assumptions in the models are gradually removed. Students first solve a steady-state CSTR model using a root-finding algorithm, which is then extended to a transient CSTR, resulting in a system of ODEs solved with a numerical integrator. Next, a steady-state PFR model is formulated as an ODE by transforming the spatial coordinate into residence time, highlighting the analogy between the transient CSTR and the steady-state PFR. Building on this foundation, students learn to transform PDEs into systems of DAEs using the finite volume method (FVM), chosen for its conservative formulation. Corresponding examples include species transport in a PFR and diffusion-reaction within catalyst pellets, resolved in both spatial and temporal domains. Finally, a polytropic PFR model integrates all previously introduced concepts such as species transport, energy transport, and intraparticle diffusion via an effectiveness-factor approach, while maintaining mass conservation. As a final step, this comprehensive model is reorganized into a system of classes to demonstrate the principles of object-oriented programming, illustrating how complex engineering models can be structured into modular, debuggable components suitable for larger software projects.

## Training Sessions

| # | Notebook | Topic |
|---|----------|-------|
| 1 | [`training_session_1.ipynb`](training_session_1.ipynb) | Introduction to Python (variables, arrays, loops, conditionals, functions, root finding, plotting) |
| 2 | [`training_session_2.ipynb`](training_session_2.ipynb) | Modeling a CSTR and a CSTR cascade |
| 3 | [`training_session_3.ipynb`](training_session_3.ipynb) | Solving DAE systems: a transient CSTR and a steady-state PFR |
| 4 | [`training_session_4.ipynb`](training_session_4.ipynb) | Spatial discretization: modeling a transient PFR |
| 5 | [`training_session_5.ipynb`](training_session_5.ipynb) | Species transport inside a catalyst pellet |
| 6 | [`training_session_6.ipynb`](training_session_6.ipynb) | Modeling a transient polytropic PFR (species, energy, and mass conservation) |
| 7 | [`training_session_7.ipynb`](training_session_7.ipynb) | Modeling a transient polytropic catalytic fixed bed reactor |
| 8 | [`training_session_8.ipynb`](training_session_8.ipynb) | Model implementation with object-oriented programming |

## Repository structure

```
.
├── training_session_1.ipynb    # Session notebooks (1-8)
├── ...
├── src/                        # Reusable Python modules (e.g. CatalystBed, Parameters, Properties)
├── figures/                    # Figures to illustrate example problems
├── pyproject.toml              # Project metadata and dependencies
└── uv.lock                     # Locked dependency versions
```

## Getting started

This project uses [uv](https://docs.astral.sh/uv/) for dependency management and requires Python 3.13 or later.

```bash
# clone the repository
git clone https://github.com/KITCatModels/CARD.git
cd card

# install dependencies into a local virtual environment
uv sync

# launch Jupyter
uv run jupyter lab
```

Alternatively, install with `pip`:

```bash
pip install -e .
jupyter lab
```

### Dependencies

- [CasADi](https://web.casadi.org/) – algebraic modeling and numerical optimization
- [SciPy](https://scipy.org/) – root finding and numerical integration
- [Matplotlib](https://matplotlib.org/) – data visualization
- [ipykernel](https://github.com/ipython/ipykernel) / [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) – Jupyter integration

See [`pyproject.toml`](pyproject.toml) for exact version constraints.

## How to cite

If you use this material, please cite the accompanying paper:

```bibtex
@article{TODO_citekey,
  author  = {TODO},
  title   = {TODO},
  journal = {TODO},
  year    = {TODO},
  doi     = {TODO}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
