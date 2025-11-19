
  

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmartfrut.github.io/HyperFEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmartfrut.github.io/HyperFEM.jl/dev/)
[![Build Status](https://github.com/MultiSimOLab/HyperFEM/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MultiSimOLab/HyperFEM/actions/workflows/ci.yml?branch=main)
[![Coverage](https://codecov.io/gh/jmartfrut/HyperFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jmartfrut/HyperFEM.jl)

# Julia Framework for Coupled Magneto-Mechanical Simulation




<div align="justify" style="margin-left: 40px; margin-right: 40px;">

**Overview**

This repository provides a general-purpose Julia framework for simulating fully coupled magneto-mechanical behavior in soft and hard magnetic materials. The code supports both 2D and 3D simulations and includes an extensive library of constitutive models, enabling users to explore a wide range of material systems and actuation mechanisms.
While the framework is broadly applicable, several examples focus on ultra-soft magnetoactive materials relevant to bioengineering and soft robotics, where accurate modeling of magneto-mechanical coupling remains a key challenge. These examples illustrate how the toolkit can be used for realistic simulation, analysis, and design of next-generation magnetoactive structures.

**Implementation**

All numerical experiments were conducted in the Julia programming language, utilizing the [Gridap](https://gridap.github.io/Gridap.jl/stable/) framework for the numerical approximation of partial differential equations (PDEs), and the [HyperFEM](https://github.com/MultiSimOLab/HyperFEM) package for modeling nonlinear physical phenomena in multifunctional materials.
Key features of the framework include:

- A modular architecture for constitutive modeling, numerical solvers, and data analysis
- Interfaces for coupled-field simulations, enabling the interaction of magnetic and mechanical effects.
- Ready-to-use scripts for benchmark problems and topology/material optimization workflows.

</div>

## Applications
- Design of soft robotic actuators with tunable magneto-mechanical responses.
- Modeling of biocompatible magneto-responsive scaffolds.
- Exploration of anisotropic behaviors in ultra-soft magnetorheological elastomers (hMREs).
- Integration into topology optimization frameworks for intelligent material design.
 
 
## Benchmarks Gallery

 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/MultiSimOLab/Adv_Mater_Ultra_Soft_Magnet/blob/main/doc/img/ex.png" width="100%">
</p>

 

## How to cite

In order to give credit to the contributors, we ask that you please reference the paper:

Carlos Perez-Garcia, Rogelio Ortigosa, Jesus Martinez-Frutos Daniel Garcia-Gonzalez. **Topology and material optimization in ultra-soft magnetoactive structures: making advantage of residual anisotropies**. Advanced Materials. 2025. DOI: 10.1002/adma.202518489

# Project funded by:
 
- European Research Council (ERC) under the European Union’s Horizon 2020 Research and Innovation Program (grant agreement no. 947723,project: 4D-BIOMAP, and grant agreement no. 101247449, project: MAGMATED), 
- Catedra UC3M-NAVANTIA-MONODON. 
- Grant PID2022-141957OA-C22 funded by MICIU/AEI/10.13039/501100011033 and by ”ERDF A way of making Europe”.
- Grant 21996/PI/22. Programme for the development of scientific and technical research by competitive groups, included in the Regional Program for the Promotion of Scientific and Technical Research of Fundación Séneca-Agencia de Ciencia y Tecnología de la Región de Murcia.


 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/MultiSimOLab/Adv_Mater_Ultra_Soft_Magnet/blob/main/doc/img/logo.png" width="100%">
</p>
 
#  Contact

Contact the project administrators [Jesús Martínez-Frutos](jesus.martinez@upct.es), [Daniel Garcia-Gonzalez](danigarc@ing.uc3m.es) for further questions about licenses and terms of use.