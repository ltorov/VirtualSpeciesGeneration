# Virtual Species Generation
Research practice by Luisa Toro Villegas
Advisor: Daniel Rojas Díaz

- Part of [Non Parametric Density Estimation For Niche Modelling](https://github.com/coberndorm/Niche-Modelling) by Camilo Oberndorfer Mejia and Miguel Valencia Ochoa.

## Table of Contents
* [General Info](#general-information)
* [Data](#data)
* [Features](#features)
* [Dependencies](#dependencies)
* [Setup](#setup)
* [Project Status](#project-status)
* [Room for Improvement](#room-for-improvement)
* [Acknowledgements](#acknowledgements)
<!-- * [License](#license) -->


## General Information
This project aims to generate virtual species as a crucial complement to species distribution models (SDMs). SDMs often require validation with real-life species data, which can be challenging and biased. Virtual species are created by defining their niche as a function of environmental variables and simulating species occurrences on a given map. This project is written entirely in Matlab and aims to explore the variable space and generate variability in the simulated species to test the accuracy of SDMs.


## Data

The environmental variables data in this project should be in .asc file format and contain the environmental variables in a given area. The data should be stored in the data folder in the following format:

- Each environmental variable should be stored in a separate .asc file.
- The .asc files should contain the values of the variable and the coordenates.

An example of the environmental variables data files is available in the data folder.


## Features
List the ready features here:
- Coefficient Method
- Beta Distributions Method
- Harmonic Functions Method


## Dependencies

This project requires the following dependencies:

- MATLAB (version R2018a or higher)
- .asc files containing the environmental variables in a given area

## Setup

To get started with this project, you can clone the repository using the following command:

`git clone https://github.com/ltorov/VirtualSpeciesGeneration.git`

To use this code, open the main.mlx file in MATLAB and follow the instructions provided.


## Project Status
Project is: _in progress_.


## Room for Improvement
Include areas you believe need improvement / could be improved. Also add TODOs for future development.

Room for improvement:
- Tuning hyperparameters.

To do:
- Further experimentation.
- Compare against similar methods in the literature.
- Use this to test the Frontier Depth method against MaxEnt and MaxLike.


## Acknowledgements

- This project was inspired by Daniel Rojas Diaz.
