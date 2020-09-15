# COVIDsaturates
> Project description: code to accompany 'The importance of saturating density dependence for predicting SARS-CoV-2 resurgence' authored by Emily Nightingale, Oliver Brady & Laith Yakob

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)
* [Features](#features)
* [Status](#status)
* [Contact](#contact)

## General info
Contact rates between people within a population increase with human density, up until a saturating level. 
Using COVID-19 associated mortality data from England, we inferred the shape of this saturating contact-rate curve.
We then included this saturating contact rate in a mathematical model of SARS-CoV-2 transmission, and compared projections with standard density-independent and density-dependent models.

## Technologies
* Python - version 3.8
* R - version 4.0.2

## Setup
For math model, following modules must be installed: matplotlib.pyplot ; numpy ; datetime ; lmfit
For statistical analysis, following environmental requirements: 

## Features
* Standard frequency-dependent COVID-19 model: 'covid_mathmodel_fd-Copy1.ipynb'
* Standard density-dependent (LINEAR density dependence) COVID-19 model: 'covid_mathmodel_linearDD-Copy1.ipynb'
* Saturating density-dependent COVID-19 model for a super-high-density (like London): 'covid_mathmodel_saturating_London-Copy1.ipynb'
* Saturating density-dependent COVID-19 model for an average-England-density: 'covid_mathmodel_saturating_Average-Copy1.ipynb'
* Saturating density-dependent COVID-19 model for a low-density (like Cornwall): 'covid_mathmodel_saturating_Cornwall-Copy1.ipynb'
* England's COVID-19 associated death data up until August 1st 2020: 'cov_data_Aug.txt'
* R script needed to perform glm to show mortality is a function of population density: 'glm_deaths_density.R'

## Status
Project is: _finished_

## Contact
Laith Yakob (https://www.lshtm.ac.uk/aboutus/people/yakob.laith) - feel free to contact me!