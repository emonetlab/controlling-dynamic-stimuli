# README

This repository contains code to go along with the paper

[Controlling and measuring dynamic odorant stimuli in the laboratory](https://jeb.biologists.org/content/early/2019/11/06/jeb.207787.abstract)

by Srinivas Gorur-Shandilya, Carlotta Martelli, Mahmut Demir and Thierry Emonet. 


# Contents


## Code to reproduce papers in figure


## Code to tune PID parameters of your Alicat MFC

## A simple model of odorant puffs, with sliders to manipulate parameters 

## Data

# Installation 


## Prequisities 

0. MATLAB 
1. [mtools](https://github.com/sg-s/srinivas.gs_mtools/)
2. [data-manager](https://github.com/sg-s/data-manager/)

## Get the code 

## Tell MATLAB where the data is 

1. Install [this toolbox]() and make sure it's on your path
2. Download and install 
2. Download necessary data from HERE and put it somewhere 
3. Tell MATLAB where the data is:

```matlab
setpref('controlling_odor','data_loc','/path/to/data/')
```

# Usage

## Reproduce figures from the paper




1. `/paper-figures/` code to generate figures for the methods paper on how to deliver dynamic odor stimuli. 
2. `natStimBuilder` a class to construct Naturalistic stimuli using MFCs 
3. 

You should be able to exactly reproduce every figure in the paper using this repository. Scripts in `paper-figures` should make figures for the paper exactly as you see them. 

# Installation 

The best way to install this repo is through my package manager: 

```
urlwrite('http://srinivas.gs/install.m','install.m'); 
install sg-s/srinivas.gs_mtools
install sg-s/data-manager   

% this url will work after the paper is accepted. 
% Before that, use git to clone this repo. 
install sg-s/how-to-deliver-odor-stimuli        
```

# Get the data and link it

Grab the data. Here is a link to download all the data you need: [INSERT LINK HERE]

Then, link the data using `data-manager`

```matlab
rehash(dataManager,'path/to/where/you/copied/the/data')
```

# make PDFs of figures

You can now navigate to `paper-figures/` and make PDFs of any figure using:

```matlab
makePDF('name-of-script')
```