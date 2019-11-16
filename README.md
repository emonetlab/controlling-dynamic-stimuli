# README

This repository contains code to go along with the paper

[Controlling and measuring dynamic odorant stimuli in the laboratory](https://jeb.biologists.org/content/early/2019/11/06/jeb.207787.abstract)

by Srinivas Gorur-Shandilya, Carlotta Martelli, Mahmut Demir and Thierry Emonet. 


# Contents


## Code to reproduce papers in figure

Code to reproduce the primary figures in the paper is contained in the `paper-figures/` folder. Once you have downloaded all the code and data (see below), you can run each script to make the figure as it appears in the paper. 

For example,

```
fig_calibration
```

will generate Figure 3 in the paper. 


## Code to tune PID parameters of your Alicat MFC

A toolbox to talk to your Alicat MFC and automatically tune PID values is included in this repository. Click [here](https://github.com/emonetlab/controlling-dynamic-stimuli/tree/master/alicat-mfc-tools) for a quick overview on this. 

## A simple model of odorant puffs, with sliders to manipulate parameters 

## Data

To make it possible to reproduce some of the figures in the paper, we have included data that went into these figures. You will have to manually download it.

[Download data here]()

# Installation 


## Prequisities 

0. MATLAB 
1. [mtools](https://github.com/sg-s/srinivas.gs_mtools/)
2. [data-manager](https://github.com/sg-s/data-manager/)

## Get the code 

If you're using git, you can 

```bash
git clone https://github.com/sg-s/srinivas.gs_mtools/
git clone https://github.com/sg-s/data-manager
git clone https://github.com/emonetlab/controlling-dynamic-stimuli
```

to get all code and dependencies. 



## Tell MATLAB where the data is 



```matlab
setpref('controlling_odor','data_loc','/path/to/data/')
```

# Usage

## Reproduce figures from the paper




# License 

GPL v3
