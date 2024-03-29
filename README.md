# Controlling and measuring dynamic odorant stimuli

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

In our paper, we described a simple model that can reproduce the rich kinetics observed in identically delivered odorant stimuli. To get a feel for the model, you can interact with the model and manipulate sliders that are linked to every parameter in the model. The interactive demo looks something like this:

![](https://user-images.githubusercontent.com/6005346/68999059-e643d100-0888-11ea-8e06-f6d97d9cca34.png)

To start this, run the `puff_explorer.m` script in the `puff-model` folder. 


## natStimBuilder

natStimBuilder is a class that will help you iteratively tune control signals to match some desired stimulus statistics, as in Figure 3. 

It's in the `src` folder. 


## Data

To make it possible to reproduce some of the figures in the paper, we have included data that went into these figures. You will have to manually download it.

[Download data here](https://github.com/emonetlab/controlling-dynamic-stimuli/releases/download/v19.11.16/data.zip)

# Installation 


## Prequisities 

0. MATLAB 
1. [mtools](https://github.com/sg-s/srinivas.gs_mtools/) (Helper functions)
2. [data-manager](https://github.com/sg-s/data-manager/) (To link data)
3. [puppeteer](https://github.com/sg-s/puppeteer/) (For interactive manipulation of models)

## Get the code 

If you're using git, you can 

```bash
git clone https://github.com/sg-s/srinivas.gs_mtools/
git clone https://github.com/sg-s/data-manager
git clone https://github.com/sg-s/puppeteer
git clone https://github.com/emonetlab/controlling-dynamic-stimuli
```

to get all code and dependencies. 



## Tell MATLAB where the data is 


Navigate to the folder where you data is, and use `dataManager` to link the data. What `dataManager` does is search these files, hashes them, and builds a hash table that maps hashes to paths. That way, the scripts you want to run can request a certain piece of data, and `dataManager` can serve that up, and the primary script doesn't have to know where the data is. 

```
cd /path/to/where/you/downloaded/data
dm = dataManager;
dm.rehash
```

# Usage

## Reproduce figures from the paper


# Bugs? Issues?

Open an issue [here](https://github.com/emonetlab/controlling-dynamic-stimuli/issues)

# License 

GPL v3
