## Abstract
This repository includes all of the code which was developed during my master's thesis.

For the thesis itself and its corresponding source code visit the repository [jpmvferreira/msc-thesis](https://github.com/jpmvferreira/msc-thesis).

## Repository Structure
```
.
├── analyzed -> includes results which are the output of processing data files, mostly made up of plots
├── aux -> includes auxiliary scripts to process data
├── config -> includes configuration files to run simplifiedmc
├── cosmology -> includes .py files that define a given cosmology to be used with gwcatalog
├── data -> includes all data files in a .csv format
├── jobs -> auxiliary script to run batch jobs
├── LICENSE.md
├── misc -> miscellaneous scripts which were used either as auxiliary programs, sketching ideas or to explore new things. generally speaking they are not commented and might be very hacky.
├── model -> includes Stan model files to run with simplifiedmc
├── output -> includes output of the MCMC's provided by simplifiedmc
└── venv -> virtual environment used throughout the work
```

## External References:
- The Stan programming language: [mc-stan.org](https://mc-stan.org/)
- The SimplifiedMC Python package: [jpmvferreira/simplifiedmc](https://github.com/jpmvferreira/simplifiedmc)
- The GWCatalog Python package: [jpmvferreira/gwcatalog](https://github.com/jpmvferreira/gwcatalog)
