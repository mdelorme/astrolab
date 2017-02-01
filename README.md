# astrolab

## Description

This repository is meant for the Introduction to Astronomy class at University of Surrey.
It mainly consists of two files setting up the necessary parameters to be able to use ds9, imexam and the fitting functions in the moons of Jupiter project and an additional Fortran file used for the final fitting part.

* `astro_lab.py` is a Python script setting up the environment to work on the Surrey students account. It also contains fitting routines and general "glue code" to link ds9 and imexam.
* `astro_lab.sh` is the main shell script that starts an ipython session while loading `astro_lab.py`.
* `fitting.f90` is the fitting subroutine presented in the project script. Use this only is you want to use Fortran to fit the final parameters in the project. You can also use Python to do the fit if you prefer so.

## Usage

For your first checkout, clone the rep in a subdirectory :

```shell
git clone https://github.com/mdelorme/astro_lab
cd astro_lab
```

If you have already cloned the repository then check before every run that all the data are up to date :

```shell
git pull
```

Then start the iptyhon session by running the main script :

```shell
./astro_lab.sh
```

You can now proceed with your project.
