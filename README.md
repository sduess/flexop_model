Pazy
- Case 2
- Gust vane


For meeting:
- overnight: beta newmark factor in closed-loop simulations
- write paper
- gust vanes m=8

add [aero] ailerons controllable flag for closed-loop simulations
- add input for user to specfiy ailerons to be controllable, i.e. allowing deflections set by network


# FLEXOP Model for SHARPy

This repository contains the [SHARPy](http://github.com/imperialcollegelondon/sharpy) version of the equivalent FLEXOP model based on the European [FLEXOP project](https://flexop.eu/) with a detailed implementation description included in [1]. 

![alt text](https://github.com/sduess/flexop_model/doc/FLEXOP_white.eps)

## Installation

Clone the repository to your local computer. It is intended for use within SHARPy so you must ensure that SHARPy is 
properly installed in your system, including the conda environment provided.

With SHARPy and its environment installed, the only required step to use this model is to run from the terminal where
you are running your scripts

```bash
source <path-to-flexop-model>/bin/flexop_vars.sh
```
 
This will append the location of the `flexop-model` folder to the system's path, such that your code will be able to find
the modules.

## Using the model

In your SHARPy model generation script, import the flexop model:

```python
from flexop_model import aircraft
```

The first initialisation step takes the form of

```python
flexop_model = FLEXOP(case_name, case_route, output)
```
followed by for an e.g. aeroelastic simulation

```python
flexop_model.init_aeroelastic(flexop_settings)
```
where the `flexop_settings` refer to various parameter settings regarding discretisation, modelling assumptions, etc. A list of these
settings is provided in the source file.

For usual aeroelastic  simulation, the input files for SHARPy (including structure, aero, and SHARPy h5 files) are generated with

```python
flexop_model.generate()
flexop_model.create_settings(settings)
```

and SHARPy can finally be run with


```python
flexop_model.run()
```

## References

[1] Duessler, S., & Palacios. A Control Design Framework for Gust Load Alleviation in more Flexible Aircraft. IFASD, 2022.

## Contact

If you have any questions and want to get in touch, 
[contact us](https://www.imperial.ac.uk/aeroelastics/people/duessler/).

If you have any questions on how to use the model or find any bugs please file an issue. 