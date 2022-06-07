## Quantum conductivity fit
Quantum conductivity fit procedure used and described in Pickem et al.'s 
'Resistivity saturation in Kondo insulators',
[Commun. Phys. 4, 226 (2021)](https://www.nature.com/articles/s42005-021-00723-z).
The fit allows extracting the size of the gap ![equation](https://latex.codecogs.com/gif.latex?\Delta)
and the scattering rate ![equation](https://latex.codecogs.com/gif.latex?\Gamma) from the low-temperature behavior of an experimental conductivity of a semiconductor or insulator.

## Script
`quantum_conductivity_fit.py` is a user script to fit experimental data with the two-level (2L)
conductivity equation

![equation](https://latex.codecogs.com/gif.latex?\sigma_{\hbox{\tiny&space;2L}}(T)\propto\frac{\beta}{\Gamma}\left[&space;\Re\Psi^\prime(z)-&space;\frac{\beta\Gamma}{2&space;\pi}\Re\Psi^{\prime\prime}(z)&space;\right].)

with

![equation](https://latex.codecogs.com/gif.latex?z=\frac{1}{2}&plus;\frac{\beta}{2\pi}(\Gamma&plus;i\frac{\Delta}{2}))

on a given dataset (here: `resistivity.dat`)

## Fit

Executing this script on the provided data results in the following graphical output
![](https://github.com/LinReTraCe/quantum_conductivity_fit/blob/main/quantum_fit_screenshot.png?raw=true)

### Acknowledgments
This project has been supported by the Austrian Science Fund (FWF) through project LinReTraCe P 30213-N36.
