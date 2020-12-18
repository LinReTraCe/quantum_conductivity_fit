## Quantum conductivity fit
Quantum conductivity fit procedure found in 'Resistivity saturation in Kondo insulators':
[arXiv:2008.05846](https://arxiv.org/pdf/2008.05846.pdf)

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
