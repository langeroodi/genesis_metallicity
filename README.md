# genesis-metallicity

non-parametric gas-phase metallicity and electron temperature estimation

- __strong-line metallicity etstimation__ when temperature-sensitive emission lines are unavailable
- __direct-method metallicity measurement__ when the temperature-sensetive [O III]4363 line is detected (can also include [O II]7320,30 if available)
- __$${\rm O}^{+}$$ electron temperature estimation__ from that of the $${\rm O}^{++}$$ zone when direct measurements of the former are unavailable
- __dust reddening correction__ of the observed emission line fluxes when multiple Balmer lines are detected

Installation
-------
``genesis_metallicity`` is pip installable:

```bash
pip install genesis_metallicity
```

Examples
-------
### strong-line metallcitiy estimation

The following is an example of the "strong-line" metallicity estimation. The emission line measurements are imported in a python dictionary, where a python list with two items is entered for each emission line: the first item corresponds to the measured line flux and the second item corresponds to the flux uncertainty. Note that providing the object ID and redshift, as done below, are optional and only meant to assist with bookkeeping.

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
input_dict['redshift'] = 9.43
input_dict['OII']      = [7.269250230638606e-20, 5.7108025103430405e-21]
input_dict['Hdelta']   = [1.592676517913214e-19, 4.405126486743839e-21]
input_dict['Hgamma']   = [2.6671788939798604e-19, 5.274969458136735e-21]
input_dict['Hbeta']    = [6.447421960287729e-19, 6.48899753642406e-21]
input_dict['O4959']    = [1.0763985148795857e-18, 8.773228322203986e-21]
input_dict['O5007']    = [3.0628160287038502e-18, 1.0895841873808477e-20]
input_dict['Hbeta_EW'] = [158.728418552416, 13.991218097105634]

galaxy = genesis_metallicity(input_dict, object=object)
print(' -> strong-line metallicity:', galaxy.metallicity)
print(' -> strong-line metallicity (nominal value):', galaxy.metallicity.n)
print(' -> strong-line metallicity (standard deviation):', galaxy.metallicity.s)
```

### direct-method metallicity measurement

If the temperature-sensitive [O III]$$\lambda$$4363 emission line is detected, the gas-phase metallicity can be measured directly. This is demonstrated in the example below, where the input emission lines remain the same as in the previous example, but with the addition of the [O III]$$\lambda$$4363 flux.

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
input_dict['redshift'] = 9.43
input_dict['OII']      = [7.269250230638606e-20, 5.7108025103430405e-21]
input_dict['Hdelta']   = [1.592676517913214e-19, 4.405126486743839e-21]
input_dict['Hgamma']   = [2.6671788939798604e-19, 5.274969458136735e-21]
input_dict['O4363']    = [7.1092219385595e-20, 5.0986852540807764e-21]
input_dict['Hbeta']    = [6.447421960287729e-19, 6.48899753642406e-21]
input_dict['O4959']    = [1.0763985148795857e-18, 8.773228322203986e-21]
input_dict['O5007']    = [3.0628160287038502e-18, 1.0895841873808477e-20]
input_dict['Hbeta_EW'] = [158.728418552416, 13.991218097105634]

galaxy = genesis_metallicity(input_dict, object=object)
print(' -> direct-method metallicity:', galaxy.metallicity)
print(' -> direct-method metallicity (nominal value):', galaxy.metallicity.n)
print(' -> direct-method metallicity (standard deviation):', galaxy.metallicity.s)
```

### electron temperature estimation

It’s often desirable to estimate the $${\rm O}^{+}$$ electron temperature (t2) from the directly measured $${\rm O}^{++}$$ electron temperature (t3). This is particularly the case at high redshifts, where the t3 can be measured directly from the [O III]$$\lambda$$4363 line, while the t2 cannot be measured directly because the [O II]7320,30 doublet is often too faint or redshifted out of coverage. In such cases, the t2 estimations are carried out automatically by ```genesis_metallicity```, and the measured electron temperatures can be accessed as shown below.

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
input_dict['redshift'] = 9.43
input_dict['OII']      = [7.269250230638606e-20, 5.7108025103430405e-21]
input_dict['Hdelta']   = [1.592676517913214e-19, 4.405126486743839e-21]
input_dict['Hgamma']   = [2.6671788939798604e-19, 5.274969458136735e-21]
input_dict['O4363']    = [7.1092219385595e-20, 5.0986852540807764e-21]
input_dict['Hbeta']    = [6.447421960287729e-19, 6.48899753642406e-21]
input_dict['O4959']    = [1.0763985148795857e-18, 8.773228322203986e-21]
input_dict['O5007']    = [3.0628160287038502e-18, 1.0895841873808477e-20]
input_dict['Hbeta_EW'] = [158.728418552416, 13.991218097105634]

galaxy = genesis_metallicity(input_dict, object=object)
print('--------------------------------------------------------------------------')
print(' -> O++ electron temperature' )
print(' -> direct-method te(OIII) [K]:', galaxy.t3)
print(' -> direct-method te(OIII) [K] (nominal value):', galaxy.t3.n)
print(' -> direct-method te(OIII) [K] (standard deviation):', galaxy.t3.s)
print('--------------------------------------------------------------------------')
print(' -> O+ electron temperature' )
print(' -> direct-method te(OII) [K]:', galaxy.t2)
print(' -> direct-method te(OII) [K] (nominal value):', galaxy.t2.n)
print(' -> direct-method te(OII) [K] (standard deviation):', galaxy.t2.s)
print('--------------------------------------------------------------------------')
```

### dust reddening correction

By default, ```genesis_metallicity``` assumes that the input emission line fluxes are the observed values without any reddening correction (this behavior can be modified by the user; see the end of this paragraph). As such, the input line fluxes are automatically corrected for dust reddening. The meaured $$A_{\rm V}$$ as well as the reddening-corrected emission line fluxes can be accessed as shown below. Note that the ```galaxy.reddening_corrected_lines``` outpout is a dictionary, where for example the reddening corrected flux and flux uncertainty of $$H\beta$$ can be accessed as ```galaxy.reddening_corrected_lines['Hbeta'].n``` and ```galaxy.reddening_corrected_lines['Hbeta'].s```, respectively. If the input emission lines are already reddening-corrected, the reddening correction can be switched off by setting ```correct_extinction=False``` in the main ```genesis_metallicity``` function call; i.e., by calling the main routine as ```genesis_metallicity(input_dict, object=object, correct_extinction=False)```.

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
input_dict['redshift'] = 9.43
input_dict['OII']      = [7.269250230638606e-20, 5.7108025103430405e-21]
input_dict['Hdelta']   = [1.592676517913214e-19, 4.405126486743839e-21]
input_dict['Hgamma']   = [2.6671788939798604e-19, 5.274969458136735e-21]
input_dict['O4363']    = [7.1092219385595e-20, 5.0986852540807764e-21]
input_dict['Hbeta']    = [6.447421960287729e-19, 6.48899753642406e-21]
input_dict['O4959']    = [1.0763985148795857e-18, 8.773228322203986e-21]
input_dict['O5007']    = [3.0628160287038502e-18, 1.0895841873808477e-20]
input_dict['Hbeta_EW'] = [158.728418552416, 13.991218097105634]

galaxy = genesis_metallicity(input_dict, object=object)
print(' -> Av:', galaxy.Av)
print(' -> reddening-corrected line fluxes:')
print(galaxy.reddening_corrected_lines)
```
