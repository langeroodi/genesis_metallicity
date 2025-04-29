# genesis-metallicity

[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/langeroodi/genesis_metallicity/blob/main/LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-2409.07455-b31b1b.svg)](https://arxiv.org/abs/2409.07455)

non-parametric gas-phase metallicity and electron temperature estimation

- __strong-line metallicity etstimation__ when temperature-sensitive emission lines are unavailable
- __direct-method metallicity measurement__ when the temperature-sensetive [O III]4363 line is detected (can also include [O II]7320,30 if available)
- __O⁺ electron temperature estimation__ from that of the O⁺⁺ zone when direct measurements of the former are unavailable
- __dust reddening correction__ of the observed emission line fluxes when multiple Balmer lines are detected

*Calibration data: The calibration data will be publicly released upon the journal publication of the associated paper. In the meantime, do not hesitate to get in touch if you are interested in using this data in your work: danial.langeroodi@nbi.ku.dk*

*Real-valued 5D and 4D images of our electron temperature and gas-phase metallicity KDEs are available in the [images/](images/) directory.*

Installation
-------
``genesis_metallicity`` is pip installable:

```bash
pip install genesis_metallicity
```

It is recommended to keep the package up to date to use the largest available calibration sample. If you have an older version installed, you can install the most recent version by:

```bash
pip uninstall -y genesis_metallicity
pip install --upgrade genesis_metallicity
```

You can check the installed version by

```python
import genesis_metallicity
print(genesis_metallicity.__version__)
```

Examples
-------
### strong-line metallcitiy estimation

The following is an example of "strong-line" metallicity estimation. The emission line measurements are imported into a Python dictionary, where a Python list with two items is entered for each emission line: the first item corresponds to the measured line flux, and the second item corresponds to the flux uncertainty. Note that providing the object ID and redshift, as done below, is optional and intended only to assist with bookkeeping.

Providing EW(Hβ) measurements is recommended but not strictly required. If EW(Hβ) measurements are not available, simply do not provide them in the dictionary below; the metallicity estimator will instead adopt an EW(Hβ)-independent calibration. This applies to other modules as well.

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

If the temperature-sensitive [O III]4363 emission line is detected, the gas-phase metallicity can be measured directly. This is demonstrated in the example below, where the input emission lines remain the same as in the previous example, but with the addition of the [O III]4363 flux.

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

It’s often desirable to estimate the O⁺ electron temperature (t2) from the directly measured O⁺⁺ electron temperature (t3). This is particularly the case at high redshifts, where the t3 can be measured directly from the [O III]4363 line, while the t2 cannot be measured directly because the [O II]7320,30 doublet is often too faint or redshifted out of coverage. In such cases, the t2 estimations are carried out automatically by ```genesis_metallicity```; the measured electron temperatures can be accessed as shown below.

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

By default, ```genesis_metallicity``` assumes that the input emission line fluxes are the observed values without any reddening correction (this behavior can be modified by the user; see the end of the next paragraph). As such, the input line fluxes are automatically corrected for dust reddening. The meaured V-band attenuation as well as the reddening-corrected emission line fluxes can be accessed as shown below. If the line of interest is included in the ```data/lines.py``` script, its flux can be included in the input dictionary for reddening correction. Otherwise, the line ID and rest-wavelength have to be added to ```data/lines.py``` by the user.

Note that the ```galaxy.reddening_corrected_lines``` outpout is a dictionary, where the reddening-corrected line fluxes are stored. For instance, the reddening-corrected flux and flux uncertainty of Hbeta can be accessed as ```galaxy.reddening_corrected_lines['Hbeta'].n``` and ```galaxy.reddening_corrected_lines['Hbeta'].s```, respectively. If the input emission lines are already reddening-corrected, the reddening correction can be switched off by setting ```correct_extinction=False``` in the main ```genesis_metallicity``` function call; i.e., by calling the main routine as ```genesis_metallicity(input_dict, object=object, correct_extinction=False)```.

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

Citation
-------

If you use ```genesis-metallicity``` in your research, please reference the [associated paper](https://ui.adsabs.harvard.edu/abs/2024arXiv240907455L/abstract):

```bibtext
@ARTICLE{2024arXiv240907455L,
       author = {{Langeroodi}, Danial and {Hjorth}, Jens},
        title = "{Genesis-Metallicity: Universal Non-Parametric Gas-Phase Metallicity Estimation}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = sep,
          eid = {arXiv:2409.07455},
        pages = {arXiv:2409.07455},
archivePrefix = {arXiv},
       eprint = {2409.07455},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv240907455L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

If you use the direct-method module of ```genesis-metallicity```, in addition to the above please also reference ```PyNeb```:

```bitext
@ARTICLE{2015A&A...573A..42L,
       author = {{Luridiana}, V. and {Morisset}, C. and {Shaw}, R.~A.},
        title = "{PyNeb: a new tool for analyzing emission lines. I. Code description and validation of results}",
      journal = {\aap},
     keywords = {methods: numerical, atomic data, Hii regions, planetary nebulae: general, ISM: abundances, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2015,
        month = jan,
       volume = {573},
          eid = {A42},
        pages = {A42},
          doi = {10.1051/0004-6361/201323152},
archivePrefix = {arXiv},
       eprint = {1410.6662},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015A&A...573A..42L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
