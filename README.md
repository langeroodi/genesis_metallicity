# genesis-metallicity

non-parametric gas-phase metallicity and electron temperature estimation

Installation
-------
``genesis_metallicity`` is pip installable:

```bash
pip install genesis_metallicity
```

Examples
-------
First; an example of the "direct-method" metallicity estimation. Each emission line is provided as a list with two entries, where the first entry is the line flux and the second entry is its uncertainty. Note that a high-significance measurement of the [O III]4363 line is provided here, which enables the direct-method measurement. The code automatically opts to measure the metallicity through the direct-method when this line is available. Defining the "object" and providing the redshift are optional, and are only meant for bookkeeping purposes.

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
```

Now, let's test the "strong-line" metallicity estimator on the same object. The input dictionary remains identical, except for removing the [O III]4363 entry entry.

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
```
