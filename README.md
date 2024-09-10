# genesis_metallicity
 non-parametric gas-phase metallicity and electron temperature estimation

Here's an example of the direct-method metallicity estimator. Note that a signal-to-noise measurement of the [O III]4363 line is provided, and as a result, the code automatically choses the to measure the metallicity through the direct-method. Each line is inputed as a python list, where the first entry is the line flux and the second entry is its uncertainty. 

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
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

Now, let's test the strong-line metallicity estimator on the same object. For this purpose, it suffices to remove the [O III]4363 entry. 

```python
from genesis_metallicity.genesis_metallicity import genesis_metallicity

object = 'JADES_3675'

input_dict = {}
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
