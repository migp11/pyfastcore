# pyfastcore

A python implementation of the context-specific model extraction method FastCore by Vlassis et al. (2014)

* Vlassis, N., Pacheco, M. P., & Sauter, T. (2014). PLoS Computational Biology, 10(1) [[Full Article]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424)

## INSTALATION

You can install pyfastcore using:
`python setup.py install`

of via pip:
`pip install pyfastcore`

## USAGE EXAMPLE
```
from cobra.test import create_test_model
from pyfastcore import Fastcore

# Loading a toy model of E. coli from cobra.test package
model = create_test_model('textbook')
# Define the list of core reactions
core_reactions = ['Biomass_Ecoli_core', 'ATPM']
# Setting the penalty of exchange fluxes to 0
penalties = {}
for r in model.exchanges:
    penalties[r.id] = 0

# Creating a fastcore solver instnace
fc_builder = Fastcore(model, core_reactions,
                      penalties=penalties,
                      default_penalty=10,
                      debug_mode=True)

# Rnunning fastcore
fc_builder.fast_core()

# checking the list of reaction in the consistent network found
consistent_subnetwork = fc_builder.consistent_subnetwork
print("Consistent subnetworksize set size", len(consistent_subnetwork))
print("Context specific core:")
print(consistent_subnetwork)

# creating a cobra model for the consistent network found
print(f"Building context-specific model for {model.id}")
cs_model = fc_builder.build_context_specific_model()

# Running and FBA using subnetwork model 
print("Running FBA on CS-model")
sol = cs_model.optimize()
print(cs_model.summary())

```


## Citation
If you use this package cite:
* Ponce-De-Leon, M. et al. (2015) Consistency Analysis of Genome-Scale Models of Bacterial Metabolism: A Metamodel Approach. PloS one, 10, e0143626.

