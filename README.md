# protnet
Protein network inference from RPPA data

# Scripts
1) run_methods.r 
- Run all 13 methods with the appropriate arguments

2) choose_best_parameter.r 
- Determine optimal parameter values and resulting edge predictions for methods that require parameter optimization
- The ground truth is the set of edges extracted from Pathway Commons using PERA.
- Computing AUPR is done by the function longPrecisionRecall. For sparse methods, this function will randomly assign orders for zero-weight edges.
