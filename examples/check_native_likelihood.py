"""Embedded Jupyter: %run -i /path/to/exostriker/examples/check_native_likelihood.py"""
import math
import numpy as np

result=fit.evaluate_native_likelihood()
reference=math.fsum(-.5*((result.observed-result.data_model)**2/result.variance
                        +np.log(result.variance)+np.log(2*np.pi)))
print(f'Native RV lnL:       {result.rv_loglik:.12f}')
print(f'Native transit lnL:  {result.transit_loglik:.12f}')
print(f'Native total lnL:    {result.loglik:.12f}')
print(f'NumPy check:         {reference:.12f}')
print(f'Native − NumPy:      {result.loglik-reference:.6g}')
for kind,name in ((0,'RV'),(1,'Transit')):
    for slot in np.unique(result.dataset[result.kind==kind]):
        count=np.count_nonzero((result.kind==kind)&(result.dataset==slot))
        print(f'{name} slot {slot}: N={count}, lnL={result.by_dataset[slot,kind]:.12f}')
