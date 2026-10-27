"""Sorted celerite inputs with predictions restored to caller sample order."""
import numpy as np


def compute_sorted(gp, times, errors):
    order = np.argsort(times, kind='stable')
    gp.compute(np.asarray(times)[order], yerr=np.asarray(errors)[order])
    return order


def predict_original(gp, residuals, training_order, times, return_var=False):
    order = np.argsort(times, kind='stable')
    inverse = np.argsort(order)
    prediction = gp.predict(np.asarray(residuals)[training_order], np.asarray(times)[order],
                            return_cov=False, return_var=return_var)
    if return_var:
        return tuple(np.asarray(value)[inverse] for value in prediction)
    return np.asarray(prediction)[inverse]
