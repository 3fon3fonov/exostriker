"""Exo-Striker axis labels: full numbers below ten integer digits."""
import math
from decimal import Decimal
import numpy as np


def tick_numbers(values, scale, spacing, logarithmic=False):
    values=np.asarray(values,dtype=float)
    if logarithmic:
        with np.errstate(over='ignore',invalid='ignore'):
            values=np.power(10.,values)
    values=values*scale
    # Log minor ticks have spacing=None; their labels use significant digits,
    # so the linear-axis decimal-place calculation is unnecessary.
    step=0. if logarithmic else abs(float(spacing)*float(scale))
    places=max(0,int(math.ceil(-math.log10(step)))) if step>0 and math.isfinite(step) else 0
    # Fractional spacings such as 0.25 need another decimal place.
    if step>0 and math.isfinite(step):
        for _ in range(12):
            unit=step*10.**min(places,308)
            if abs(unit-round(unit))<=1e-8*max(1.,abs(unit)):break
            places+=1
    strings=[]
    for value in values:
        if not math.isfinite(value):
            text=str(value)
        elif abs(value)>=1e9:
            text=np.format_float_scientific(value,unique=True,trim='-')
        elif logarithmic:
            text=format(Decimal(format(value,'.15g')),'f')
        else:
            text=format(value,f'.{places}f')
        if '.' in text and 'e' not in text:
            text=text.rstrip('0').rstrip('.')
        if text=='-0':text='0'
        strings.append(text)
    return strings
