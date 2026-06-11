import numpy as np

def periodogram(t,x,Pmin,Pmax,n=1000):

    P=np.linspace(Pmin,Pmax,n)
    freqs=1/P
    power = []

    x=x-np.mean(x)
    t=t-min(t)

    for f in freqs:
        s = np.sum(x * np.exp(-2j*np.pi*f*t))
        power.append(np.abs(s)**2)

    power=np.array(power)/max(power)
    periods = 1/freqs
    return periods,power