import numpy as np
import pyqtgraph as pg
import sys 

class CustomAxisItem(pg.AxisItem):

    def tickStrings(self, values, scale, spacing):
        if self.logMode:
            return self.logTickStrings(values, scale, spacing)

        places = max(0, np.ceil(-np.log10(spacing*scale)))
        strings = []
        for v in values:
            vs = v * scale
            vstr = ("%%0.%df" % places) % vs
            strings.append(vstr)
        return strings

    def logTickStrings(self, values, scale, spacing):
       # estrings = ["%0.1g"%x for x in 10 ** np.array(values).astype(float) * np.array(scale)]
        estrings = ["%1g"%x for x in 10 ** np.array(values).astype(float)]
        if sys.version_info < (3, 0):
            # python 2 does not support unicode strings like that
            return estrings
        else:  # python 3+
            convdict = {"0": "⁰",
                        "1": "¹",
                        "2": "²",
                        "3": "³",
                        "4": "⁴",
                        "5": "⁵",
                        "6": "⁶",
                        "7": "⁷",
                        "8": "⁸",
                        "9": "⁹",
                        }
            dstrings = []
            for e in estrings:
                if e.count("e"):
                    v, p = e.split("e")
                    sign = "⁻" if p[0] == "-" else ""
                    pot = "".join([convdict[pp] for pp in p[1:].lstrip("0")])
                    if v == "1":
                        v = ""
                    else:
                        v = v + "·"
                    dstrings.append(v + "10" + sign + pot)
                else:
                    dstrings.append(e)
            return dstrings



