import re
import numpy as np


def rd(p):
    L = open(p).read().split("\n")
    hdr = [l.split() for l in L if l.strip() and not l.startswith("#") and len(l.split()) == 3
           and re.fullmatch(r"[-\d.eE+]+", l.split()[0]) and l.split()[2].isdigit()]
    nt = int(hdr[0][2]); ph0 = float(hdr[1][0]); ph1 = float(hdr[1][1]); npf = int(hdr[1][2])
    rows = [l.split() for l in L if len(l.split()) == 4]
    d = np.array([[float(x) for x in r] for r in rows[-nt * npf:]]).reshape(nt, npf, 4)
    et = d[..., 0] + 1j * d[..., 1]; ep = d[..., 2] + 1j * d[..., 3]
    return np.linspace(0, 180, nt), np.linspace(ph0, ph1, npf), et, ep
