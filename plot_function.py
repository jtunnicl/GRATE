#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt


x, y = np.loadtxt(sys.argv[1], skiprows=0, delimiter=",", unpack=True)
plt.plot(x, y, linewidth=2)
plt.xlabel("XS.Width")
plt.ylabel("Qb_cap")
plt.grid(True)
if len(sys.argv) > 2:
    xlim = int(sys.argv[2])
    plt.xlim(0, xlim)
plt.tight_layout()
plt.show()
