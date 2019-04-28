#!/usr/bin/env python

import os
import sys

import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    print("Usage: plot_function.py FILENAME")
    sys.exit(1)

if not os.path.exists(sys.argv[1]):
    print("Error: input file does not exist: %s" % sys.argv[1])
    sys.exit(1)

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
