#!/usr/bin/env python

import numpy as np
import sys

nlam = sys.argv[1]; nlam = float (nlam)
a = np.linspace(0, 1, num = nlam, dtype = float)
for x in range(len(a)):
    print("{:.8f}".format(a[x])),

