##########################################
# generate equally spaced lambda windows
function gen_lambdas {
        cat <<EOF > gen_lambda.py
#!/usr/bin/env python

import numpy as np
import sys

nlam = sys.argv[1]; nlam = float (nlam)
a = np.linspace(0, 1, num = nlam, dtype = float)
for x in range(len(a)):
    print("{:.8f}".format(a[x])),

EOF

local nlambda=$1
local lams=$(python2.7 gen_lambda.py ${nlambda})
echo $lams
rm -rf gen_lambda.py
}

