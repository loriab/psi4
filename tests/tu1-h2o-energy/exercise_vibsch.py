import qcelemental as qcel
from qcelemental.models.v2 import Vibrations
from pathlib import Path
import json


vibrec = json.loads(Path("input.h2o.54437.vibrec").read_text())
# The mu array ordering can be non-deterministic with multi-threaded LAPACK, so check max value
omega = json.loads(vibrec["omega"])["data"]
#compare_values(1000., max(mu_data), "weird H mass from file", atol=1.e-6)

print(type(omega[0]))
print(omega)

Vibrations(omega=omega)

