from subprocess import call
from time import time

time_id = str(int(time()) - 1479729127)

call_list = []

call_list.extend(["./Exe_GeneralThreeRegion"])
call_list.extend(["--ID", time_id])  # sim id
call_list.extend(["--CRL", "0.25"])  # cortical rest length
call_list.extend(["--CSC", "1e6"])   # cortical spring constant
call_list.extend(["--TRL", "0.01"])  # Transmembrane rest length
call_list.extend(["--TSC", "1e7"])   # Transmembrane spring constant
call_list.extend(["--AD", "2.0"])   # Adhesion modifier
call_list.extend(["--DI", "0.02"])   # Interaction dist for cell-cell forces
call_list.extend(["--RM", "10"])     # ReMesh frequency
call_list.extend(["--TS", "5000"])   # Numer of time steps

print call_list
call(call_list)
