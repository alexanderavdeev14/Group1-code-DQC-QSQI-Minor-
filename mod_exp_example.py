# Import stuff
import os
from getpass import getpass
from qiskit import BasicAer
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import execute
from quantuminspire.qiskit import QI
from quantuminspire.credentials import save_account
import numpy as np
import random as rand
from math import gcd
import math
import random

from DistributedFunctions import *
from AdderFunctionsDistr import *

# If you use token leave it as it is
QI_EMAIL = os.getenv('QI_EMAIL')
QI_PASSWORD = os.getenv('QI_PASSWORD')
QI_URL = os.getenv('API_URL', 'https://api.quantum-inspire.com/')

# Enter your token if needed and de-comment line below
# save_account('')


# Authenication procedure. Plug it in all your expermiments to establish connection between you and QI
authentication = get_authentication()
QI.set_authentication(authentication, QI_URL)
qi_backend = QI.get_backend('QX single-node simulator')


SetCounters()


# Specify the numbers you want to multiply
a_number = 1
x_number = 2

# Specify the mod value
N_number = 65535

# max n bits per register
n_max = 16

# Checks weather N is larger than a
if N_number < a_number or N_number < x_number:
    raise NameError("N must be larger than a")

# Converts N_number to n_string

N_string = str(bin(N_number))[2:]
n = len(N_string)
print(n)

# Converts the decimal representation to binary
a_string = binary_string(a_number, n)
x_string = binary_string(x_number, n)


# The input quantum bits
ones_state = list(range(n))
zero_b = list(range(n, 2*n))
a = list(range(2*n, 3*n))
carry_register = list(range(3*n, 4*n))
x_qu = list(range(4*n, 5*n))
N_reg = list(range(5*n, 6*n))
temp_reg = 6*n
ob = 6*n + 1



n_tot = 100
n_c = 49
c_name = "Mod"

qc, loc, clr = GenerateQCircuit(n_tot, n_c, c_name)

add_classical = ClassicalRegister(n_tot)
qc.add_register(add_classical)

# Creates the binary value of x_number in quantum circuit for qubits x
for x in range(n):
    if int(x_string[x]) == 1:
        qc.x(x_qu[n-x-1])

# Creates the binary value of N_number in quantum circuit for qubits N_reg
for x in range(n):
    if int(N_string[x]) == 1:
        qc.x(N_reg[n-x-1])

qc.x(ones_state[0])

# Controlled a*x mod N
loc = mod_exp(a, ones_state, zero_b, ob, x_qu, carry_register, a_number, N_reg, N_string, temp_reg, qc, clr, n_c, loc)


NormalCounter, ECounter, OneCounter = GetCounters()
print("number of 'normal' 2 - qubit gates")
print(NormalCounter)
print("number of 'normal' 1 - qubit gates")
print(OneCounter)
print("number of entanglement actions")
print(ECounter)

