#!/usr/bin/env python
# coding: utf-8

# In[17]:


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

# Prepare variables
N = 127
n = 7

# Prepare register qubits numbers
input_1 = list(range(7))
input_2 = list(range(7,14))

a_reg = list(range(14,21))
x_reg = list(range(21,28))
n_reg = list(range(28,35))
c_reg = list(range(35,42))
ob = 43
temp = 44

while True:

    a_number = N-1
    print("Random guess:", a_number)
    g = gcd(a_number, N)
    print("GCD:", g)
    # When we got our guess right first time
    if g != 1 or N == 1:
        factor1 = g
        factor2 = int(N / g)
        print('N is factorized by ', factor1, ' and ', factor2)
        break
    else:
        print("iam qunatum")
        # generate circuit
        n_c = 4
        n_tot = 45
        qc_name = "Distributed Shor"
        qc, loc, clr = GenerateQCircuit(n_tot, n_c, qc_name)

        add_classical = ClassicalRegister(n_tot)
        qc.add_register(add_classical)

        # Prep x register
        for i in range(n):
            loc = ApplyGate(x_reg[i], x_reg[i], qc, clr, n_c, loc, "h")

        # Prep N_register
        n_string = binary_string(N, n)
        loc = prep_register(n_string, n_reg, qc, clr, n_c, loc)

        # Prep input 1 to value 1
        loc = ApplyGate(input_1[0], input_1[0], qc, clr, n_c, loc, "x")


        # Modular exponentiation
        mod_exp(a_reg, input_1, input_2, ob, x_reg, c_reg, a_number, n_reg, n_string, temp, qc, clr, n_c, loc)

        # Inverse QFT
        InverseQFT(x_reg, qc, clr, n_c, loc)

        # Measure all
        measure_qubits(qc, n_c, loc, input_2, x_reg)
        
        break
        # QI backend non local
        # qi_job = execute(qc, backend=qi_backend, shots=1)
        #qi_result = qi_job.result()
        #histogram = qi_result.get_counts(qc)
        #print("Results from the non-local QI backend\n")
        #print('State\tCounts')
        #[print('{0}\t{1}'.format(state, counts)) for state, counts in histogram.items()]

        #for state, counts in histogram.items():
            #start_bit = len(state)
            # result_out = state[start_bit-n-2:start_bit-2]

        #print(state)
        #state_adjusted = ReturnState(loc, state, n_c)
        #print(state_adjusted)
        #r_tmp = state_adjusted[6:9]
        #r_tmp = r_tmp[::-1]
        #r = int(r_tmp, 2)



        # Reset the circuit

        # we want to find a r that a**r/2 = 1 Mod N
        #if r % 2 != 0:
            #print(" r mod 2 is not 0")
            # this continue sets it back to the while loop
            #continue
        #elif a_number ** (int(r / 2)) % N == -1 % N:
            #print("a**r/2 mod N =-1", a_number ** (int(r / 2)) % N)
            #continue
        #else:
            # here we found our correc  t r
            #factor1 = gcd(a_number ** int(r / 2) + 1, N)
            #factor2 = gcd(a_number ** int(r / 2) - 1, N)
            #if factor1 == N or factor2 == N:
                #continue
            #else:
                #print('N is factorized by', factor1, 'and', factor2)
                #break

            
NormalCounter, ECounter, OneCounter = GetCounters()
print("number of 'normal' 2-qubit gates")
print(NormalCounter)
print("number of 'normal' 1-qubit gates")
print(OneCounter)
print("number of entanglement actions")
print(ECounter)


# In[ ]:




