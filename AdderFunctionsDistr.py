from quantuminspire.credentials import load_account, get_token_authentication, get_basic_authentication
from getpass import getpass
from DistributedFunctions import ApplyGate, FindLoc
import math
import numpy as np

def get_authentication():
    """ Gets the authentication for connecting to the Quantum Inspire API."""
    token = load_account()
    if token is not None:
        return get_token_authentication(token)
    else:
        if QI_EMAIL is None or QI_PASSWORD is None:
            print('Enter email:')
            email = input()
            print('Enter password')
            password = getpass()
        else:
            email, password = QI_EMAIL, QI_PASSWORD
        return get_basic_authentication(email, password)


def binary_string(number, n):

    strng = str(bin(number))[2:]

    diff = n - len(strng)

    for _ in range(diff):
        strng = "0" + strng

    return strng


def str_same_length(*args):
    strings = []

    # Loads strings
    for strng in args:
        strings.append(strng)

    # Determines the length
    n = len(strings)

    largest = 0
    index = 0

    for i in range(n):
        curr_str_length = len(strings[i])
        if curr_str_length > largest:
            largest = curr_str_length
            index = i

    for i in range(n):
        if index != i:
            curr_str_length = len(strings[i])
            for j in range(largest - curr_str_length):
                strings[i] = "0" + strings[i]

    return strings


def prep_register(string, register, qc, clr, n_c, loc):

    n = len(string)

    for i in range(n):
        if string[i] == "1":
            loc = ApplyGate(register[n - 1 - i], register[n - 1 - i],  qc, clr, n_c, loc, "x")

    return loc


def modInverse(a, m):
    m0 = m
    y = 0
    x = 1

    if (m == 1):
        return 0

    while (a > 1):
        # q is quotient
        q = a // m

        t = m

        # m is remainder now, process
        # same as Euclid's algo
        m = a % m
        a = t
        t = y

        # Update x and y
        y = x - q * y
        x = t

        # Make x positive
    if (x < 0):
        x = x + m0

    return x


def q_sum(a, b, c, qc, clr, n_c, loc):
    """ Sums the a, b and c and stores the result in c
        The following truth table:
            INPUT   OUTPUT
            A B C   A B C
            0 0 0   0 0 0
            0 0 1   0 0 1
            0 1 0   0 1 1
            0 1 1   0 1 0
            1 0 0   1 0 1
            1 0 1   1 0 0
            1 1 0   1 1 0
            1 1 1   1 1 1

        Arguments:
            a - first input qubit
            b - second input qubit
            c - third input line and out put line for the result
            qc - quantum circuit in question
    """
# CNOT of first line and third line
    loc = ApplyGate(a, c, qc, clr, n_c, loc, 'cnot')

# CNOT of the second line and the third line
    loc = ApplyGate(b, c, qc, clr, n_c, loc, 'cnot')
    
    return loc


# In[ ]:


def toffoli(a, b, c, qc, clr, n_c, loc):
    """ Toffoli gate made only out of 2 qubit gates

        Arguments:
            a - first control line qubit
            b - second control line qubit
            c - target qubit
            qc - quantum circuit in question

    """

    loc = ApplyGate(c, c, qc, clr, n_c, loc, 'h')
    loc = ApplyGate(b, c, qc, clr, n_c, loc, 'cnot')
    loc = ApplyGate(c, c, qc, clr, n_c, loc, 'tdg')

    loc = ApplyGate(a, c, qc, clr, n_c, loc, 'cnot')
    loc = ApplyGate(c, c, qc, clr, n_c, loc, 't')
    loc = ApplyGate(b, c, qc, clr, n_c, loc, 'cnot')

    loc = ApplyGate(c, c, qc, clr, n_c, loc, 'tdg')
    loc = ApplyGate(a, c, qc, clr, n_c, loc, 'cnot')
    loc = ApplyGate(b, b, qc, clr, n_c, loc, 't')

    loc = ApplyGate(c, c, qc, clr, n_c, loc, 't')
    loc = ApplyGate(a, b, qc, clr, n_c, loc, 'cnot')
    loc = ApplyGate(c, c, qc, clr, n_c, loc, 'h')

    loc = ApplyGate(a, a, qc, clr, n_c, loc, 't')
    loc = ApplyGate(b, b, qc, clr, n_c, loc, 'tdg')
    loc = ApplyGate(a, b, qc, clr, n_c, loc, 'cnot')
    
    return loc


# In[ ]:


def carry(q0, q1, q2, q3, qc, clr, n_c, loc):
    """ Calculates the resulting carry of the addition q0 + q1 + q2 and stores it in q3
        Truth table:
        Q0	Q1	Q2 	Q3		Q0	Q1 	Q2	Q3
        0	0	0	0		0	0	0	0
        0	0	0	1		0	0	0	1
        0	0	1	0		0	0	1	0
        0	0	1	1		0	0	1	1
        0	1	0	0		0	1	1	0
        0	1	0	1		0	1	1	1
        0	1	1	0		0	1	0	1
        0	1	1	1		0	1	0	0
        1	0	0	0		1	0	0	0
        1	0	0	1		1	0	0	1
        1	0	1	0		1	0	1	1
        1	0	1	1		1	0	1	0
        1	1	0	0		1	1	1	1
        1	1	0	1		1	1	1	0
        1	1	1	0		1	1	0	1
        1	1	1	1		1	1	0	0

        Arguments:
            q0 - the first inbut qubit
            q1 - the second input qubit
            q2 - the third input qubit
            q3 - the output qubit
            qc - the quantum circuit in question
    """
    # Tofolli of the second , third and fourth line
    loc = toffoli(q1, q2, q3, qc, clr, n_c, loc)

    # CNOT of the second and third line
    loc = ApplyGate(q1, q2, qc, clr, n_c, loc, 'cnot')

    # Tofolli of first, third and fourth line
    loc = toffoli(q0, q2, q3, qc, clr, n_c, loc)
    
    return loc


# In[ ]:


def carry_inverse(q0, q1, q2, q3, qc, clr, n_c, loc):
    """ Calculates the reverse operation of carry
        Truth table:
        Q0	Q1	Q2	Q3		Q0	Q1	Q2	Q3
        0	0	0	0		0	0	0	0
        0	0	0	1		0	0	0	1
        0	0	1	0		0	0	1	0
        0	0	1	1		0	0	1	1
        0	1	0	0		0	1	1	1
        0	1	0	1		0	1	1	0
        0	1	1	0		0	1	0	1
        0	1	1	1		0	1	0	1
        1	0	0	0		1	0	0	0
        1	0	0	1		1	0	0	1
        1	0	1	0		1	0	1	1
        1	0	1	1		1	0	1	0
        1	1	0	0		1	1	1	1
        1	1	0	1		1	1	1	0
        1	1	1	0		1	1	0	1
        1	1	1	1		1	1	0	0

        Arguments:
            q0 - the first inbut qubit
            q1 - the second input qubit
            q2 - the third input qubit
            q3 - the output qubit
            qc - the quantum circuit in question
    """
    # Tofolli of first, third and fourth line
    loc = toffoli(q0, q2, q3, qc, clr, n_c, loc)

    # CNOT of the second and third line
    loc = ApplyGate(q1, q2, qc, clr, n_c, loc, 'cnot')

    # Tofolli of the second , third and fourth line
    loc = toffoli(q1, q2, q3, qc, clr, n_c, loc)
    
    return loc


def measure_qubits(qc, n_c, loc, *args):
    n_row = math.sqrt(n_c)

    for arg in args:
        n = len(arg)
        for j in range(n):
            measured = False
            mod_c = arg[j] % n_c
            comp_no_c = int((arg[j] - mod_c) / n_c)
            cloc_i, cloc_j = FindLoc(arg[j], loc[comp_no_c])
            q_measure = int((n_c + 1) * comp_no_c + n_row * cloc_i + cloc_j)

            qc.measure(q_measure, q_measure + len(loc) - comp_no_c)


def adder(a, b, ob, c, qc, clr, n_c, loc):

    # The length of register a
    n = len(a)

    # Raises error if the sizes of registers a and c are not of the same length
    if n != len(c):
        raise NameError("The register a and c must be of equal size")

    # Calculate all carries until c_n
    for x in range(n-1):
        loc = carry(c[x], a[x], b[x], c[x+1], qc, clr, n_c, loc)

    # Calculates the overflow carry and stores it in b_n+1 bit
    loc = carry(c[n-1], a[n-1], b[n-1], ob, qc, clr, n_c, loc)

    # Restores the b_n after carry operation
    loc = ApplyGate(a[n-1], b[n-1], qc, clr, n_c, loc, 'cnot')

    # Sums the result of a_n + b_n + c_n and stores it in b_n
    loc = q_sum(c[n-1], a[n-1], b[n-1], qc, clr, n_c, loc)

    # Here the c_n...c_0 registers are restored i.e. set to all zeros and restores the b_n-1 register after it was per-
    # sturbed by the carry operation, and the sum of the c_n-1...c_0 + a_n-1...a_0 + b_n-1...b0 is calculated and stored
    # in b_n-1.
    for x in range(n-1):
        loc = carry_inverse(c[n-2-x], a[n-2-x], b[n-2-x], c[n-1-x], qc, clr, n_c, loc)

        loc = q_sum(c[n-2-x], a[n-2-x], b[n-2-x], qc, clr, n_c, loc)
        
    return loc


def adder_inverse(a, b, ob, c, qc, clr, n_c, loc):
    """ Subtract a from b and stores it in b
        if a >= b then b = a-b
        if a < b then b = 2^n+1 - (b-a) where n is the size of register a. In this case the b_n+1 will be always 1.
    """

    # The length of register a
    n = len(a)

    # Raises error if the sizes of registers a and c are not of the same length
    if n != len(c):
        raise NameError("The register a and c must be of equal size")

    for x in range(n-1):
        loc = q_sum(c[x], a[x], b[x], qc, clr, n_c, loc)
        loc = carry(c[x], a[x], b[x], c[x+1], qc, clr, n_c, loc)

    loc = q_sum(c[n - 1], a[n - 1], b[n - 1], qc, clr, n_c, loc)

    loc = ApplyGate(a[n-1], b[n-1], qc, clr, n_c, loc, 'cnot')

    loc = carry_inverse(c[n - 1], a[n - 1], b[n - 1], ob, qc, clr, n_c, loc)

    for x in range(n-1):
        loc = carry_inverse(c[n-2-x], a[n-2-x], b[n-2-x], c[n-1-x], qc, clr, n_c, loc)
    
    return loc


def adder_mod(a, b, ob, c, n_reg, n_string, t, qc, clr, n_c, loc):
    """ Computes the a + b mod N and stores it in register b

        Arguments:
            a - the fist value from the summation (n qubits)
            b - the second value from the summation (n + 1 qubits, b_n+1 is the overflow bit)
            c - the carry register (n qubits)
            n_reg - the qubit representation of N
            n_string - string that represents the binary value of N. i.e. if N = 3 then n_string = "011"
            t - the temporary qubit where the inverse of the overflow bit b_n+1 will be stored
            qc - quantum circuit in question
    """

    # Computes the size of register a
    n = len(a)

    # Computes the a + b and stores it in b
    loc = adder(a, b, ob, c, qc, clr, n_c, loc)

    # Computes N - b if N is smaller than b the overflow bit in register b (b_(n+1)) will be set
    loc = adder_inverse(n_reg, b, ob, c, qc, clr, n_c, loc)

    # If the overflow bit is set the N register will be substracted from the result
    loc = ApplyGate(ob, ob, qc, clr, n_c, loc, 'x')
    loc = ApplyGate(ob, t, qc, clr, n_c, loc, 'cnot')
    loc = ApplyGate(ob, ob, qc, clr, n_c, loc, 'x')

    # Changes N register to all zero register if t is set
    for x in range(n):
        if int(n_string[x]) == 1:
            loc = ApplyGate(t, n_reg[n - x - 1], qc, clr, n_c, loc, 'cnot')

    # Computes a+b mod N and stores it in b
    loc = adder(n_reg, b, ob, c, qc, clr, n_c, loc)

    # Restores the N register back to it's original state
    for x in range(n):
        if int(n_string[x]) == 1:
            loc = ApplyGate(t, n_reg[n - x - 1], qc, clr, n_c, loc, 'cnot')

    # Computes (a+b) mod N - a and stores it in b
    loc = adder_inverse(a, b, ob, c, qc, clr, n_c, loc)

    # Rests t qubit to is's original state
    loc = ApplyGate(ob, t, qc, clr, n_c, loc, 'cnot')

    # Computes the result (a+b) mod N
    loc = adder(a, b, ob, c, qc, clr, n_c, loc)
    
    return loc


def adder_mod_inverse(a, b, ob, c, n_reg, n_string, t, qc, clr, n_c, loc):
    """ This code is the reverse of the adder_mod: input - (a + b) mod N, output - b

        Arguments:
            a - the fist value from the summation (n qubits)
            b - the second value from the summation (n + 1 qubits, b_n+1 is the overflow bit)
            c - the carry register (n qubits)
            n_reg - the qubit representation of N
            n_string - string that represents the binary value of N. i.e. if N = 3 then n_string = "011"
            t - the temporary qubit where the inverse of the overflow bit b_n+1 will be stored
            qc - quantum circuit in question
    """

    # Computes the size of register a
    n = len(a)

    loc = adder_inverse(a, b, ob, c, qc, clr, n_c, loc)

    loc = ApplyGate(ob, t, qc, clr, n_c, loc, 'cnot')

    loc = adder(a, b, ob, c, qc, clr, n_c, loc)

    # Changes N register to all zero register if t is set
    for x in range(n):
        if int(n_string[x]) == 1:
            loc = ApplyGate(t, n_reg[n - x - 1], qc, clr, n_c, loc, 'cnot')

    loc = adder_inverse(n_reg, b, ob, c, qc, clr, n_c, loc)

    # Changes N register to all zero register if t is set
    for x in range(n):
        if int(n_string[x]) == 1:
            loc = ApplyGate(t, n_reg[n - x - 1], qc, clr, n_c, loc, 'cnot')

    loc = ApplyGate(ob, ob, qc, clr, n_c, loc, 'x')

    loc = ApplyGate(ob, t, qc, clr, n_c, loc, 'cnot')

    loc = ApplyGate(ob, ob, qc, clr, n_c, loc, 'x')

    loc =  adder(n_reg, b, ob, c, qc, clr, n_c, loc)

    # Computes (a+b) mod N - a and stores it in b
    loc = adder_inverse(a, b, ob, c, qc, clr, n_c, loc)
    
    return loc


def cntr_mult_mod(a, x, zero_b, ob, control_bit, carry_register, a_number, n_reg, n_string, t, qc, clr, n_c, loc):
    """Controlled modular multiplier : if control_bit is set then  b = a * x mod N  else b = x

        Arguments:
            a - one of the factors to multiply
            zero_b - register where the result will be stored
            control_bit - control bit
            carry_register - register where the carries of each addition of function adder will be saved
            a_string - binary string of the number a i.e. if a = 1 a_string = "10"
            n_reg - register where the number N is saved
            n_string - the same as a_string but for number N
            cb_temp_reg - temporary classical register for the overflow bit memory of adder_mod
            t - temporary quantum register
            qc - quantum circuit in question
    """

    # Computes the size of register a
    n = len(a)

    if len(x) != n:
        raise NameError("The x register must have the same value as the a register")

    # Convert n_string to integer

    N_number = int(n_string, 2)

    # The loop where the n stages of adder_mod are executed
    for i in range(n):
        # Calculate the number to add at the ith stage
        a_number_stage = (a_number * 2**i) % N_number

        # Convert the number to the string
        a_number_stage_string = binary_string(a_number_stage, n)

        # Conditionally selects the a or 0
        for j in range(n):
            if int(a_number_stage_string[j]) == 1:
                loc = toffoli(control_bit, x[i], a[n - j - 1], qc, clr, n_c, loc)

        # The adder mod operation
        loc = adder_mod(a, zero_b, ob, carry_register, n_reg, n_string, t, qc, clr, n_c, loc)

        # Conditionally selects the a or 0
        for j in range(n):
            if int(a_number_stage_string[j]) == 1:
                loc = toffoli(control_bit, x[i], a[n - j - 1], qc, clr, n_c, loc)

    for i in range(n):
        # Inverts the control qubit
        loc = ApplyGate(control_bit, control_bit, qc, clr, n_c, loc, 'x')

        # copies the x in y register if c qubit is 0 initially
        loc = toffoli(control_bit, x[i], zero_b[i], qc, clr, n_c, loc)

        # Inverts the control qubit back to it's original state
        loc = ApplyGate(control_bit, control_bit, qc, clr, n_c, loc, 'x')

    return loc


def cntr_mult_mod_inverse(a, x, zero_b, ob, control_bit, carry_register, a_number, n_reg, n_string, t, qc, clr, n_c, loc):
    """Controlled modular multiplier : if control_bit is set then  b = a * x mod N  else b = x

        Arguments:
            a - one of the factors to multiply
            zero_b - register where the result will be stored
            control_bit - control bit
            carry_register - register where the carries of each addition of function adder will be saved
            a_string - binary string of the number a i.e. if a = 1 a_string = "10"
            n_reg - register where the number N is saved
            n_string - the same as a_string but for number N
            cb_temp_reg - temporary classical register for the overflow bit memory of adder_mod
            t - temporary quantum register
            qc - quantum circuit in question
    """

    # Computes the size of register a
    n = len(a)

    # Convert a_string to integer
    N_number = int(n_string, 2)

    for i in range(n):
        # Inverts the control qubit
        loc = ApplyGate(control_bit, control_bit, qc, clr, n_c, loc, 'x')

        # copies the x in y register if c qubit is 0 initially
        loc = toffoli(control_bit, zero_b[i], x[i], qc, clr, n_c, loc)

        # Inverts the control qubit back to it's original state
        loc = ApplyGate(control_bit, control_bit, qc, clr, n_c, loc, 'x')

    # The loop where the n stages of adder_mod are executed
    for i in range(n):
        # Calculate the number to add at the ith stage
        a_number_stage = (a_number * 2**(n - i - 1)) % N_number

        # Convert the number to the string
        a_number_stage_string = binary_string(a_number_stage, n)

        # Conditionally selects the a or 0
        for j in range(n):
            if int(a_number_stage_string[j]) == 1:
                loc = toffoli(control_bit, x[n - i - 1], a[n - j - 1], qc, clr, n_c, loc)

        # The adder mod operation
        loc = adder_mod_inverse(a, zero_b, ob, carry_register, n_reg, n_string, t, qc, clr, n_c, loc)

        # Conditionally selects the a or 0
        for j in range(n):
            if int(a_number_stage_string[j]) == 1:
                loc = toffoli(control_bit, x[n - i - 1], a[n - j - 1], qc, clr, n_c, loc)


        
    return loc


def mod_exp(a, ones_state, zero_b, ob, x_input, carry_register, a_number, n_reg, n_string, t, qc, clr, n_c, loc):

    # Computes the size of register a
    n = len(x_input)

    n_number = int(n_string, 2)

    for i in range(n):

        # The ith stage of controlled modular multiplier
        a_number_stage = (a_number ** 2 ** i) % n_number

        # Calculates inverse mod of a_number_stage mod N
        a_number_stage_inverse = (modInverse(a_number_stage, n_number)) % n_number

        if (i % 2) == 0:
            cntr_mult_mod(a, ones_state, zero_b, ob, x_input[i], carry_register, a_number_stage, n_reg,
                          n_string, t, qc, clr, n_c, loc)

            # The ith stage of controlled modular multiplier inverse
            cntr_mult_mod_inverse(a, zero_b, ones_state, ob, x_input[i], carry_register, a_number_stage_inverse,
                                  n_reg, n_string, t, qc, clr, n_c, loc)

        else:

            cntr_mult_mod(a, zero_b, ones_state, ob, x_input[i], carry_register, a_number_stage, n_reg,
                          n_string, t, qc, clr, n_c, loc)

            # The ith stage of controlled modular multiplier inverse
            cntr_mult_mod_inverse(a, ones_state, zero_b, ob, x_input[i], carry_register,
                                  a_number_stage_inverse, n_reg, n_string, t, qc, clr, n_c, loc)

    return loc

def InverseQFT(n_quibits, qc, clr, n_c, loc):
    n = len(n_quibits)
    for i in range(n):
        loc = ApplyGate(n-1-i, n-1-i, qc, clr, n_c, loc, 'h')
        for j in range(n-1-i):
            print(n-1-i, n-2-j-i)
            loc = ApplyGate(n-1-i, n-2-j-i, qc, clr, n_c, loc, 'cphase', 2 * np.pi / (2 ** ((j - i) + 2)))

    return loc
