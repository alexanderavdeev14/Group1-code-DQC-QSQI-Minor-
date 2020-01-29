import os
from getpass import getpass
from quantuminspire.credentials import load_account, get_token_authentication, get_basic_authentication

from qiskit.tools.visualization import circuit_drawer, plot_histogram, plot_state_city
from qiskit import BasicAer
from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import execute

from quantuminspire.qiskit import QI
from quantuminspire.credentials import save_account

import math
import numpy as np



# Function establishes entanglment between a control bit and a GHZ qubit on another computer using a cat-entangler
def single_entangler(c, loc_indexc, loc_indext, cb, qc, n_c):
    """ cat-entangler of control bit c and n qubits GHZ state """
    # c - control qubit (int)
    # loc_indexc - computer that the control qubit is located on (int)
    # loc_indext - computer that the target qubit is located on (int)
    # cb - classical bit where the outcome of first bit of the GHZ will be stored (classical register) 
    # qc - quantum circuit in question 
    # n_c - number of qubits per computer (squared integer)

    # Find the qubit number of the GHZ state
    ghz_c = n_c * (loc_indexc + 1) + loc_indexc
    ghz_t = n_c * (loc_indext + 1) + loc_indext

    # Harmard gate  at the first line of GHZ
    qc.h(ghz_c)
    # CNOT gate between GHZ states
    qc.cx(ghz_c, ghz_t)
    # CNOT of control bit and first ghz bit
    qc.cx(c, ghz_c)
    # Measure the result of the first line of GHZ state
    qc.measure(ghz_c, cb)
    qc.barrier()
    # Controlled x gate based on measured bit of the GHZ
    qc.x(ghz_c).c_if(cb, 1)
    qc.x(ghz_t).c_if(cb, 1)


# Function terminates entanglement between a control qubit and a GHZ qubit from another computer,
# While resetting the control bit to it's original state
def single_disentangler(c, loc_indext, cb, qc, n_c):
    # c - control bit that is reset to original state (int)
    # qc - quantum circuit (int)
    # cb - classical bit where the outcome of the measurement will be stored
    # loc_indext - computer that the target qubit is located on
    # n_c - number of qubits per computer

    # Find the qubit number of the GHZ state
    ghz_t = n_c * (loc_indext + 1) + loc_indext
    qc.h(ghz_t)
    qc.measure(ghz_t, cb)
    qc.barrier()
    qc.z(c).c_if(cb, 1)
    qc.x(ghz_t).c_if(cb, 1)


# Function generates a square matrix with indices counting from 0 to n-1
def sqtop(n):
    # n - number of qubits in the matrix  (squared int)
    if math.sqrt(n) % 1 == 0:
        rootn = int((math.sqrt(n)))
        numb = np.linspace(0, n - 1, n)
        shape = (rootn, rootn)
        top = numb.reshape(shape)
        top = top.astype(int)
        return top
    else:
        raise CircuitError('n_c is not a squared integer')
        return


# Swap entries in the location matrix
def swloc(loc1_row, loc1_col, loc2_row, loc2_col, loc):
    if ((loc1_row == loc2_row) and (loc1_col == loc2_col)):
        raise CircuitError("Unable to swap identical locations")
        return
    else:
        a = loc[loc1_row, loc1_col]
        b = loc[loc2_row, loc2_col]
        loc[loc1_row, loc1_col] = b
        loc[loc2_row, loc2_col] = a
        return loc


# returns the row and column location of a qubit
def FindLoc(qbit, loc):
    # qbit - qubit that is located (int)
    # loc - location matrix that contains the target qubit

    a = np.where(loc == qbit)
    # check for duplicate value
    if np.size(a) > 2:
        display(np.size(a))
        raise CircuitError('Error: duplicate qubit location')
    else:
        row_q = a[0][0]
        column_q = a[1][0]
        return row_q, column_q


# Function finds the location of qubits and moves them such that they are nearest neighbours
# returns the new location matrix after moving the qubits
def NNMove(c, t, n_c, qc, loc, loc_index=0):
    # c - first qubit (int)
    # t - second qubit (int)
    # n_c - number of qubits per computer (squared int)
    # qc - quantum circuit
    # loc - location matrix that contains the two qubits
    # loc_index - the index of the set of location matrix corresponding to the matrix 'loc'

    # find positions of the qbits in matrix indices i and j
    c_i, c_j = FindLoc(c, loc)
    t_i, t_j = FindLoc(t, loc)

    # find the differences between the indices
    difi = c_i - t_i
    difj = c_j - t_j

    # Give an error for same bits
    if difi == 0 and difj == 0:
        raise CircuitError('Bits are the same')
        return loc

    # move the bits untill they are neighbours.
    else:
        while not ((difi == 0 and abs(difj) == 1) or (abs(difi) == 1 and difj == 0)):

            # The rows

            # If the difference between control and target row (difi) is positive,
            # move the target bit down
            if (difi > 1) or ((difi == 1) and abs(difj) > 0):
                loc = swapqubit(t, loc, 'down', qc, n_c, loc_index)
                difi = difi - 1

            # If the difference between control and target row (difi) is negative,
            # move the control bit down
            elif (difi < -1) or ((difi == -1) and abs(difj) > 0):
                loc = swapqubit(c, loc, 'down', qc, n_c, loc_index)
                difi = difi + 1

            # the colums

            # If the difference between control and target column (difj) is positive,
            # move the target bit right
            elif (difj > 1) or ((difj == 1) and abs(difi) > 0):
                loc = swapqubit(t, loc, 'right', qc, n_c, loc_index)
                difj = difj - 1

            # If the difference between control and target column (difj) is negative,
            # move the conrol bit right
            elif (difj < -1) or ((difj == -1) and abs(difi) > 0):
                loc = swapqubit(c, loc, 'right', qc, n_c, loc_index)
                difj = difj + 1
            else:
                raise CircuitError('location difference conditions not met')

        return loc


# Move a qubit in a specified direction
# Function returns the new location matrix after the swap gate is applied
def swapqubit(t, loc, direction, qc, n_c, loc_index=0):
    # t - target qubit that is moved (int)
    # loc - locaton matrix that contains the target qubit
    # direction - direction that the qubit is moved in ('up', 'down', 'left' or 'right')
    # qc - quantum circuit
    # n_c - number of qubits per computer (squared integer)
    # loc_index - the index of the set of location matrix corresponding to the matrix 'loc'
    
    #initialize global counter
    global NormalCounter
    
    # determine qubit location in its matrix
    loc_i, loc_j = FindLoc(t, loc)

    n = int(math.sqrt(n_c))
    # Calculate target qubit in the square topology that needs to be swapped
    t_swap = loc_index * n_c + loc_i * n + loc_j + loc_index

    # t_swap1 is the qubit in the topology that t_swap is swapped with
    if direction == 'up':
        t_swap1 = loc_index * n_c + (loc_i - 1) * n + loc_j + loc_index
        qc.swap(t_swap, t_swap1)
        loc = swloc(loc_i, loc_j, loc_i - 1, loc_j, loc)
    elif direction == 'down':
        t_swap1 = loc_index * n_c + (loc_i + 1) * n + loc_j + loc_index
        qc.swap(t_swap, t_swap1)
        loc = swloc(loc_i, loc_j, loc_i + 1, loc_j, loc)
    elif direction == 'right':
        t_swap1 = loc_index * n_c + loc_i * n + loc_j + 1 + loc_index
        qc.swap(t_swap, t_swap1)
        loc = swloc(loc_i, loc_j, loc_i, loc_j + 1, loc)
    elif direction == 'left':
        t_swap1 = loc_index * n_c + loc_i * n + loc_j - 1 + loc_index
        qc.swap(t_swap, t_swap1)
        loc = swloc(loc_i, loc_j, loc_i, loc_j - 1, loc)
    else:
        raise CircuitError("ERROR nonvalid direction")
    # return with the new location matrix
    NormalCounter = NormalCounter + 3
    return loc


# Move a target bit to a destination of choice #MoveQbitTo(t,dest,qc,loc)
# Function returns the new location matrix after moving the target qubit
def MoveQbitTo(t, dest_i, dest_j, qc, n_c, loc, loc_index=0):
    # t - target qubit that is moved (int)
    # dest_i - row that the target qubit is moved to 
    # dest_j - column that the target qubit is moved to
    # qc - quantum circuit
    # n_c - number of qubits per computer (squared integer)
    # loc - locaton matrix that contains the target qubit 
    # loc_index - the index of the set of location matrix corresponding to the matrix 'loc'

    # find target location first
    t_i, t_j = FindLoc(t, loc)

    dif_i = dest_i - t_i
    dif_j = dest_j - t_j

    # move target bit until it is has arrived at the destination
    while (dif_i != 0) or (dif_j != 0):
        if dif_i > 0:
            loc = swapqubit(t, loc, 'down', qc, n_c, loc_index)
            dif_i = dif_i - 1
        elif dif_i < 0:
            loc = swapqubit(t, loc, 'up', qc, n_c, loc_index)
            dif_i = dif_i + 1
        elif dif_j > 0:
            loc = swapqubit(t, loc, 'right', qc, n_c, loc_index)
            dif_j = dif_j - 1
        elif dif_j < 0:
            loc = swapqubit(t, loc, 'left', qc, n_c, loc_index)
            dif_j = dif_j + 1
        else:
            break

    return loc


# Function generates a qubit topology
# returns a quantum circuit (qc).
# returns location matrices indicating where the algorithm qubits are currently located in the topology (loc)
# returns a set of classical registers for the ghz qubit measurement (clr)
# The size of the quantum circuit that is returned is larger than n_tot due to the extra ghz
# qubits that are generated per computer
def GenerateQCircuit(n_tot, n_c, Cname):
    # n_tot - total number of qubits needed in the  (int) 
    # n_c - number of qubits per computer (squared int) (e.g. 4 or 9 or 16 etc.)
    # Cname - name of the squatum circuit (e.g. 'QFT')

    # Check if n_c is a squared integer
    if math.sqrt(n_c) % 1 != 0:
        raise CircuitError('n_c is not a squared integer')

    # Calculate the number of computers that is necessary to accomodate for all the qubits that are used
    mod = n_tot % n_c
    if mod > 0:
        counter = (n_tot - mod) / n_c + 1
    else:
        counter = n_tot / n_c

    counter = int(counter)

    # Define a set of location matrices
    loc = []
    loc_new = sqtop(n_c)
    loc.append(loc_new)

    # Set the first quantum- and classical register
    Qregister = QuantumRegister(n_c, 'q0')
    GHZregister = QuantumRegister(1, 'qghz0')
    # initialize the quantum circuit
    qc = QuantumCircuit(Qregister, GHZregister, name=Cname)

    # Add quantumregisters and classical registers to the initialized quantum circuit
    if counter > 0:
        for x in range(1, counter):
            # generate new quantum registers
            qname = "q{}".format(x)
            Qregister = QuantumRegister(n_c, qname)
            ghzname = "qghz{}".format(x)
            GHZregister = QuantumRegister(1, ghzname)
            cname = "c{}".format(x)

            # Add the new registers to the already existing quantum circuit
            qc.add_register(Qregister)
            qc.add_register(GHZregister)

            # append the vector of location matrices
            loc_new = sqtop(n_c) + n_c * x
            loc.append(loc_new)

            # initialize a set of classical registers
    clr = []
    for n in range(0, counter):
        # add new classical register to the quantum circuit and append the tot the set of classical registers
        clr.append(ClassicalRegister(1))
        qc.add_register(clr[n])

    return qc, loc, clr


# Function applies a gate to a qubit within an algorithm. It searches where in the topology the qubit is located.
# For controlled gates, the function will then make sure the control and targer qubit are nearest neighbours.
# If the control and target qubit are on different computer, the function will set up entanglement between these computers,
# and after the gate is applied, this entanglement is terminated and the control qubit is reset to its original state.
# The function returns the new location matrices after the gate has been applied.
def ApplyGate(c, t, qc, clr, n_c, loc, gate, phase=0):
    # c - control qubit (int)
    # t - target qubit (int)
    # if a single (non-controlled) gate needs to be applied, set te control and target qubit to be equal (c=t)
    # qc - quantum circuit
    # clr - set of classical registers
    # n_c - number of qubits per computer (squared integer)
    # loc - set of location matrices correspoding to the quatum circuit
    # gate - gate that is applied to the circuit (available gates: 'h', 'x', 'y', 'z', 't', 'tdg', 'cnot', 'cphase')
    # phase - phase that is applied during a cphase gate (default phase=0) 
    
    global NormalCounter
    global ECounter
    global OneCounter
    
    # calculate the number of qubits per row in a computer
    n_row = math.sqrt(n_c)

    if c < 0 or t < 0:
        raise CircuitError("non-valid target or control qubit")

    elif n_row % 1 != 0:
        raise CircuitError("n_c is not a squared integer")

    # single (non-controlled) gate if control(c) and target(t) are the same
    elif c == t:
        # find which computer the qubit is located on
        mod_c = c % n_c
        comp_no_c = int((c - mod_c) / n_c)

        # Find the location of the target qubit
        cloc_i, cloc_j = FindLoc(c, loc[comp_no_c])

        # Determine the qubit in the circuit that the gate needs to be applied to
        c_index = int(n_c * comp_no_c + comp_no_c + n_row * cloc_i + cloc_j)

        # We can apply the gate now
        if gate == 'x':
            qc.x(c_index)
        elif gate == 'y':
            qc.y(c_index)
        elif gate == 'z':
            qc.z(c_index)
        elif gate == 'h':
            qc.h(c_index)
        elif gate == 'tdg':
            qc.tdg(c_index)
        elif gate == 't':
            qc.t(c_index)
        else:
            raise GateError('Non-existent gate')
        
        OneCounter = OneCounter + 1
        
        return loc

    # find location of control and target
    # number of the computer that the control qubit is located on.
    mod_c = c % n_c
    comp_no_c = int((c - mod_c) / n_c)
    # number of the computer that the target qubit is located on.
    mod_t = t % n_c
    comp_no_t = int((t - mod_t) / n_c)

    # find the qubit locations withing the computers that they are located in
    cloc_i, cloc_j = FindLoc(c, loc[comp_no_c])
    tloc_i, tloc_j = FindLoc(t, loc[comp_no_t])

    # Check if qubits are located on the same computer
    if comp_no_c == comp_no_t:
        # move control and target qubit such that they are in neares neighbour position
        loc[comp_no_c] = NNMove(c, t, n_c, qc, loc[comp_no_c], comp_no_c)

        # Find the new locations of the control and target qubit
        cloc_i, cloc_j = FindLoc(c, loc[comp_no_c])
        tloc_i, tloc_j = FindLoc(t, loc[comp_no_t])

        # from the new location matrices, determine the qubit in the circuit that the gate needs to be applied to
        c_index = int(n_c * comp_no_c + comp_no_c + n_row * cloc_i + cloc_j)
        t_index = int(n_c * comp_no_t + comp_no_c + n_row * tloc_i + tloc_j)

        # We now have ensured that the control and target qubits are nearest neighbours,
        # so we can apply the gate
        if gate == 'cnot':
            qc.cx(c_index, t_index)
        elif gate == 'cphase':
            qc.cu1(phase, c_index, t_index)
        else:
            raise GateError('Non-existent gate')
            
        NormalCounter = NormalCounter + 1

    else:
        # move qubits to the locations within their computer where they can interact with
        # First, calculate the destination coordinates
        dest_i = math.sqrt(n_c) - 1
        dest_j = math.sqrt(n_c) - 1

        # Use the MoveQbitTo function to move the bit to the right position
        loc[comp_no_c] = MoveQbitTo(c, dest_i, dest_j, qc, n_c, loc[comp_no_c], comp_no_c)
        loc[comp_no_t] = MoveQbitTo(t, dest_i, dest_j, qc, n_c, loc[comp_no_t], comp_no_t)

        # Find the new locationx of the control and target qubits
        tloc_i, tloc_j = FindLoc(t, loc[comp_no_t])
        cloc_i, cloc_j = FindLoc(c, loc[comp_no_c])

        # from the new location matrices, determine the qubit in the circuit that the gate needs to be applied to
        t_index = int(n_c * comp_no_t + comp_no_t + n_row * tloc_i + tloc_j)
        c_index = int(n_c * comp_no_c + comp_no_c + n_row * cloc_i + cloc_j)



        # Next, establish an entanglement between the GHZ qubits of the two computers
        single_entangler(c_index, comp_no_c, comp_no_t, clr[comp_no_c], qc, n_c)
        
        # The control qubit is now the ghz qubit on the computer of the target qubit
        c_ghz = (comp_no_t + 1) * n_c + comp_no_t

        # Apply gate
        # Now the target qubit can interact with it's ghz qubit
        if gate == 'cnot':
            qc.cx(c_ghz, t_index)
        elif gate == 'cphase':
            qc.cu1(phase, c_ghz, t_index)
        else:
            raise GateError('Non-existent gate')

        # Disentagle the computers after the gate has been applied

        single_disentangler(c_index, comp_no_t, clr[comp_no_t], qc, n_c)
        
        NormalCounter = NormalCounter + 1
        ECounter = ECounter + 1

    return loc


# Function returns the values of a measured state adjusted for the location matrix
def ReturnState(loc, state, n_c):
    # loc - set of location matrices
    # state_set - set of measured states, not adjusted for location
    # state_adjusted - set of measured states, adjusted for location
    # n_c - number of qubits per computer
    
    
    state_new = state.replace(" ", "")

    state_set = list(state_new)

    state_int = np.zeros(len(state_set))

    for ii in range(0, len(state_set)):
        state_int[ii] = int(state_set[ii])

    state_set = state_int
    new_state = state_set[::-1]
    print(new_state)
    
    state_adjusted = np.zeros(len(loc) * n_c)

    # Find the number of qubits in one row
    n_row = int(math.sqrt(n_c))

    for x in range(0, len(loc)):
        for y in range(0, n_row):
            for z in range(0, n_row):
                # find the qubit inhabiting this location
                qubit_here = loc[x][y][z]

                # find corresponding index in the state_set vector
                #value = state_set[z + y * n_row + x * (n_c + 1)]
                value = new_state[z + y*n_row + x*n_c+len(loc)]

                # Write this value to the corresponding adjusted state vector index
                state_adjusted[qubit_here] = int(value)

    return state_adjusted

def SetCounters():
    global NormalCounter
    NormalCounter = 0
    global ECounter
    ECounter = 0
    global OneCounter
    OneCounter = 0
    return

def GetCounters():
    global NormalCounter
    global ECounter
    global OneCounter
    return NormalCounter, ECounter, OneCounter