from qiskit import QuantumCircuit as QC
from qiskit import transpile
from qiskit.circuit.library import iSwapGate
import numpy as np
from itertools import product
from qiskit_aer import QasmSimulator


# 1.1 Define amplitudes of discretized Schrödinger Wave here
psi = np.array([0.1, 0.1, 0.1 ,0.1, 0.1, 0.1 ,0.7, 0.9,0.7 ,0.1, 0.1, 0.1 ,0.1, 0.1, 0.1, 0.1, 0.1, 0.1], dtype=complex) # complex wavefunction
psi /= np.linalg.norm(psi)  # Normalize the wavefunction
print(psi)

# 1.2 Define unitary time transformation matrix here:
# a -> probability that the particle moves to the left/right, equal prob
# b -> probability that the particle stays in the same position
b = 0.8660254037844386
a = 0.49999999999999994j
U = np.matrix([[1, 0, 0, 0],
               [0, b, a, 0],
               [0, a, b, 1],
               [0, 0, 0, 1]])

# 1.3 Number of time steps:
time_steps = np.pow(10, 4)

# 1.4 Number of measurement shots:
shots = np.pow(2, 13)

# 2. State preparation algo
def flip_qubits(qc : QC, qubits : list[int]) -> QC:
    for q in qubits:
        print(f"Flipping q{q}")
        qc.x(q)
    return qc

def state_preparation(qc : QC) -> QC:
    controls = [0]
    
    # Init first qubit
    theta1 = 2*np.arcsin(psi[0]).real
    qc.ry(theta1, 0)

    # # Prepare second qubit
    theta2 = 2*np.arcsin(psi[1]/np.sqrt(1-psi[0]**2)).real
    qc.cry(theta2, 0, 1, ctrl_state=0)

    # Subsequent qubits
    for i in range(1, qc.num_qubits - 1): # i+1 -> target
        controls.append(i)
        print("Controls:", controls)

        theta = 2*np.arcsin(
            psi[i+1]/np.sqrt(1-sum([abs(psi[j])**2 for j in range(i+1)]))
            ).real
        
        print("Theta-loop:", theta)

        theta2 = 2*np.arcsin(psi[1]/np.sqrt(1-psi[0]**2)).real
        print("Theta 2-manual:", theta2)
        theta3 = 2*np.arcsin(psi[2]/np.sqrt(1-psi[0]**2-psi[1]**2)).real
        print("Theta 3-manual:", theta3)

        flip_qubits(qc, controls)
        qc.mcry(theta, controls, qc.qubits[i+1])
        print(f"Applied mcry with theta={theta} on qubit {i+1} with controls {controls}")
        flip_qubits(qc, controls)

    return qc


# 3. Create & prepare state circuit
qc = QC(len(psi), len(psi))
print("Initial circuit: |000>")
qc = state_preparation(qc)


# 4. Define Unitary transformation matrix
def time_transform(circuit: QC, U: np.matrix):
    theta = np.arccos(U[1,1].real)
    # theta = np.pi/2  # Example fixed angle for demonstration
    n = range(circuit.num_qubits - 1)
    for i in n:
        circuit.rxx(-theta, circuit.qubits[i], circuit.qubits[i+1])
        circuit.ryy(-theta, circuit.qubits[i], circuit.qubits[i+1])
    return circuit

# 5. Apply time transformation
for i in range(time_steps):
    qc = time_transform(qc, U)

# 6. Measurement
measure_list = []
for i in range(0, qc.num_qubits):
    measure_list.append(i)
qc.measure(measure_list, measure_list)
print(measure_list)

# 7. Transpile and simulate on Qasm Simulator
# Note.. as reference:
# My Macbook pro M1 Max 2023 can simulate up to 12 qubits on this circuit before it gets super hot... wouldn't risk more
backend = QasmSimulator(method='statevector')
job = backend.run(qc, shots=shots)
result = job.result()

counts = result.get_counts()
print(counts)

# 8. Process results
# --- total number of shots
shots = sum(counts.values())

# --- all 5-qubit basis states (2^5 = 32)
states = [''.join(bits) for bits in product('01', repeat=qc.num_qubits)]

# --- convert counts → probabilities
pvec = np.array([counts.get(s, 0) / shots for s in states])

# --- make a matching NumPy array of state labels
states_arr = np.array(states)

# --- eliminate zero components
mask = pvec > 0
pvec_nonzero = pvec[mask]
states_nonzero = states_arr[mask]

# --- (optional) compute amplitudes from probabilities
amps_nonzero = np.sqrt(pvec_nonzero)

# --- pretty print results
print("Non-zero states:", states_nonzero)
print("Probabilities:", pvec_nonzero)
print("Amplitudes:", amps_nonzero)

# --- optional: combine into a dictionary
filtered = dict(zip(states_nonzero, pvec_nonzero))
filtered_amps = dict(zip(states_nonzero, amps_nonzero))
print("Dictionary (probabilities):", filtered)
print("Dictionary (amplitudes):", filtered_amps)


# 9. Calculate dispersion factor 
# a physical Schrödinger free particle wave packet should dispers over time
# test physical accuracy by calculating dispersion factor here
def dispersion_factor(initial_psi: np.ndarray, final_amps: np.ndarray) -> float:
    """
    Compare initial and final amplitude distributions to quantify dispersion.
    Returns a number in [0, 1]:
        0 = identical distributions (no spread)
        1 = fully different / maximally dispersed
    """
    # --- Normalize both wavefunctions
    p_init = np.abs(initial_psi)**2
    p_init /= p_init.sum()
    
    p_final = np.abs(final_amps)**2
    p_final /= p_final.sum()
    
    # --- Compute cosine similarity of the two probability vectors
    numerator = np.sum(p_init * p_final)
    denom = np.sqrt(np.sum(p_init**2) * np.sum(p_final**2))
    similarity = numerator / denom if denom != 0 else 0.0
    
    # --- Dispersion factor (1 - similarity)
    dispersion = 1 - similarity
    return dispersion

print(f"Dispersion factor: {dispersion_factor(psi, amps_nonzero)}")