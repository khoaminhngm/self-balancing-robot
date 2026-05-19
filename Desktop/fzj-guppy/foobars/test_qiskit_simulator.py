from qiskit import QuantumCircuit
from qiskit_aer import QasmSimulator
import numpy as np

qc = QuantumCircuit(2, 2) # Create a quantum circuit with 2 qubits and 2 classical bits
qc.h(0)             # Put qubit 0 into superposition
qc.cx(0, 1)         # Entangle qubit 1 with qubit 0
qc.measure([0, 1], [0, 1])  # Measure both qubits

backend = QasmSimulator(method='statevector')
job = backend.run(qc, shots=np.pow(2,25))
result = job.result()

counts = result.get_counts()
print(counts)
