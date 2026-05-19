import numpy as np
from qiskit import QuantumCircuit as QC


class TimeEvolution:
    def __init__(self, qc: QC, psi: np.ndarray, b: float, a: float):
        self.psi = psi
        self.U = np.matrix([[1, 0, 0, 0],
                             [0, b, a, 0],
                             [0, a, b, 1],
                             [0, 0, 0, 1]])
        self.qc = qc

    def time_transform(self):
        theta = np.arccos(self.U[1,1].real)
        # theta = np.pi/2  # Example fixed angle for demonstration with iSwap
        n = range(self.qc.num_qubits - 1)
        N = self.qc.num_qubits
        for i in n:
            self.qc.rxx(-theta, i, i+1)
            self.qc.ryy(-theta, i, i+1)

        # <-- this is the boundary condition: wrap the chain
        self.qc.rxx(-theta, N-1, 0)
        self.qc.ryy(-theta, N-1, 0)
        return self.qc

    def apply_time_evolution(self, time_steps: int) -> QC:
        for i in range(time_steps):
            self.qc = self.time_transform()
        return self.qc