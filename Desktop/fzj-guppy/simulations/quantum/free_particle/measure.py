from qiskit import QuantumCircuit as QC
import numpy as np

class Measurement:
    def __init__(self, qc: QC):
        self.qc = qc
        self.measure_list = []

    def measure_all(self):
        """Measure all qubits in the circuit."""
        for i in range(0, self.qc.num_qubits):
            self.measure_list.append(i)
        self.qc.measure(self.measure_list, self.measure_list)