from qiskit import QuantumCircuit as QC
import numpy as np

class State_Preparator:
    def __init__(self, psi: np.array, qc: QC):
        self.psi = psi / np.linalg.norm(psi) # normalize
        self.qc = qc

    def flip_qubits(self, qubits: list[int]) -> QC:
        for q in qubits:
            # print(f"Flipping q{q}")
            self.qc.x(q)
        return self.qc
    
    def prepare(self):
        controls = [0]
    
        # Init first qubit
        theta1 = 2*np.arcsin(self.psi[0]).real
        self.qc.ry(theta1, 0)

        # # Prepare second qubit
        theta2 = 2*np.arcsin(self.psi[1]/np.sqrt(1-self.psi[0]**2)).real
        self.qc.cry(theta2, 0, 1, ctrl_state=0)

        # Subsequent qubits
        for i in range(1, self.qc.num_qubits - 1): # i+1 -> target
            controls.append(i)
            # print("Controls:", controls)

            theta = 2*np.arcsin(
                self.psi[i+1]/np.sqrt(1-sum([abs(self.psi[j])**2 for j in range(i+1)]))
                ).real
            
            # print("Theta-loop:", theta)

            theta2 = 2*np.arcsin(self.psi[1]/np.sqrt(1-self.psi[0]**2)).real
            # print("Theta 2-manual:", theta2)
            theta3 = 2*np.arcsin(self.psi[2]/np.sqrt(1-self.psi[0]**2-self.psi[1]**2)).real
            # print("Theta 3-manual:", theta3)

            self.flip_qubits(controls)
            self.qc.mcry(theta, controls, self.qc.qubits[i+1])
            # print(f"Applied mcry with theta={theta} on qubit {i+1} with controls {controls}")
            self.flip_qubits(controls)

        return self.qc