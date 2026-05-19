from qiskit_aer import StatevectorSimulator
import numpy as np
from itertools import product
from qiskit import QuantumCircuit as QC
from qiskit.quantum_info import Statevector

class QuantumSimulator:
    def __init__(self, qc: QC, shots: int, psi: np.ndarray):
        self.qc = qc
        self.psi = psi / np.linalg.norm(psi)
        self.shots = shots
        self.backend = StatevectorSimulator(precision='double')
        # self.backend.clear_options()  # uncomment if you previously set options

    def _statevector_circuit(self):
        qc_sv = self.qc.copy()
        qc_sv.remove_final_measurements()
        qc_sv.save_statevector()
        return qc_sv

    def run(self):
        """Return full complex statevector (no backend, no measurement)."""
        # Work on a copy so original circuit isn't destroyed
        qc_sv = self.qc.copy()

        # Remove all measurement instructions (new Qiskit syntax)
        qc_sv.remove_final_measurements(inplace=True)
        qc_sv.data = [inst for inst in qc_sv.data if inst.operation.name != "measure"]

        # Convert directly to statevector (no transpile, no backend)
        sv = Statevector.from_instruction(qc_sv)
        return np.asarray(sv.data)

    def site_amps_onehot(self, statevector: np.ndarray, endian: str = "little") -> np.ndarray:
        """
        Convert full statevector → site amplitudes assuming one-hot position encoding.
        Returns a vector of length num_sites representing ψ(x).
        """
        n = self.qc.num_qubits
        site_amps = np.zeros(n, dtype=complex)

        for j in range(n):
            idx = (1 << j) if endian == "little" else (1 << (n - 1 - j))
            site_amps[j] = statevector[idx]

        return site_amps

    def post_process(self, statevector: np.ndarray, tol: float = 1e-12):
        """
        Return only basis states with non-zero amplitudes.
        Useful for debugging only.
        """
        n = self.qc.num_qubits
        states = [''.join(bits) for bits in product('01', repeat=n)]
        amps = np.asarray(statevector)

        mask = np.abs(amps) > tol
        return np.array(states)[mask], amps[mask]

    def dispersion_factor(self, final_site_amps: np.ndarray):
        """
        Compare initial ψ(x) to final ψ(x):
        0 = unchanged
        1 = maximally spread.
        """
        p_init = np.abs(self.psi)**2
        p_init /= p_init.sum()

        p_final = np.abs(final_site_amps)**2
        if p_final.sum() == 0:
            return 1.0
        p_final /= p_final.sum()

        num = np.sum(p_init * p_final)
        den = np.sqrt(np.sum(p_init**2) * np.sum(p_final**2))

        similarity = (num / den) if den != 0 else 0.0
        return 1.0 - similarity