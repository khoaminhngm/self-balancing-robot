from state_prep import State_Preparator
from time_evol import TimeEvolution
from measure import Measurement
from quantum_simulator import QuantumSimulator
from plotter import WavePlotter

import numpy as np
from qiskit import QuantumCircuit

# 1. Define simulation parameters
# 1.1 Define amplitudes of discretized Schrödinger Wave here
psi = np.array([0.1, 0.1, 0.1, 0.7, 0.9, 0.7, 0.1, 0.1, 0.1], dtype=complex) # complex wavefunction
psi /= np.linalg.norm(psi)  # Normalize the wavefunction
print(psi)

# 1.2 Define unitary time transformation matrix here:
# a -> probability that the particle moves to the left/right, equal prob
# b -> probability that the particle stays in the same position
b = 0.8660254037844386
a = 0.49999999999999994j

# 1.3 Number of time steps:
time_steps = 2

# 1.4 Number of measurement shots:
shots = np.pow(2, 13)

# 2. Init Instances of Modules & Quantum Circuit
qc = QuantumCircuit(len(psi), len(psi))

state_preparator = State_Preparator(psi, qc)
time_evolution = TimeEvolution(qc, psi, b, a)
measurement = Measurement(qc)
quantum_simulator = QuantumSimulator(qc, shots, psi)
# plotter = WavePlotter(num_sites=len(psi))


for t in range(time_steps):
    # 3. Prepare Initial State
    qc = state_preparator.prepare()
    # 4. Apply Time Evolution
    qc = time_evolution.apply_time_evolution(t)
    # 5. Measurement
    qc = measurement.measure_all()
    # 6. Run Simulation
    sv = quantum_simulator.run()
    # 7. Post-process Results
    site_amps = quantum_simulator.site_amps_onehot(sv)         # complex per-site amps (one-hot lattice)
    disp = quantum_simulator.dispersion_factor(site_amps)  
    # Plot the results at each time step
    # plotter.update(site_amps)   


# # 3. Prepare Initial State
# qc = state_preparator.prepare()
# print("State preparation complete.")

# # 4. Apply Time Evolution
# qc = time_evolution.apply_time_evolution()
# print("Time evolution complete.")

# # 5. Measurement
# qc = measurement.measure_all()
# print("Measurement setup complete.")

# # 6. Run Simulation
# counts = quantum_simulator.run()

# # 7. Post-process Results
# final_amps = quantum_simulator.post_process(counts)
# print("Post-processing complete.")
# print(f"\n\nFinal Amps: {final_amps}")

# # 8. Compute Dispersion Factor
# dispersion = quantum_simulator.dispersion_factor(final_amps)
# print(f"\n\n Dispersion Factor: {dispersion}")