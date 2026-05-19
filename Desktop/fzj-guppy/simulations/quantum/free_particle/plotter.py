# plotter.py
import numpy as np
import matplotlib.pyplot as plt

class WavePlotter:
    def __init__(self, num_sites: int):
        self.num_sites = num_sites
        self.x = np.arange(num_sites)

        plt.ion()
        self.fig, self.ax = plt.subplots()

        # Start with zero wave
        self.line, = self.ax.plot(self.x, np.zeros(num_sites), '-o', linewidth=2, markersize=6)

        self.ax.set_xlim(-0.5, num_sites - 0.5)
        self.ax.set_ylim(-1.0, 1.0)  # gets updated dynamically
        self.ax.set_xlabel("Lattice Position")
        self.ax.set_ylabel("Amplitude (Re[ψ])")
        self.ax.set_title("Wavefunction Evolution (Smooth Amplitude Curve)")

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def update(self, site_amps: np.ndarray):
        """
        site_amps: complex array of length num_sites
        (output of QuantumSimulator.site_amps_onehot(statevector))
        """

        # Display **signed real amplitude**
        displayed = np.real(site_amps)

        # Update line data
        self.line.set_ydata(displayed)

        # Rescale axis so wave stays visible
        max_amp = max(1e-9, np.max(np.abs(displayed)))
        self.ax.set_ylim(-1.1 * max_amp, 1.1 * max_amp)

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()