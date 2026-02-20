# plot_energy_pressure_vs_step.py
# Reads a Monte Carlo log and plots E*/N and P* vs step

import re
import numpy as np
import matplotlib.pyplot as plt

# ---- set your file path here ----
fname = "indep-1.0.1.Monte_Carlo.log"
# Example (your uploaded file path in this session):
# fname = "/mnt/data/indep-1.0.1.Monte_Carlo.log"

# Regex to capture: Step, E*/N, P*
pat = re.compile(
    r"Step:\s*(\d+).*?E\*/N:\s*([-\d\.Ee+]+).*?P\*:\s*([-\d\.Ee+]+)"
)

steps = []
e_per_n = []
p_star = []

with open(fname, "r", errors="ignore") as f:
    for line in f:
        m = pat.search(line)
        if m:
            steps.append(int(m.group(1)))
            e_per_n.append(float(m.group(2)))
            p_star.append(float(m.group(3)))

if not steps:
    raise RuntimeError("No data found. Check the file format or the regex pattern.")

steps = np.array(steps)
e_per_n = np.array(e_per_n)
p_star = np.array(p_star)

# ---- plots ----
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(9, 6))

ax[0].plot(steps, e_per_n, linewidth=1.5)
ax[0].set_ylabel("E*/N")
ax[0].grid(True)

ax[1].plot(steps, p_star, linewidth=1.5)
ax[1].set_xlabel("Step")
ax[1].set_ylabel("P*")
ax[1].grid(True)

plt.tight_layout()
plt.savefig("energy_pressure_vs_step.png", dpi=300)