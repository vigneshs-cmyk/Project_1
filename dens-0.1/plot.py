# plot_energy_pressure_vs_step.py
# Reads a Monte Carlo log and plots E*/N and P* vs step

import re
import numpy as np
import matplotlib.pyplot as plt

# ---- set your file path here ----
dens = 0.3
runs = range(1, 5)

plt.rcParams.update({
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.titlesize": 16
})

pat = re.compile(
    r"Step:\s*(\d+).*?E\*/N:\s*([-\d\.Ee+]+).*?P\*:\s*([-\d\.Ee+]+)"
)


def read_series(fname):
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
        raise RuntimeError(f"No data found in {fname}. Check the file format or regex pattern.")

    return np.array(steps), np.array(e_per_n), np.array(p_star)


# ---- plots ----
plt.figure(1, figsize=(12, 6))
for i in runs:
    fname = rf"C:\Users\vigne\OneDrive\Desktop\CHE496\Project_1\dens-{dens}\indep-{i}.{dens}.Monte_Carlo.log"
    print(f"Processing file: {fname}")
    steps, e_per_n, _ = read_series(fname)
    plt.plot(steps, e_per_n, linewidth=1, label=f"Run-{i}, <E*/N>={np.mean(e_per_n):.3f}")

plt.xlabel("Step")
plt.title(f"E*/N vs Step rho*= {dens}")
plt.xlim(0, 5e6)
plt.ylabel("E*/N")
plt.grid(True)
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

plt.savefig(f"E_vs_steps_{dens}.png", dpi=300)

plt.figure(2, figsize=(12, 6))
for i in runs:
    fname = rf"C:\Users\vigne\OneDrive\Desktop\CHE496\Project_1\dens-{dens}\indep-{i}.{dens}.Monte_Carlo.log"
    steps, _, p_star = read_series(fname)
    plt.plot(steps, p_star, linewidth=1, label=f"Run-{i}, <P*>={np.mean(p_star):.3f}")

plt.xlabel("Step")
plt.ylabel("P*")
plt.title(f"P* vs Step rho*= {dens}")
plt.xlim(0, 5e6)
plt.grid(True)
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()

plt.savefig(f"P_vs_steps_{dens}.png", dpi=300)


