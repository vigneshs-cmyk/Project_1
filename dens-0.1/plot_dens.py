import re
from pathlib import Path
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt

df = pd.read_excel("heat_capacity_exp.xlsx")
density_exp_cv = df["rho*"].to_numpy()
cv_exp = df["Cv/N"].to_numpy()
cv_exp_err = df["Cv/N_err"].to_numpy()

dg = pd.read_excel("Energy_exp.xlsx")
density_exp = dg["rho*"].to_numpy()
E_exp = dg["E*/N"].to_numpy()
E_exp_err = dg["E*/N_err"].to_numpy()
P_exp = dg["P*"].to_numpy()
P_exp_err = dg["P*_err"].to_numpy()
mu_exp = dg["mu_xs"].to_numpy()
mu_exp_err = dg["mu_xs_err"].to_numpy()


DENSITIES = np.round(np.linspace(0.1, 0.9, 9), 1)
RUNS = range(1, 5)

plt.rcParams.update(
    {
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.titlesize": 16,
    }
)

number = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
step_pat = re.compile(r"Step:\s*(\d+)", re.IGNORECASE)
e_pat = re.compile(rf"E\*/N:\s*({number})", re.IGNORECASE)
p_pat = re.compile(rf"P\*:\s*({number})", re.IGNORECASE)
cv_pat = re.compile(rf"(?:Cv\*/N_xs|Cv_xs):\s*({number})", re.IGNORECASE)
mu_pat = re.compile(rf"(?:Mu\*_xs|mu_xs):\s*({number})", re.IGNORECASE)


def read_series(log_path: Path):
    steps, e_vals, p_vals, cv_vals, mu_vals = [], [], [], [], []

    with log_path.open("r", errors="ignore") as handle:
        for line in handle:
            if "Step:" not in line:
                continue

            m_step = step_pat.search(line)
            m_e = e_pat.search(line)
            m_p = p_pat.search(line)
            m_cv = cv_pat.search(line)
            m_mu = mu_pat.search(line)

            if not (m_step and m_e and m_p and m_cv and m_mu):
                continue

            steps.append(int(m_step.group(1)))
            e_vals.append(float(m_e.group(1)))
            p_vals.append(float(m_p.group(1)))
            cv_vals.append(float(m_cv.group(1)))
            mu_vals.append(float(m_mu.group(1)))

    if not steps:
        raise RuntimeError(f"No valid data found in {log_path}")

    return (
        np.asarray(steps),
        np.asarray(e_vals),
        np.asarray(p_vals),
        np.asarray(cv_vals),
        np.asarray(mu_vals),
    )


def mean_and_sem(values):
    values = np.asarray(values, dtype=float)
    mean_val = np.mean(values)
    if values.size <= 1:
        return mean_val, np.nan
    sem_val = np.std(values, ddof=1) / np.sqrt(values.size)
    return mean_val, sem_val


def main():
    project_root = Path(__file__).resolve().parent.parent

    e_mean, e_sem = [], []
    p_mean, p_sem = [], []
    cv_mean, cv_sem = [], []
    mu_mean, mu_sem = [], []

    for dens in DENSITIES:
        dens_dir = project_root / f"dens-{dens:.1f}"

        e_runs, p_runs, cv_runs, mu_runs = [], [], [], []
        for run in RUNS:
            log_file = dens_dir / f"indep-{run}.{dens:.1f}.Monte_Carlo.log"
            print(f"Processing: {log_file}")

            _, e_vals, p_vals, cv_vals, mu_vals = read_series(log_file)

            e_runs.append(np.mean(e_vals))
            p_runs.append(np.mean(p_vals))
            cv_runs.append(np.mean(cv_vals))
            mu_runs.append(np.mean(mu_vals))

        em, es = mean_and_sem(e_runs)
        pm, ps = mean_and_sem(p_runs)
        cvm, cvs = mean_and_sem(cv_runs)
        mum, mus = mean_and_sem(mu_runs)

        e_mean.append(em)
        e_sem.append(es)
        p_mean.append(pm)
        p_sem.append(ps)
        cv_mean.append(cvm)
        cv_sem.append(cvs)
        mu_mean.append(mum)
        mu_sem.append(mus)

    output_dir = Path(__file__).resolve().parent

    figure_specs = [
        (e_mean, e_sem, r"$E^{*}/N$", "Average Reduced Energy", "E_vs_density_sem.png"),
        (p_mean, p_sem, r"$P^{*}$", "Average Reduced Pressure", "P_vs_density_sem.png"),
        (cv_mean, cv_sem, r"$C_{v,XS}^{*}$", "Average Excess Heat Capacity", "Cv_xs_vs_density_sem.png"),
        (mu_mean, mu_sem, r"$\mu_{XS}^{*}$", "Average Excess Chemical Potential", "mu_xs_vs_density_sem.png"),
    ]

    for mean_vals, sem_vals, y_label, title, file_name in figure_specs:
        fig = plt.figure(figsize=(8, 5))
        plt.errorbar(DENSITIES, mean_vals, yerr=sem_vals, marker=".", capsize=4, linewidth=1.4, label="Simulation")  # Error bars for simulation data
        if y_label == r"$C_{v,XS}^{*}$":
            plt.errorbar(density_exp_cv, cv_exp, yerr=cv_exp_err, marker='.',capsize=4, linewidth=1.4, label=r"Experimental")  # Example error bar for experimental data
            plt.legend()
        elif y_label == r"$E^{*}/N$":
            plt.errorbar(density_exp, E_exp, yerr=E_exp_err, marker='.',capsize=4, linewidth=1.4, label=r"Experimental")  # Example error bar for experimental data
            plt.legend()
        elif y_label == r"$P^{*}$":
            plt.errorbar(density_exp, P_exp, yerr=P_exp_err, marker='.',capsize=4, linewidth=1.4, label=r"Experimental")  # Example error bar for experimental data
            plt.legend()
        elif y_label == r"$\mu_{XS}^{*}$":
            plt.errorbar(density_exp, mu_exp, yerr=mu_exp_err, marker='.',capsize=4, linewidth=1.4, label=r"Experimental")  # Example error bar for experimental data
            plt.legend()
        plt.xlabel(r"Density, $\rho^*$")
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid(True)
        plt.tight_layout()

        out_path = output_dir / file_name
        fig.savefig(out_path, dpi=300)
        print(f"Saved: {out_path}")
        plt.close(fig)


if __name__ == "__main__":
    main()
