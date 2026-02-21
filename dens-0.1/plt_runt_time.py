import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


REAL_TIME_PATTERN = re.compile(r"^\s*real\s+([0-9]*\.?[0-9]+)", re.IGNORECASE)


def read_real_time_seconds(timing_file: Path) -> float:
	with timing_file.open("r", encoding="utf-8", errors="ignore") as handle:
		for line in handle:
			match = REAL_TIME_PATTERN.search(line)
			if match:
				return float(match.group(1))
	raise RuntimeError(f"Could not find 'real' runtime in {timing_file}")


def density_from_folder(folder_name: str) -> float:
	match = re.fullmatch(r"dens-(\d*\.?\d+)", folder_name)
	if not match:
		raise ValueError(f"Invalid density folder name: {folder_name}")
	return float(match.group(1))


def get_timing_files(dens_dir: Path) -> list[Path]:
	indep_files = sorted(dens_dir.glob("indep-*.Monte_Carlo_timing.dat"))
	if indep_files:
		return indep_files
	return sorted(dens_dir.glob("Monte_Carlo_timing.dat"))


def main() -> None:
	script_dir = Path(__file__).resolve().parent
	project_root = script_dir.parent

	density_dirs = sorted(
		(path for path in project_root.glob("dens-*") if path.is_dir()),
		key=lambda path: density_from_folder(path.name),
	)

	rows = []
	print("Density summary (seconds):")

	for dens_dir in density_dirs:
		density = density_from_folder(dens_dir.name)
		timing_files = get_timing_files(dens_dir)

		if not timing_files:
			print(f"  rho*={density:.1f} -> no timing files found, skipping")
			continue

		runtimes = np.array([read_real_time_seconds(path) for path in timing_files], dtype=float)
		mean_val = float(np.mean(runtimes))
		std_val = float(np.std(runtimes, ddof=1)) if runtimes.size > 1 else 0.0
		rows.append((density, mean_val, std_val))

		print(f"  rho*={density:.1f} | n={runtimes.size} | mean={mean_val:.3f} s | std={std_val:.3f} s")

	if not rows:
		raise RuntimeError("No timing data found in any dens-* folder.")

	data = np.array(rows, dtype=float)
	densities = data[:, 0]
	mean_runtime = data[:, 1]
	std_runtime = data[:, 2]

	fig, ax = plt.subplots(figsize=(8, 5))
	ax.errorbar(
		densities,
		mean_runtime,
		yerr=std_runtime,
		marker="o",
		linestyle="-",
		linewidth=1.5,
		capsize=4,
	)
	ax.set_xlabel(r"Reduced Density, $\rho^*$")
	ax.set_ylabel("Runtime (s)")
	ax.set_title("Runtime vs Reduced Density")
	ax.grid(True)
	ax.legend()
	fig.tight_layout()

	output_plot = script_dir / "runtime_vs_density.png"
	fig.savefig(output_plot, dpi=300)
	plt.close(fig)
	print(f"Saved plot: {output_plot}")

	output_table = script_dir / "runtime_vs_density_summary.csv"
	np.savetxt(
		output_table,
		data,
		delimiter=",",
		header="density,mean_runtime_s,std_runtime_s",
		comments="",
	)
	print(f"Saved table: {output_table}")


if __name__ == "__main__":
	main()
