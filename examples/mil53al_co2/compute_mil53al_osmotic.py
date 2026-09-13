#!/usr/bin/env python3
"""Published-parameter osmotic model for CO2 adsorption in MIL-53(Al).

This is a reproducible two-state benchmark, not an atomistic RASPA run.  The
parameters are transcribed from Boutin et al. (arXiv:1904.11921, Table 1).
The source writes K = K' exp(-DeltaH_ads / RT) with DeltaH_ads as a negative
adsorption enthalpy.  Here ``h_ads`` stores the positive magnitude, so the
implemented expression is K = K' exp(h_ads / RT).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


R_KJ = 8.314462618e-3  # kJ mol^-1 K^-1
OUT = Path(__file__).resolve().parent

# Host thermodynamics from the published MIL-53(Al) osmotic model.
HOST_DH = 15.0  # kJ mol^-1, lp - np convention
HOST_DS = 0.074  # kJ mol^-1 K^-1, lp - np convention

# CO2 parameters from Table 1: K = K0 exp(h_ads / RT), Nmax = a - b T.
CO2 = {
    "np": {"h_ads": 38.8, "K0": 1.1e-5, "a": 2.3, "b": 0.0002},
    "lp": {"h_ads": 26.0, "K0": 8.5e-5, "a": 16.8, "b": 0.025},
}


def phase_parameters(phase: str, temperature: float) -> tuple[float, float]:
    """Return Langmuir K (bar^-1) and Nmax (molecules per unit cell)."""
    p = CO2[phase]
    rt = R_KJ * temperature
    k = p["K0"] * math.exp(p["h_ads"] / rt)
    nmax = p["a"] - p["b"] * temperature
    return k, nmax


def uptake(phase: str, pressure_bar: float, temperature: float) -> float:
    """Rigid-host Langmuir uptake in molecules per unit cell."""
    k, nmax = phase_parameters(phase, temperature)
    return nmax * k * pressure_bar / (1.0 + k * pressure_bar)


def omega_ads(phase: str, pressure_bar: float, temperature: float) -> float:
    """Adsorption grand potential relative to the empty phase, kJ per u.c."""
    k, nmax = phase_parameters(phase, temperature)
    return -R_KJ * temperature * nmax * math.log1p(k * pressure_bar)


def delta_omega(pressure_bar: float, temperature: float) -> float:
    """Omega(lp) - Omega(np), in kJ per unit cell."""
    delta_g_host = HOST_DH - temperature * HOST_DS
    return delta_g_host + omega_ads("lp", pressure_bar, temperature) - omega_ads(
        "np", pressure_bar, temperature
    )


def roots_for_temperature(temperature: float) -> list[float]:
    """Find all sign changes on a logarithmic pressure grid."""
    grid = np.logspace(-10, 4, 3000)
    values = np.array([delta_omega(float(p), temperature) for p in grid])
    roots: list[float] = []
    for i in range(len(grid) - 1):
        if values[i] == 0.0:
            roots.append(float(grid[i]))
            continue
        if values[i] * values[i + 1] > 0.0:
            continue
        lo, hi = float(grid[i]), float(grid[i + 1])
        flo = delta_omega(lo, temperature)
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            fmid = delta_omega(mid, temperature)
            if flo * fmid <= 0.0:
                hi = mid
            else:
                lo, flo = mid, fmid
        roots.append(math.sqrt(lo * hi))
    return roots


def write_phase_boundary() -> None:
    rows = []
    for temperature in range(195, 344):
        roots = roots_for_temperature(float(temperature))
        rows.append(
            {
                "temperature_K": temperature,
                "lp_to_np_eq_bar": roots[0] if len(roots) == 2 else "",
                "np_to_lp_eq_bar": roots[-1] if len(roots) >= 1 else "",
                "number_of_crossings": len(roots),
            }
        )
    with (OUT / "phase_boundary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def make_plot() -> None:
    temperatures = np.arange(195.0, 344.0, 1.0)
    low, high = [], []
    for temperature in temperatures:
        roots = roots_for_temperature(float(temperature))
        low.append(roots[0] if len(roots) == 2 else np.nan)
        high.append(roots[-1] if roots else np.nan)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    ax0.plot(temperatures, low, color="#2f6f9f", label="lp -> np")
    ax0.plot(temperatures, high, color="#b24c35", label="np -> lp (gate opening)")
    ax0.set_yscale("log")
    ax0.set_xlabel("Temperature (K)")
    ax0.set_ylabel("Equilibrium pressure (bar)")
    ax0.set_title("MIL-53(Al) + CO2: osmotic phase boundary")
    ax0.grid(True, which="both", alpha=0.25)
    ax0.legend(frameon=False)

    pressure = np.logspace(-7, 1, 500)
    for temperature, color in [(223.0, "#b24c35"), (254.0, "#2f6f9f"), (298.0, "#555555")]:
        ax1.plot(
            pressure,
            [delta_omega(float(p), temperature) for p in pressure],
            label=f"{temperature:.0f} K",
            color=color,
        )
    ax1.axhline(0.0, color="black", lw=0.8)
    ax1.set_xscale("log")
    ax1.set_xlabel("CO2 pressure (bar)")
    ax1.set_ylabel("Delta Omega = Omega_lp - Omega_np (kJ/u.c.)")
    ax1.set_title("Free-energy crossing criterion")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(frameon=False)

    fig.savefig(OUT / "mil53al_co2_osmotic_phase_boundary.png", dpi=220)
    plt.close(fig)


def print_summary() -> None:
    print("MIL-53(Al) + CO2 published-parameter osmotic benchmark")
    print("host convention: DeltaG_host = G_lp - G_np = 15 - 0.074*T kJ/u.c.")
    for temperature in (204.0, 223.0, 254.0, 273.0, 298.0):
        roots = roots_for_temperature(temperature)
        k_np, n_np = phase_parameters("np", temperature)
        k_lp, n_lp = phase_parameters("lp", temperature)
        print(
            f"T={temperature:.0f} K  K_np={k_np:.5g}  K_lp={k_lp:.5g}  "
            f"Nmax_np={n_np:.4g}  Nmax_lp={n_lp:.4g}  "
            f"equilibrium crossings (bar)={', '.join(f'{r:.6g}' for r in roots) or 'none'}"
        )


if __name__ == "__main__":
    write_phase_boundary()
    make_plot()
    print_summary()
