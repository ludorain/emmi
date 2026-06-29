#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt


def autofocus(filename, half_window=3, top_fraction=None):
    # Input file: x  y  ey
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    ey = data[:, 2]

    # Rough maximum
    imax = np.argmax(y)

    # Local fit window
    i0 = max(0, imax - half_window)
    i1 = min(len(x) - 1, imax + half_window)

    xfit = x[i0:i1 + 1]
    yfit = y[i0:i1 + 1]

    # Optional: fit only points near the top
    if top_fraction is not None:
        mask = yfit > top_fraction * np.max(yfit)
        if np.count_nonzero(mask) >= 3:
            xfit = xfit[mask]
            yfit = yfit[mask]

    # Quadratic fit: y = a x^2 + b x + c
    a, b, c = np.polyfit(xfit, yfit, 2)

    # Vertex of parabola
    if a < 0:
        zmax = -b / (2 * a)
        ymax_fit = a * zmax**2 + b * zmax + c
    else:
        # Fallback if the fit is not concave
        zmax = x[imax]
        ymax_fit = y[imax]

    # Smooth curve for plotting
    xquad = np.linspace(xfit.min(), xfit.max(), 400)
    yquad = a * xquad**2 + b * xquad + c

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(
        x, y, yerr=ey,
        fmt='o',
        markersize=8,
        linewidth=2,
        capsize=3,
        label='Data'
    )
    ax.plot(xquad, yquad, '--', linewidth=2, label='Quadratic fit')
    ax.axvline(zmax, linestyle=':', linewidth=2, label=f'Best focus = {zmax:.4f} mm')

    ax.set_xlabel('z (mm)')
    ax.set_ylabel('laplace variance')
    ax.legend()
    fig.tight_layout()
    fig.savefig("autofocus.png", dpi=150)

    print(f"rough maximum at x = {x[imax]}")
    print(f"quadratic maximum at x = {zmax}")
    print(f"fitted maximum y = {ymax_fit}")
    print("Saved plot to autofocus.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} data.txt")
        sys.exit(1)

    autofocus(sys.argv[1], half_window=3, top_fraction=0.8)
    
