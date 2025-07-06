"""
make_density_animation.py
-------------------------

Create an MP4 animation that visualises an SPH-style density field.
Each particle is rendered as a 2-D Gaussian “blob” on a fixed pixel grid
and the blobs are summed to form a smooth surface.

Assumes the Go simulation has written files named
    positions_0000.csv, positions_0001.csv, …
in the working directory, where every file contains lines of

    index,x,y        (no header row)

Dependencies
============
pip install numpy pandas matplotlib
and have FFmpeg on your PATH (https://ffmpeg.org/).

Run
===
python make_density_animation.py        # produces density.mp4
"""
from pathlib import Path
import glob
import numpy as np
import pandas as pd

from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# -----------------------------------------------------------------------------
# Load snapshots
# -----------------------------------------------------------------------------
files = sorted(glob.glob("out/positions_*.csv"))
if not files:
    raise FileNotFoundError("No positions_*.csv files found.")

frames = [
    pd.read_csv(
        f,
        usecols=[1, 2],                 # only x and y
        header=None,
        names=["x", "y"],
        dtype=np.float64,               # force float64, errors → NaN
        on_bad_lines="skip",            # ignore malformed rows
    ).to_numpy()
    for f in files
]
data = np.stack(frames)                      # shape = (F, P, 2)
F, P, _ = data.shape
print(f"Loaded {F} frames, {P} particles per frame.")

# -----------------------------------------------------------------------------
# Grid and kernel parameters
# -----------------------------------------------------------------------------
xmin, ymin = 0, 0 #data.min(axis=(0, 1))
xmax, ymax = 10, 10 #data.max(axis=(0, 1))

grid_size = 200                      # number of pixels along x and y
sigma     = (xmax - xmin) / 25.0     # std-dev of the Gaussian kernel

x_lin = np.linspace(xmin, xmax, grid_size)
y_lin = np.linspace(ymin, ymax, grid_size)
X, Y  = np.meshgrid(x_lin, y_lin)    # shape = (grid_size, grid_size)

# -----------------------------------------------------------------------------
# Animation setup
# -----------------------------------------------------------------------------
fig, ax = plt.subplots()
im = ax.imshow(
    np.zeros_like(X),
    extent=[xmin, xmax, ymin, ymax],
    origin="lower",
    cmap="Blues",
    animated=True,
)

# ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_ylim(0, 10)
ax.set_xlim(0, 10)

def init():
    im.set_data(np.zeros_like(X))
    return im,

def update(frame):
    """
    Fast Gaussian accumulation using separability + BLAS.
    Memory: 2 × (grid_size · P) instead of grid_size² · P.
    """
    pos = data[frame]                       # shape (P, 2)

    # Pre-compute 1-D kernels for *all* particles at once
    # gx, gy both shape = (grid_size, P)
    gx = np.exp(-((x_lin[:, None] - pos[:, 0]) ** 2) / (2.0 * sigma * sigma))
    gy = np.exp(-((y_lin[:, None] - pos[:, 1]) ** 2) / (2.0 * sigma * sigma))

    # Sum over particles via matrix multiplication: (grid × P) @ (P × grid)
    density = gy @ gx.T                     # shape = (grid_size, grid_size)

    im.set_data(density)
    im.set_clim(vmin=0.0, vmax=density.max())   # consistent colour scale

    # iso = 0.6 * density.max()         # 60 % of peak works well for 200–500 pts
    # mask = density >= iso
    # im.set_data(np.where(mask, density, np.nan))   # draw NaN as fully transparent
    # im.set_clim(vmin=iso, vmax=density.max())      # ice-blue shell

    # if pos.shape[0] > 1:
    #     # crude: mean distance to k-nearest neighbours, k = 4
    #     k = 4
    #     dists = np.partition(
    #         np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2), k, axis=1
    #     )[:, 1:k+1]                     # skip self-distance at index 0
    #     sigma_frame = 0.4 * dists.mean()  # 0.4× looks “surface-like”
    # else:
    #     sigma_frame = sigma             # fallback for single particle

    # denom = 2.0 * sigma_frame * sigma_frame
    # gx = np.exp(-((x_lin[:, None] - pos[:, 0]) ** 2) / denom)
    # gy = np.exp(-((y_lin[:, None] - pos[:, 1]) ** 2) / denom)
    # density = gy @ gx.T

    # # im.set_data(density)
    # # im.set_clim(vmin=0.0, vmax=density.max())   # consistent colour scale

    # iso = 0.6 * density.max()         # 60 % of peak works well for 200–500 pts
    # mask = density >= iso
    # im.set_data(np.where(mask, density, np.nan))   # draw NaN as fully transparent
    # im.set_clim(vmin=iso, vmax=density.max())      # ice-blue shell

    return im,


anim = FuncAnimation(
    fig,
    update,
    frames=F,
    init_func=init,
    # blit=True,
    # interval=30,         # ~25 fps
)

pbar = tqdm(total=F, desc="Encoding frames")

def _progress(frame_number, total):
    pbar.update(frame_number - pbar.n)

# -----------------------------------------------------------------------------
# Save animation
# -----------------------------------------------------------------------------
output = Path("density.mp4")
print(f"Writing {output} ...")
writer = FFMpegWriter(fps=30)#, codec="libx264", bitrate=1800)
anim.save(output, writer=writer, dpi=150, progress_callback=_progress)
pbar.close()
print("Done.")
