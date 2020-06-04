#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("log_file", type=str, nargs='+', help="Input log file(s)")
args = p.parse_args()

logs = [pd.read_csv(f, delim_whitespace=True) for f in args.log_file]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_file)]

fig, axes = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    # Compute smoothed velocity
    log['velocity'] = np.gradient(log['y'], log['time'])

    # Remove outliers
    v_low = log["velocity"].abs().quantile(0.1)
    v_hi = log["velocity"].abs().quantile(0.9)
    mask = np.logical_or(log['velocity'].abs() < 0.5 * v_low,
                         log['velocity'].abs() > 2 * v_hi)
    log['velocity'].mask(mask, np.nan, inplace=True)

    log.plot('time', 'y', ax=axes[0, 0], label=f'y{i}')
    log.plot('time', 'velocity', ax=axes[0, 1], label=f'velocity{i}')
    log.plot('time', 'max(E)', ax=axes[1, 0], label=f'max(E){i}')
    log.plot('time', 'x.2', ax=axes[1, 1], label=f'radius{i}')
plt.legend()
plt.show()