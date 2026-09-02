# A Mesh-Based Knee Contact Pressure Model for Tracking and Predictive Simulations

This repository contains the code for tracking and predictive simulations formulated as optimal control problems using a mesh-based knee contact pressure model. The code accompanies the following paper:

-

---

## Tracking Simulations

The main script for running the tracking simulations is `OCP_GC/RunAllSimulations.m`.

To reproduce the results presented in the article, the user can configure the following parameters:

- **`movements`**: `ngait_og1`, `ngait_og5`, `bouncy1`, `bouncy4`, `mtpgait3`, `mtpgait9`, `ngait_tm_fast1`, or `ngait_tm_set1`. The selected movements should be specified in a cell array.

- **`damping_Pairs`**: `[0 0]`, `[0.001 1000]`, or `[0.0005 500]`, representing pairs of damping parameters for passive forces and moments at the knee.

- **`kmax_list`**: `1e3`, `1e4`, or `"Max"`, representing the $k_{\mathrm{max}}$ parameter. `kmax_list` is a cell array. If more than one value is specified, all combinations will be run.

- **`kpress_list`**: `1e3`, `1e4`, or `1e5`, representing the $k_{\mathrm{press}}$ parameter. `kpress_list` is a vector. If more than one value is specified, all combinations will be run.

- **`kcheck_list`**: `1e2`, `1e3`, or `1e4`.

- **`tifaces_values` and `femfaces_values`**: accept the pairs $49 \times 171$, $75 \times 258$, or $100 \times 342$, representing the different mesh resolutions.

The optimal control problems (OCPs) for the tracking simulations are formulated in `OCP_GC/TrackSim_3D_GC_v2.m`.

## Predictive Simulations

The main script for running the predictive simulations is `OCP_GC/PredSim_3D_GC_v2.m`.

To reproduce the results presented in the article while varying only the weight penalizing the maximum knee contact pressure, the user only needs to set `W.minpressures` to `0`, `0.03`, `0.3`, or `3`.

The results of the predictive simulations can be visualized in [this video](Results/SupplementaryVideo.mp4).

