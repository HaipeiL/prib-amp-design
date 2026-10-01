# Model weights

Frozen model weights are published with this repository for reproducibility (Git LFS is recommended). Expected files and SHA256 values:

| file | SHA256 |
|---|---|
| `m2_step_0010000.pt` | `fdc17a1052dcd00678bbe317dbd7e0d29ee36ba36e5cc6221707f07ffb8e2141` |
| `pmic_seed_20260903_step_008000.pt` | `a690836e8f9fe68e318a012fa35bb1c80615d608ea1173466ac6e49479eb181a` |
| `pmic_seed_20260917_step_005000.pt` | `5eaf2024023c53fef9d5ba58ce3d375fc448a3eaad62cf3365ba6dd8f2ac1ea3` |
| `pmic_seed_20261001_step_007000.pt` | `2881231b3a631ba874f0a9d6a287b864f9c99a2d882130a8238878cdf58e223b` |

Run `scripts/export_weights_from_cluster.sh` on the PRIB Cluster42 workspace to copy and verify the exact frozen files into this directory.
