# Reference trajectories

Long channel runs (walls in y, periodic x and z; H = 32, R = 4, Re = 0.48,
start on the centreline at 30 degrees, soft wall potential eps = 0.2, hc = 2,
planar motion). They let the tutorial notebook show long-time behaviour without
waiting for the runs. Regenerate them with:

| file | command | steps | wall time (1 core) |
|---|---|---|---|
| `channel_neutral_traj.dat` | `./bin/lbm_squirmer inputs/channel_neutral.in steps=40000` | 40 000 | ~8 min |
| `channel_puller_traj.dat`  | `./bin/lbm_squirmer inputs/channel_puller.in`  | 90 000 | ~17 min |
| `channel_pusher_traj.dat`  | `./bin/lbm_squirmer inputs/channel_pusher.in`  | 30 000 | ~6 min |

Columns: `t X Y Z ex ey ez Ux Uy Uz Wx Wy Wz Fwall_y Fwall_z` (lattice units,
X is unwrapped along the periodic x-direction).
