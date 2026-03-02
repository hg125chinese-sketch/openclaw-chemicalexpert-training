# DRD2 predictors (scaffold split)
N_total=5067; train=4053 valid=507 test=507

| Model | RMSE | MAE | R^2 | RMSE/std(y) | Params |
|---|---:|---:|---:|---:|---:|
| RF (Morgan, 500 trees) | 0.586 | 0.449 | 0.567 | 0.658 | n/a |
| GCN (3-layer, hidden=128) | 0.821 | 0.673 | 0.150 | 0.922 | 59009 |
