#!/bin/bash
# Sequential production sweeps (defaults: 256x256, 500 cycles, 200 steps/cycle).
cd /home/kashaf/duffing-parameter-space
for g in 0.001 0.150 0.250 0.300; do
    echo "=== starting gamma=$g at $(date -Is) ==="
    .venv/bin/python duffing_parameter_space.py --gamma "$g" --backend numba \
        --out "results/duffing_gamma${g}_n256.png" > "logs/gamma${g}.log" 2>&1
    echo "=== gamma=$g exit=$? at $(date -Is) ==="
done
echo "ALL DONE $(date -Is)"
