#!/bin/bash
cd /home/kashaf/duffing-parameter-space
set -- "0.250 1.0 1.6 0.18 0.30" "0.300 1.0 1.6 0.21 0.30"
for spec in "$@"; do
    read -r g o0 o1 g0 g1 <<< "$spec"
    echo "=== zoom gamma=$g Omega[$o0,$o1] G[$g0,$g1] at $(date -Is) ==="
    .venv/bin/python duffing_parameter_space.py --gamma "$g" --n 256 --cycles 500 --spc 200 \
        --backend numba --omega "$o0" "$o1" --G "$g0" "$g1" \
        --out "results/zoom/duffing_gamma${g}_zoom_n256.png" > "logs/zoom_gamma${g}.log" 2>&1
    echo "=== gamma=$g exit=$? at $(date -Is) ==="
done
echo "ZOOM DONE $(date -Is)"
