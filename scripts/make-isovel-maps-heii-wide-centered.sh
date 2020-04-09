ranges='-060-020'
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py heii $vrange maps20
done
