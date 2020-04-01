ranges='-030-010 -050-030 -070-050'
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py heii $vrange maps20
done
