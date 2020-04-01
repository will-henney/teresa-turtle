ranges='-020+000 -040-020 -060-040 -080-060 -080+000'
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py heii $vrange maps20
done
