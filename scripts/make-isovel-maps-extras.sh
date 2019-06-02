ranges='+040+050 -100-090 -110-100'
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py ha $vrange
    python scripts/turtle-spectral-map.py nii $vrange
done
