ranges='+040+060 +020+040 +000+020 -020+000 -040-020 -060-040 -080-060 -100-080 -120-100'
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py ha $vrange maps20
    python scripts/turtle-spectral-map.py nii $vrange maps20
done
