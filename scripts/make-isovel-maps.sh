ranges='+030+040 +020+030 +010+020 +000+010 -010+000 -020-010 -030-020 -040-030 -050-040 -060-050 -070-060 -080-070 -090-080'
# blueranges='-030-010 -050-030 -070-050'
# farblueranges='-090-070 -110-090 -130-110'
# for vrange in $redranges $blueranges $farblueranges; do
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py ha $vrange
    python scripts/turtle-spectral-map.py nii $vrange
done
