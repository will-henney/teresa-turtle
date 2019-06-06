ranges='+040+060 +020+040 +000+020 -020+000 -040-020 -060-040 -080-060 -100-080 -120-100'
# blueranges='-030-010 -050-030 -070-050'
# farblueranges='-090-070 -110-090 -130-110'
# for vrange in $redranges $blueranges $farblueranges; do
for vrange in $ranges; do
    python scripts/turtle-spectral-map.py oiii $vrange maps-oiii-20
done
