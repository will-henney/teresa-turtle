files='0128'
D=~/Dropbox/Papers/LL-Objects/NGC6210/Temporada2015
for f in ${files[*]}; do
    solve-field --ra 251.1230 --dec 23.7999 --radius 0.1 --scale-units arcsecperpix --scale-low 0.35 --scale-high 0.4 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/spm${f}o_b.fits --odds-to-solve 15 --quad-size-min 0.03 --code-tol 0.03 -v
done
