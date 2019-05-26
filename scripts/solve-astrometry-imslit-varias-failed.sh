unset files
declare -a files=( 'obj1002_bcr' 'obj1006_bcr' 'spm293_bcr' )
D=~/Dropbox/Papers/LL-Objects/NGC6210/Varias-temporadas
for f in ${files[*]}; do
    solve-field --ra 251.1230 --dec 23.7999 --radius 0.1 --scale-units arcsecperpix --scale-low 0.5 --scale-high 0.8 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/${f}.fits --odds-to-solve 1e3 --quad-size-min 0.01 --code-tol 0.03
done
