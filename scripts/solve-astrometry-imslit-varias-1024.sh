unset files
declare -a files=( 'spm003_b' 'spm004_b' 'spm017_b' 'spm019_b' 'spm021_b' 'spm111_cr' )
D=~/Dropbox/Papers/LL-Objects/NGC6210/Varias-temporadas
for f in ${files[*]}; do
    solve-field --ra 251.1230 --dec 23.7999 --radius 0.1 --scale-units arcsecperpix --scale-low 0.3 --scale-high 0.4 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/${f}.fits --odds-to-solve 1e6 --quad-size-min 0.03 --code-tol 0.01
done
