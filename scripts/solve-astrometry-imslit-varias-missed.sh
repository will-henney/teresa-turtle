unset files
declare -a files=( 'spm116-2_bcr' )
D=~/Dropbox/Papers/LL-Objects/NGC6210/Varias-temporadas
for f in ${files[*]}; do
    solve-field -v --ra 251.1230 --dec 23.7999 --radius 0.1 --scale-units arcsecperpix --scale-low 0.3 --scale-high 0.7 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/${f}.fits --odds-to-solve 1e2 --code-tol 0.03
done
