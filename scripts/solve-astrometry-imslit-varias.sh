unset files
declare -a files=( 'obj1002_bcr ' 'obj1006_bcr ' 'obj1009_bcr ' 'obj1012_bcr ' 'obj1017_bcr ' 'spm110_u' 'spm115_bcr ' 'spm118_bcr ' 'spm121_bcr ' 'spm124_bcr ' 'spm127_bcr ' 'spm134_bcr ' 'spm165_bcr ' 'spm223_bcr ' 'spm229_bcr ' 'spm318_bcr ' 'spm328_bcr ' 'spm331_bcr ' 'spm414_bcr ' 'spm236_bcr ' 'spm238_bcr ' 'spm242_bcr ' 'spm251_bcr ' 'spm256_bcr ' 'spm293_bcr ' 'spm296_bcr ' 'spm302_bcr ' 'spm303_bcr ' 'spm305_bcr ' 'spm306_bcr ' 'spm350_bcr ' 'spm358_bcr ' 'spm407_bcr' 'spm111_bcr' 'spm116_bcr' 'spm003_b' 'spm004_b' 'spm017_b' )
D=~/Dropbox/Papers/LL-Objects/NGC6210/Varias-temporadas
for f in ${files[*]}; do
    solve-field --ra 251.1230 --dec 23.7999 --radius 0.1 --scale-units arcsecperpix --scale-low 0.5 --scale-high 0.8 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/${f}.fits --odds-to-solve 1e6 --quad-size-min 0.01 --code-tol 0.03
done
