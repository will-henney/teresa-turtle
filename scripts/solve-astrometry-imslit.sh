unset files
declare -a files=( 'spm0030' 'spm0035' 'spm0036' 'spm0041' 'spm0042' 'spm0047' 'spm0048' 'spm0053' 'spm0056' 'spm0061' 'spm0109' 'spm0114' 'spm0115' 'spm0120' 'spm0121' 'spm0123' 'spm0128' 'spm0172' 'spm0177' 'spm0178' 'spm0183' 'spm0185' 'spm0190' )
D=~/Dropbox/Papers/LL-Objects/NGC6210/Temporada2015
for f in ${files[*]}; do
    solve-field --ra 251.1230 --dec 23.7999 --radius 1.0 --scale-units arcsecperpix --scale-low 0.35 --scale-high 0.4 --dir . --new-fits '%s-wcs.fits' --no-tweak --overwrite $D/${f}o_b.fits
done
