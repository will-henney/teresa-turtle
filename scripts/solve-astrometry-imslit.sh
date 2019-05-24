unset files
declare -a files=( 'spm0030' 'spm0035' 'spm0036' 'spm0041' )
echo ${files[*]}
for f in ${files[*]}; do
    echo $f
done
