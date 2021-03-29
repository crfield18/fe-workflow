for i in *parm7; do
	base=$(basename ${i}); ext="${base##*.}"; base="${base%.*}"
	parmutils-Organize.py -p ${base}.parm7 -c  ${base}.rst7 -m ":L1" -o "${base}_com" -e1 '!:WAT,K+,Na+,Cl-' -o1 "${base}_com_dry"

	cat << EOF > cpp.in
parm ${base}.parm7
trajin ${base}.rst7
trajout ${base}.pdb
go
quit
EOF
#	cpptraj < cpp.in
done

