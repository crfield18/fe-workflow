for i in 1h1q.parm7 1h1r.parm7 1oiu.parm7; do
	base=$(basename ${i}); ext="${base##*.}"; base="${base%.*}"
	parmutils-Organize.py -p ${base}.parm7 -c  ${base}.rst7 -m ":L1" -o "${base}_com" -e1 '!:WAT,K+,Na+,Cl-' -o1 "${base}_com_dry" -e2 ':L1' -o2 "${base}_lig"

	cat << EOF > cpp.in
parm ${base}.parm7
trajin ${base}.rst7
trajout ${base}.pdb
go
quit
EOF
#	cpptraj < cpp.in
done

