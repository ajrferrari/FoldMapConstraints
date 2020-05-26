Help()
{
echo "Extract constraints file to be used in refinement of a protein model"
echo
echo "Syntax: extract_constraints.sh database_folder query_protein_model mapconstraints_folder lovoalign_executable function_to_define_constraints"
echo
echo "Functions available FLAT_HARMONIC LINEAR_PENALTY USOG"
}

while getopts ":h" option; do
	case $option in
		h) # display Help
		Help
		exit;;
	esac
done



database=$(readlink -f $1)
query=$(readlink -f $2)
mapconstraints=$(readlink -f $3)
reference=$database/reference.pdb
clean_pdbs=$(readlink -f $mapconstraints/src/clean_pdb.py)
extract_constraints=$(readlink -f $mapconstraints/src/extract_constraints.py)
lovoalign=$( readlink -f $4 ) 


# Clean query pdb

python2.7 $clean_pdbs $query A

mv *_A.pdb query.pdb

# Align query to reference pdb

$lovoalign -p1 query.pdb -p2 $reference -o query_aligned.pdb

# Run extract constraints

python3.7 $extract_constraints -database $database/ -query query_aligned.pdb -func $5

# Remove unnecessary files
rm query*pdb *fasta
