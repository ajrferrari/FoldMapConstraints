Help()
{
echo "Extract position files from reference pdb after checking for Euclidean distances and dihedral angles patterns over list of PDBs"
echo "A database folder is created containing such files and the reference pdb for future use to generate constraints for refinement for protein models"
echo 
echo "Syntax: extract_database.sh reference_pdb pdblist mapconstraints_folder lovoalign_executable"
echo
}

while getopts ":h" option; do
	case $option in
		h) # display Help
		Help
		exit;;
	esac
done


ref=$(readlink -f $1)
pdblist=$(readlink -f $2)
mapconstraints=$(readlink -f $3)
clean_pdbs=$(readlink -f $mapconstraints/src/clean_pdb.py)
lovoalign=$( readlink -f $4 )
extract_database=$(readlink -f $mapconstraints/src/extract_database.py)



# Clean pdbs

mkdir cleaned_and_aligned

cd cleaned_and_aligned

python2.7 $clean_pdbs $ref A

mv *_A.pdb reference.pdb

for i in $(cat $pdblist); do python2.7 $clean_pdbs $i A; done

readlink -f *_A.pdb > pdbs.lst


# Align pdbs to reference pdb

for i in $(cat pdbs.lst); do $lovoalign -p1 $i -p2 reference.pdb -o $(echo $i | sed 's/_A.pdb/_aligned.pdb/g'); done

rm *_A.pdb *fasta

readlink -f  *aligned.pdb > pdbs.lst


# Create database

python3.7 $extract_database -ref reference.pdb -pdbs_list pdbs.lst

# Move database directory and remove files

mv database ..

rm -r cleaned_and_aligned

echo "Database creation completed"


