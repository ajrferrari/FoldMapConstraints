ref=$(readlink -f $1)
pdblist=$(readlink -f $2)
clean_pdbs=$(readlink -f /home/allan/softwares/rosetta_src_2019.35.60890_bundle/tools/protein_tools/scripts/clean_pdb.py)
lovoalign=$(readlink -f /mnt/d/softwares/lovoalign/bin/lovoalign)
extract_database=$(readlink -f /mnt/d/work/2020/POSTDOC/scripts/extract_database.py)



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


