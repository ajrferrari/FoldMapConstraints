database=$(readlink -f $1)
query=$(readlink -f $2)
reference=$database/reference.pdb
clean_pdbs=$(readlink -f /home/allan/softwares/rosetta_src_2019.35.60890_bundle/tools/protein_tools/scripts/clean_pdb.py)
lovoalign=$(readlink -f /mnt/d/softwares/lovoalign/bin/lovoalign)
extract_constraints=$(readlink -f /mnt/d/work/2020/POSTDOC/scripts/map_constraints_reduced/extract_constraints.py)


# Clean query pdb

python2.7 $clean_pdbs $query A

mv *_A.pdb query.pdb

# Align query to reference pdb

$lovoalign -p1 query.pdb -p2 $reference -o query_aligned.pdb

# Run extract constraints

python3.7 $extract_constraints -database $database/ -query query_aligned.pdb
