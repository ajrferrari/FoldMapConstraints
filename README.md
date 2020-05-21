# MapConstraints

Two modules are available: extract database and extract constraints.

## Extract Database

Given a reference PDB file and a list of related PDBs, extract database will search for conserved residue pair positions (relative to the reference PDB) and avarage values (Cb-Cb Euclidean distance and Ca-Cb-Cb-Ca dihedral angle) among all PDBs.


Some details:

  The module will:
  1) clean all unecessary information from PDBs using clean_pdb.py (src/clean_pdb.py); 
  2) output aligned PDB files to the reference using LovoAlign (https://www.ime.unicamp.br/~martinez/lovoalign/home.html);
  3) for each residue pair in the reference PDB with residue relative position greater than 10:
  - try to find a Cb atom in every other PDB file so that that maximum distance of the aligned Cb atoms is less than 3 Å;
  - if a pair of Cb atoms is found, compute the Euclidean distance between them and the Ca-Cb-Cb-Ca dihedral angle;
  - create a list containing the computed Euclidean distance for each PDB file;
  - create a list containing the computed dihedral angle for each PDB file;
  4) for each residue pair in the reference PDB file check if Euclidean distance were computed for all PDBs and if the difference between the maximum and minimum value is lower than a threshold (i.e. 2 Å) and the avarage distance amoung the PDBs;
  5) do the same as 4) for the dihedral angles;
  6) write .pos and .avg files containing the pair positions (relative to the reference PDB) and avarage values for each interval of values. For example: Euclidean_1_2.pos contains pair positions for which the difference between the maximum and minimum values computed for CB-CB distance is between 1 and 2 Å.
  
  
## Extract Constraints

Given a query PDB model, the reference PDB and the database created above, extract constraints will align the query PDB model to the reference PDB, search the analogous residue pairs found to be conserved in terms of Euclidean distance or Dihedral angle and create a constraint file to be used with Rosetta. 
   
