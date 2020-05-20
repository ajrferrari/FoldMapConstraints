from functions import *

parser = argparse.ArgumentParser(description='Extract database from aligned pdbs')
parser.add_argument('-ref', '--ref', type=str, required=True, help='Reference pdb')
parser.add_argument('-pdbs_list', '--pdbs_list', type=str, required=True, help='File containing list of pdbs')
args = parser.parse_args()


def main():

    start = time()
    
    print("Reading pdb files ...")
    files_list = read_pdbs_files(args.ref, args.pdbs_list)
    
    print("Calculating Euclidean distances for %s ..." %files_list[0])
    Euclidean_distances_ref = list_of_Euclidean_distances(files_list[0])
    
    print("Calculating Dihedral angles for %s ..." %files_list[0])
    Dihedral_angles_ref = list_of_dihedrals(Euclidean_distances_ref, files_list[0])
    
    print("Calculating Euclidean distances and Dihedrals angles for all %d pdbs ... " % ( len(files_list) ) )
    Distance_dict, Dihedral_dict = get_dictionaries_of_distances_and_dihedrals(Dihedral_angles_ref, files_list, files_list[0])
    
    print("Writting database files ... ")
    write_database(Distance_dict, Dihedral_dict)

    os.system('cp %s %s' % (files_list[0], 'database') )
	
    print("Finished in %f minutes" % ( (time()-start)/60 ))

    

if __name__ == '__main__':
    main()
