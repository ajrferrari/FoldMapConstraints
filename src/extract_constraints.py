from functions import *

parser = argparse.ArgumentParser(description='Extract constraints from aligned pdbs')
parser.add_argument('-database', '--database', type=str, required=True, help='database directory')
parser.add_argument('-query', '--query_pdb', type=str, required=True, help='query pdb file')
parser.add_argument('-func', '--func', type=str, required=True, help='Function to describe constrainsts: FLAT_HARMONIC and LINEAR_PENALTY available')
args = parser.parse_args()

def main():

    start = time()

    database_directory = args.database
    reference = database_directory + 'reference.pdb'
    query = args.query_pdb
    func = args.func

    database_files = glob.glob(database_directory + '*pos')

    print('Writting constraints file ...')

    for pos_file in database_files:
        write_to_xl_file(reference, query, pos_file, func)


    print("Finished in %f minutes" % ( (time()-start)/60 ))

    

if __name__ == '__main__':
    main()
