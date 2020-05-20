from functions import *

parser = argparse.ArgumentParser(description='Extract constraints from aligned pdbs')
parser.add_argument('-database', '--database', type=str, required=True, help='database directory')
parser.add_argument('-query', '--query_pdb', type=str, required=True, help='query pdb file')
args = parser.parse_args()

def main():

    start = time()

    database_directory = args.database
    reference = database_directory + 'reference.pdb'
    query = args.query_pdb

    database_files = glob.glob(database_directory + '*pos')

    print('Writting constraints file ...')

    for pos_file in database_files:
        write_to_xl_file(reference, query, pos_file)


    print("Finished in %f minutes" % ( (time()-start)/60 ))

    

if __name__ == '__main__':
    main()
