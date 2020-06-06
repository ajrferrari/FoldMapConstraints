''' Functions to run extract database and extract constraints '''

import numpy as np
import math
from time import time
import os
import argparse
from itertools import combinations
import glob as glob


def read_pdbs_files(ref, pdbs_list):
    ''' Read pdbs and dump them into a list'''
    files_list = []
    files_list.append(ref)
    with open(pdbs_list, 'r') as f:
        lines = f.readlines()
        for line in lines:
            files_list.append(line[:-1])
    return files_list

def pdb_length(pdbfile):
    ''' Get number of residues from pdbfile'''
    import subprocess
    return int(subprocess.check_output(['tail', '-2', pdbfile])[23:26])


def generate_combinations(pdbfile):
    """ Generate combinations of pairs of numbers from 1 to the protein size """
    return [x for x in combinations(range(1,pdb_length(pdbfile)+1), 2)]


def get_coord_dict_from_pdbfile(pdbfile):
    """ Create dictionary from pdbfile"""
    pdb_dict = {}
    for residue_number in range(1, pdb_length(pdbfile)+1):
        pdb_dict[residue_number] = {}
    with open('%s'%pdbfile) as f:
        for line in f:
            if 'TER' not in line and 'ATOM' in line:
                line = line.rstrip()
                atom_number = line[7:12]
                atom_name = str(line[13:16].replace(" ",""))
                residue_number = int(line[22:26])
                x_coor, y_coor, z_coor = (float(line[32:38]), float(line[40:46]), float(line[48:54]))
                pdb_dict[residue_number][atom_name] = (x_coor, y_coor, z_coor)
    return pdb_dict

def get_res_name(pdbfile, residue_number):
    """ Get residue name"""
    with open('%s'%pdbfile) as f:
        for line in f:
            if 'TER' not in line and 'ATOM' in line:
                if int(line[22:26]) == int(residue_number):
                    return str(line[17:20])


def get_map_res_num_to_res_name(pdbfile):
    """ Create dictionary of residues mapping residue number to residue name from pdbfile"""
    res_number_to_residue_name_dict = {}
    for residue_number in range(1, pdb_length(pdbfile)+1):
        res_number_to_residue_name_dict[residue_number] = get_res_name(pdbfile, residue_number)
    return res_number_to_residue_name_dict


def check_if_residue_is_not_GLY(i, map_res_num_to_res_name):
    if map_res_num_to_res_name[i] !='GLY' :
        return True
    return False


def get_res_coord(pdbdict, residue_number, atom_name):
    """ Get x_coor, y_coor, z_coor of a specific atom, residue and pdbfile"""
    res_coord = pdbdict[residue_number][atom_name]                  
    return np.array(res_coord)


def calc_distance(p1, p2):
    """ Calculate Euclidean distance between two coordinates and return it as a numpy array"""
    return np.linalg.norm(p1-p2)


def list_of_Euclidean_distances(pdbfile, max_distance=20, residue_gap=8):
    ''' Get list of residue pairs and Euclidean distances with Euclidean distance shorter than max_distance=20A'''
    distances_CB_CB = []
    pairs_residues = generate_combinations(pdbfile) # Get residue pairs
    
    pdbdict = get_coord_dict_from_pdbfile(pdbfile) # Get mapping to coordinates
    
    map_res_num_to_res_name = get_map_res_num_to_res_name(pdbfile) # Get mapping to residue name
    
    for i, j in pairs_residues:
        if abs(i-j) >= residue_gap:
            if check_if_residue_is_not_GLY(i, map_res_num_to_res_name) and check_if_residue_is_not_GLY(j, map_res_num_to_res_name):
                p1 = get_res_coord(pdbdict, i, 'CB')
                p2 = get_res_coord(pdbdict, j, 'CB')
                distance = calc_distance(p1,p2)
                if distance <= max_distance:
                    distances_CB_CB.append([i, j, distance])
    return distances_CB_CB


def calc_q_vectors(p1, p2, p3, p4):
    """ Function to calculate q vectors"""
    q1 = np.subtract(p2,p1)
    q2 = np.subtract(p3,p2)
    q3 = np.subtract(p4,p3)
    return q1, q2, q3

def calc_cross_vectors(q1,q2,q3):
    """  Function to calculate cross vectors """
    q1_x_q2 = np.cross(q1,q2)
    q2_x_q3 = np.cross(q2,q3)
    return q1_x_q2, q2_x_q3


def calc_normals(q1_x_q2, q2_x_q3):
    """ Function to calculate normal vectors to planes"""
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))
    return n1, n2

def calc_orthogonal_unit_vectors(n2,q2):
    """ Function to calculate orthogocal unit vectors"""
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2,q2)))
    u2 = np.cross(u3,u1)
    return u1, u2, u3

def calc_dihedral_angle(n1,u1,u2,u3):
    """ Function to calculate dihedral angle"""
    cos_theta = np.dot(n1,u1)
    sin_theta = np.dot(n1,u2)
    theta = -math.atan2(sin_theta,cos_theta)
    return theta


def list_of_dihedrals(list_of_Euclidean_distances, pdbfile):
    """ Return a list of calculated dihedral angle for residue pairs found with Euclidean distance shorter than max_distance=20A"""
    
    dihe_CA_CB_CB_CA = []
    
    pdbdict = get_coord_dict_from_pdbfile(pdbfile)
    
    for i in range(len(list_of_Euclidean_distances)):
        p1 = get_res_coord(pdbdict, list_of_Euclidean_distances[i][0], 'CA')
        p2 = get_res_coord(pdbdict, list_of_Euclidean_distances[i][0], 'CB')
        p3 = get_res_coord(pdbdict, list_of_Euclidean_distances[i][1], 'CB')
        p4 = get_res_coord(pdbdict, list_of_Euclidean_distances[i][1], 'CA')

        q1, q2, q3 = calc_q_vectors(p1,p2,p3,p4)
        q1_x_q2, q2_x_q3 = calc_cross_vectors(q1, q2, q3)
        n1, n2 = calc_normals(q1_x_q2, q2_x_q3)
        u1, u2, u3 = calc_orthogonal_unit_vectors(n2, q2)
        theta = calc_dihedral_angle(n1, u1, u2, u3)
        dihe_CA_CB_CB_CA.append([list_of_Euclidean_distances[i][:2][0], list_of_Euclidean_distances[i][:2][1], theta])
    return dihe_CA_CB_CB_CA



def get_closest_residue(i, pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length, dtol=3.0):
    """ Find the jth closest residue from the ith residue from the reference pdb """
    d = 10.
    res_i = 0
     
    coord_ref = get_res_coord(pdb_ref_dict, i, 'CB')
                           
    for j in range(1, length+1):
        if check_if_residue_is_not_GLY(j, map_res_num_to_res_name):
            try:
                coord_query = get_res_coord(pdb_query_dict, j, 'CB') 
            except:
                continue
            d_tmp = calc_distance(coord_ref, coord_query)
            if d_tmp < d:
                d = d_tmp
                res_i = j
        if d < 1.0:
            return res_i
    if d < dtol:
        return res_i
    else:
        return None 



def get_dictionaries_of_distances_and_dihedrals(list_of_positions, pdblist, reference):
    """ Create two dictionaries: Euclidean distances and dihedrals. list_of_positions: Dihedrals_angles or Eulidean_distances
    computed above. """


    distance_dict = {}
    dihedral_dict = {}

    for pos in range(len(list_of_positions)):

            distance_dict[pos] = []
            dihedral_dict[pos] = []

    for pdb in pdblist:


        pdb_ref_dict = get_coord_dict_from_pdbfile(reference)
        pdb_query_dict = get_coord_dict_from_pdbfile(pdb)
        length = pdb_length(pdb)
        map_res_num_to_res_name = get_map_res_num_to_res_name(pdb)

        for pos in range(len(list_of_positions)):

            res1 = get_closest_residue(list_of_positions[pos][0], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length)
            res2 = get_closest_residue(list_of_positions[pos][1], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length)
            if res1 == None or res2 == None:
                distance_dict[pos].append(None) 
                dihedral_dict[pos].append(None)	
            else:
                try:
                    # Calculate Euclidean distance
                    p1 = get_res_coord(pdb_query_dict, res1, 'CB')
                    p2 = get_res_coord(pdb_query_dict, res2, 'CB')
                    distance_dict[pos].append(calc_distance(p1,p2))
                except:
                    distance_dict[pos].append(None)
                    
                try:
                    # Calculate theta dihedral angle
                    p1 = get_res_coord(pdb_query_dict, res1, 'CA')
                    p2 = get_res_coord(pdb_query_dict, res1, 'CB')
                    p3 = get_res_coord(pdb_query_dict, res2, 'CB')
                    p4 = get_res_coord(pdb_query_dict, res2, 'CA')

                    q1, q2, q3 = calc_q_vectors(p1,p2,p3,p4)
                    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1, q2, q3)
                    n1, n2 = calc_normals(q1_x_q2, q2_x_q3)
                    u1, u2, u3 = calc_orthogonal_unit_vectors(n2, q2)
                    theta = calc_dihedral_angle(n1, u1, u2, u3)

                    dihedral_dict[pos].append(theta)
                    
                except:
                    dihedral_dict[pos].append(None)
                
    return distance_dict, dihedral_dict



def write_database(distance_dict, dihedral_dict):
    ''' Write database files containing reference pdb positions of all conserved distances and dihedrals'''
    try: 
        os.mkdir('database')
    except:
        print(' ... database directory found')
    
    distance = 1    
    while distance <= 20:
        # Get list of positions
        positions_to_Euclidean = []
        avarage_value_Euclidean = []
        for key in distance_dict:
            if distance_dict[key].count(None) == 0:
                avarage = sum(distance_dict[key])/len(distance_dict[key])
                minimum = min(distance_dict[key])
                maximum = max(distance_dict[key])
                if (maximum - minimum > distance-1) and (maximum - minimum <= distance)  :
                    positions_to_Euclidean.append(key)
                    avarage_value_Euclidean.append(avarage)
        with open('database/' + 'Euclidean_' + str(distance-1) + '_' + str(distance) + '.pos', 'w') as f:
            f.writelines("%s\n" % line for line in positions_to_Euclidean)
        with open('database/' + 'Euclidean_' + str(distance-1) + '_' + str(distance) + '.avg', 'w') as f:
            f.writelines("%s\n" % line for line in avarage_value_Euclidean)
        distance += 1
    
    angle = 0.25
    while angle <= 1:
        #Get list of dihedrals
        positions_to_Dihedrals = []
        avarage_value_Dihedrals = []
        for key in dihedral_dict:
            if dihedral_dict[key].count(None) == 0:
                angles = np.array(dihedral_dict[key])
                angles = (2*np.pi + angles ) * (angles < 0) + (angles) * (angles >0)
                avarage = sum(angles)/len(angles)
                if avarage > np.pi:
                    avarage = avarage - 2*np.pi
                minimum = min(angles)
                maximum = max(angles)
                if (maximum - minimum > angle - 0.25) and (maximum - minimum <= angle):
                    positions_to_Dihedrals.append(key)
                    avarage_value_Dihedrals.append(avarage)
        with open('database/' + 'Dihedral_' + str(angle-0.25) + '_' + str(angle) + '.pos', 'w') as f:
            f.writelines("%s\n" % line for line in positions_to_Dihedrals)
        with open('database/' + 'Dihedral_' + str(angle-0.25) + '_' + str(angle) + '.avg', 'w') as f:
            f.writelines("%s\n" % line for line in avarage_value_Dihedrals)
        angle+=0.25


def read_pos_and_avg_files(pos_file):
    
    if 'Euclidean' in pos_file:
        positions_to_Euclidean = []
        avarage_value_Euclidean = []

        with open(pos_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line[:-1]
                positions_to_Euclidean.append(int(line))
        with open(pos_file.replace('pos', 'avg'), 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line[:-1]
                avarage_value_Euclidean.append(float(line))

        return positions_to_Euclidean, avarage_value_Euclidean

    elif 'Dihedral' in pos_file: 
        positions_to_Dihedrals = []
        avarage_value_Dihedrals = []

        with open(pos_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line[:-1]
                positions_to_Dihedrals.append(int(line))
        with open(pos_file.replace('pos', 'avg'), 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line[:-1]
                avarage_value_Dihedrals.append(float(line))

        return positions_to_Dihedrals, avarage_value_Dihedrals
    else:
        return None
		


def write_to_xl_file(ref, query, pos_file, func_type, dtol=3.0):
    """ 1. Get list of positions that satisfy conditions of ddistance and ddihedral 
        2. Map closest residue of each position pair from reference residue pair from previous list
        3. Write xl file """
    
     
    pdb_ref_dict = get_coord_dict_from_pdbfile(ref)
    pdb_query_dict = get_coord_dict_from_pdbfile(query)
    length = pdb_length(query)
    map_res_num_to_res_name = get_map_res_num_to_res_name(query)
    

    Euclidean_distances = list_of_Euclidean_distances(ref)
    Dihedral_angles = list_of_dihedrals(Euclidean_distances, ref)
    
    positions, avarage_value = read_pos_and_avg_files(pos_file)
    
    tol = float(pos_file.split('/')[-1].split('_')[-1].replace('.pos', ''))/2

    if func_type != 'LINEAR_PENALTY' and func_type != 'FLAT_HARMONIC' and func_type != 'USOG':
        print('Function not recognized')
        exit()


    with open('constraints.map', 'w') as f:
        
        if 'Euclidean' in pos_file:

            i=0
            for pos in positions:
                res1 = get_closest_residue(Euclidean_distances[pos][0], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length, dtol)
                res2 = get_closest_residue(Euclidean_distances[pos][1], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length, dtol)
                if (res1 != None and res2 != None) and (abs(res1-res2) >= 10):
                    if func_type == 'LINEAR_PENALTY':
                        f.write("AtomPair CB %i CB %i LINEAR_PENALTY %f -1.0 %f 1.0\n" % (res1, res2, avarage_value[i], tol ))
                    elif func_type == 'FLAT_HARMONIC': 
                        f.write("AtomPair CB %i CB %i FLAT_HARMONIC %f 1.0 %f\n" % (res1, res2, avarage_value[i], tol ))
                    elif func_type == 'USOG':
                        f.write("AtomPair CB %i CB %i USOGFUNC 1 %f 1.0 1.0 \n" % (res1, res2, avarage_value[i] ))
                i+=1

        elif 'Dihedral' in pos_file:
            i=0        
            for pos in positions:
                res1 = get_closest_residue(Euclidean_distances[pos][0], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length)
                res2 = get_closest_residue(Euclidean_distances[pos][1], pdb_ref_dict, pdb_query_dict, map_res_num_to_res_name, length)				
                if (res1 != None and res2 != None) and (abs(res1-res2) >= 10):
                    #if func_type == 'LINEAR_PENALTY':
                     #   f.write("Dihedral CA %i CB %i CB %i CA %i LINEAR_PENALTY %f -1.0 %f 1.0 \n" % ( res1, res1, res2, res2, avarage_value[i], tol))
                    #elif func_type == 'FLAT_HARMONIC':
                    #    f.write("Dihedral CA %i CB %i CB %i CA %i FLAT_HARMONIC %f 1.0 %f\n" % ( res1, res1, res2, res2, avarage_value[i], tol))
                    #elif func_type == 'USOG':
                    #    f.write("Dihedral CA %i CB %i CB %i CA %i USOGFUNC 1 %f 1.0 1.0 \n" % ( res1, res1, res2, res2, avarage_value[i] ))
                    if tol == 0.25:
                        f.write("Dihedral CA %i CB %i CB %i CA %i CIRCULARHARMONIC %f 0.5 \n" % ( res1, res1, res2, res2, avarage_value[i] ))
                    if tol == 0.50:
                        f.write("Dihedral CA %i CB %i CB %i CA %i CIRCULARHARMONIC %f 0.75 \n" % ( res1, res1, res2, res2, avarage_value[i] ))
                    if tol == 0.75:
                        f.write("Dihedral CA %i CB %i CB %i CA %i CIRCULARHARMONIC %f 1.0 \n" % ( res1, res1, res2, res2, avarage_value[i] ))
                    if tol == 1.0:
                        f.write("Dihedral CA %i CB %i CB %i CA %i CIRCULARHARMONIC %f 1.25" % ( res1, res1, res2, res2, avarage_value[i] ))

                i+=1

