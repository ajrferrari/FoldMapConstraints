import sys, os, json
import numpy as np
import time
import random
import tempfile
import string
import argparse

parser = argparse.ArgumentParser(description='Generate constraints from npz file')
parser.add_argument('-fa', '--fasta', type=str, required=True, help='Fasta file')
parser.add_argument('-npz', '--npz', type=str, required=False, help='NPZ file with encapsulated constraints')
parser.add_argument('-gen', '--gen_rst', action='store_true', default=False, help='Generate constraints in tmp folder')
parser.add_argument('--rst_file', type=str, required=False, help='Constraints list file containing output from gen_rst')
parser.add_argument('-p', '--pcut', type=float, required=False, default=0.15, help='Constraints probability cutoff')
parser.add_argument('-write', '--write_rst', action='store_true', default=False, help='Write constraint file')
args = parser.parse_args()

def read_fasta(file):
    fasta=""
    with open(file, "r") as f:
        for line in f:
            if(line[0] == ">"):
                continue
            else:
                line=line.rstrip()
                fasta = fasta + line;
    return fasta
	
def gen_rst(npz_file, fasta_file):
    
    tmpdir = ''
    for _ in range(10):
        tmpdir += random.choice(string.ascii_letters)

    tmpdir = './tmp/' + tmpdir

    os.mkdir('tmp/')
    os.mkdir(tmpdir)
    
    print('temp folder:  ', tmpdir)
    
    npz = np.load(npz_file)
    
    seq = read_fasta(fasta_file)

    dist,omega,theta,phi = npz['dist'],npz['omega'],npz['theta'],npz['phi']

    # dictionary to store Rosetta restraints
    rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [], 'rep' : []}

    ########################################################
    # assign parameters
    ########################################################
    PCUT  = 0.05 # params['PCUT']
    PCUT1 = 0.5  # params['PCUT1']
    EBASE = -0.5 # params['EBASE']
    EREP  = [10.0, 3.0, 0.5] # params['EREP']
    DREP  = [0.0, 2.0, 3.5] # params['DREP']
    PREP  = 0.1 # params['PREP']
    SIGD  = 10.0 # params['SIGD']
    SIGM  = 1.0 # params['SIGM']
    MEFF  = 0.0001 # params['MEFF']
    DCUT  = 19.5 # params['DCUT']
    ALPHA = 1.57 # params['ALPHA']

    DSTEP = 0.5 # params['DSTEP']
    ASTEP = np.deg2rad(15.0) #np.deg2rad(params['ASTEP'])

    # seq = params['seq']

    
    ########################################################
    # dist: 0..20A
    ########################################################
    nres = dist.shape[0]
    bins = np.array([4.25+DSTEP*i for i in range(32)])
    prob = np.sum(dist[:,:,5:], axis=-1)
    bkgr = np.array((bins/DCUT)**ALPHA)
    attr = -np.log((dist[:,:,5:]+MEFF)/(dist[:,:,-1][:,:,None]*bkgr[None,None,:]))+EBASE
    repul = np.maximum(attr[:,:,0],np.zeros((nres,nres)))[:,:,None]+np.array(EREP)[None,None,:]
    dist = np.concatenate([repul,attr], axis=-1)
    bins = np.concatenate([DREP,bins])
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    nbins = 35
    step = 0.5
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir+"/%d.%d.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(dist[a,b])+'\n')
                f.close()
            rst_line = 'AtomPair %s %d %s %d SPLINE TAG %s 1.0 %.3f %.5f'%('CB',a+1,'CB',b+1,name,1.0,step)
            rst['dist'].append([int(a),int(b),float(p),rst_line])
    print("dist restraints:  %d"%(len(rst['dist'])))


    ########################################################
    # omega: -pi..pi
    ########################################################
    nbins = omega.shape[2]-1+4
    bins = np.linspace(-np.pi-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(omega[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    omega = -np.log((omega+MEFF)/(omega[:,:,-1]+MEFF)[:,:,None])
    omega = np.concatenate([omega[:,:,-2:],omega[:,:,1:],omega[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b>a:
            name=tmpdir+"/%d.%d_omega.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.5f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.5f'*nbins%tuple(omega[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral CA %d CB %d CB %d CA %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,b+1,name,1.0,ASTEP)
            rst['omega'].append([int(a),int(b),float(p),rst_line])
    print("omega restraints: %d"%(len(rst['omega'])))


    ########################################################
    # theta: -pi..pi
    ########################################################
    prob = np.sum(theta[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    theta = -np.log((theta+MEFF)/(theta[:,:,-1]+MEFF)[:,:,None])
    theta = np.concatenate([theta[:,:,-2:],theta[:,:,1:],theta[:,:,1:3]],axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir+"/%d.%d_theta.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(theta[a,b])+'\n')
                f.close()
            rst_line = 'Dihedral N %d CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,a+1,b+1,name,1.0,ASTEP)
            rst['theta'].append([int(a),int(b),float(p),rst_line])
            #if a==0 and b==9:
            #    with open(name,'r') as f:
            #        print(f.read())
    print("theta restraints: %d"%(len(rst['theta'])))


    ########################################################
    # phi: 0..pi
    ########################################################
    nbins = phi.shape[2]-1+4
    bins = np.linspace(-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
    prob = np.sum(phi[:,:,1:], axis=-1)
    i,j = np.where(prob>PCUT)
    prob = prob[i,j]
    phi = -np.log((phi+MEFF)/(phi[:,:,-1]+MEFF)[:,:,None])
    phi = np.concatenate([np.flip(phi[:,:,1:3],axis=-1),phi[:,:,1:],np.flip(phi[:,:,-2:],axis=-1)], axis=-1)
    for a,b,p in zip(i,j,prob):
        if b!=a:
            name=tmpdir+"/%d.%d_phi.txt"%(a+1,b+1)
            with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(phi[a,b])+'\n')
                f.close()
            rst_line = 'Angle CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,name,1.0,ASTEP)
            rst['phi'].append([int(a),int(b),float(p),rst_line])

    print("phi restraints:   %d"%(len(rst['phi'])))

    return rst
	
	
def write_rst(rst, fasta_file, pcut, sep1=1):
    
    seq = read_fasta(fasta_file)
    sep2 = len(seq)
    
    array = []
    
    array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
    array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
    array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
    array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]  
    
    random.shuffle(array)
        
    with open('constraints.rst', 'w') as f:
        for line in array:
            f.write(line+'\n')
    
    print('Generated constraints.map file with %d constraints' %(len(array)))
	
def main():
	
	start = time.time()
	
	if args.gen_rst:
		print('Generating tmp folder with constraints files ...')
		rst = gen_rst(args.npz, args.fasta)
		j = json.dumps(rst, indent=4)
		with open('rst_list.txt', 'w') as f:
			f.write(j)		 
		
	if args.write_rst and args.gen_rst:
		print('Writting constraint file ...')
		write_rst(rst, args.fasta, args.pcut)
	elif args.write_rst:
		try:
			with open(args.rst_file, 'r') as f:
				rst = json.load(f)
		except:
			print('Requires constraints list file: --rst_file option')
		
		print('Writting constraints file ... ')
		write_rst(rst, args.fasta, args.pcut)	
	
	# print('Generating tmp folder with constraints files ...')
	# rst = gen_rst(args.npz, args.fasta)
	
	# print('Writting constraint file ...')
	# write_rst(rst, args.fasta, args.pcut)
	
	print("Finished in %f minutes" % ( (time.time()-start)/60 ))
	

if __name__ == '__main__':
    main()
