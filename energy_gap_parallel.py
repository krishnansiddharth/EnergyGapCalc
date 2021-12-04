from gridData import Grid
import numpy as np
import sys
import os
import timeit
import multiprocessing

#a function to find a numbers closest value in a given array
def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


#the function averages over all potential points in a sphere
def potential(R,Ravg,V): # Ravg- information about the sidechain geometr, V-the potential grid
    X,Y,Z = R     #X,Y and Z vectors for the box
    POS = Ravg
    x,i = find_nearest(X,POS[0])
    y,j = find_nearest(Y,POS[1])
    z,k = find_nearest(Z,POS[2])
    R = [X[i],Y[j],Z[k]]
    
    V_atom = V[i][j][k]
  
    return V_atom
print(sys.argv)
state = sys.argv[1]
V_in  = sys.argv[2]
Geom_in = sys.argv[3]
Res_type= sys.argv[4]
Charged_res = sys.argv[5]
Neutral_res = sys.argv[6]
filename    = sys.argv[7]
avg_time = 20.0 #in ps


#storing file with charges 
if Res_type=="W":
	ox_charges = "/home/sk87/Programs/DNAPolymerase2021/Phi29_CTPR8_PDBs/analysisscritps/Residue_geom/Charged_rid35_0.02ns_0_0/Charge.txt"
	red_charges = "/home/sk87/Programs/DNAPolymerase2021/Phi29_CTPR8_PDBs/analysisscritps/Residue_geom/Neutral_rid35_0.02ns_0_0/Charge.txt"
elif Res_type=="Y":
	ox_charges = "Y_charges2/Charged.txt"
	red_charges = "Y_charges2/Neutral.txt"
else:
	raise ValueError("Invalid residue type")	

global Atoms,Q_red,Q_ox

OX = np.genfromtxt(ox_charges,delimiter="\t",dtype=str,skip_header=1)
Atoms_ox = OX[:,0]
Q_ox = [float(q) for q in OX[:,1]]

Red = np.genfromtxt(red_charges,delimiter="\t",dtype=str,skip_header=1)
Atoms_red = Red[:,0]
Q_red    = [float(q) for q in Red[:,1]]



#extracting 
if state == "Neutral":
	ResID = Neutral_res 
	Atoms = Atoms_red
	
else:
	ResID = Charged_res
     
	Atoms = Atoms_ox


#extracting position of the atom


def QPhi_atom(atom_num):	
	atom = Atoms[atom_num]
	q_ox = Q_ox[atom_num]
	q_red= Q_red[atom_num]
	del_q = q_ox - q_red
	POS=np.genfromtxt(Geom_in+ "/"+"Pos_rid"+str(ResID)+"_"+atom.strip()+".dat") ####
	x = POS[:,1][1:]
	y = POS[:,2][1:]
	z = POS[:,3][1:]
	t = POS[:,0][1:]
	V_ = []
	

#*****************#
	
	for i in range(len(t)-1): #iterating over all dx files	
		dxfile = V_in+"/"+filename+"_"+str(avg_time)+"_pot_"+str(int(t[i]))+"_"+str(int(t[i+1]-1))+".dx"

		g = Grid(dxfile)
		V = g.grid
		xi,yi,zi = g.origin
		nx,ny,nz = V.shape
		dx,dy,dz=g.delta
		X = np.linspace(xi,xi+(nx-1)*dx,nx)
		Y = np.linspace(yi,yi+(ny-1)*dy,ny)
		Z = np.linspace(zi,zi+(nz-1)*dz,nz)

		pos = [x[i],y[i],z[i]]
		pot = potential([X,Y,Z],pos,V) #averaged potential
		V_.append(pot)
		
			#f.write(str(V_[i])+"\t"+str(t[i])+"\n")
	#f.close()
	X_atom = np.multiply(V_,del_q)
	return X_atom 

#extracting frames data	
POS=np.genfromtxt(Geom_in+ "/"+"Pos_rid"+str(ResID)+"_"+Atoms[0].strip()+".dat") ####
t = POS[:,0][1:]
######
start = timeit.default_timer()

a_pool = multiprocessing.Pool(18)
X_atoms = np.array(a_pool.map(QPhi_atom,range(len(Atoms))))
Egap   = X_atoms.sum(axis=0)		
stop = timeit.default_timer()
print("Run Time ="+ str(stop-start))
		
f = open("X_"+state+str(ResID)+"_"+str(avg_time)+"ns_cent"+str(ResID)+"_"+str(int(t[0]))+"_"+str(int(t[-1]))+".dat","w+")
f.write("X(kT)\tframes(10ps)\n")
for i in range (len(Egap)):
	f.write(str(Egap[i])+"\t"+str(t[i])+"\n")
f.close()

