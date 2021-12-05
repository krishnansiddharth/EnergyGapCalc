import numpy as np
from scipy.fft import fft, fftfreq
import os
import fnmatch

def X_time(file,ddcfreq,Dt):#dt in femtoseconds
    dt = Dt*1e-6
    filename  = file
    data      =  np.genfromtxt(filename)
    V =  data[1:,0] *0.0259#converting to eV
    frame_data      =  data[1:,1]
    t = ddcfreq*frame_data*dt #in nanosecond
    return t,V

def cum_avg(X):
    cumsum = np.cumsum(X)
    denominator = np.arange(1,len(X)+1)
    return cumsum/denominator

def auto_corr(X):
    delX = X - cum_avg(X)
    corr = np.correlate(delX,delX,mode="full")
    return corr[int(corr.size/2):]/corr[int(corr.size/2)]

def Energy_gap(acceptor,donor):
    ta,Xa = append_files(acceptor)
    td,Xd = append_files(donor)
    
    if len(ta)==len(td):
        X     = np.array(Xa)-np.array(Xd)
    
        return ta,X
    else:
        print( "Error: Different Trajectories for donor and acceptor")
        return NaN

def append_files(path):
    dirname = path
    dirname_split= dirname.split("_")
    for i in dirname_split: 
        if 'ps' in i : 
            sampling_time = float(i.strip('ps'))
            time_per_step = 4
            sampling_steps = int(sampling_time/4*1000)

        else:
            continue
   
    t = np.array([])
    X = np.array([])
    
    for files in os.listdir(path):
        if files.endswith('.dat'):
            t_ox,X_ox  = X_time(path+"/"+files,sampling_steps,time_per_step)

            
            t=np.append(t,t_ox)
            X=np.append(X,X_ox)
        else:
            continue    
    t_sorted,X_sorted = (list(i) for i in zip(*sorted(zip(t,X)))) #linking the 2 lists and sorting wrt to Pair
    return t_sorted,X_sorted

def lambda_var(X):
    return X.std()**2/(2*0.0259)

def lambda_st(X1,X2):
    return np.abs(X1.mean()-X2.mean())/2
def reorg(X1,X2):
    l_st = lambda_st(X1,X2)
    l_var= (lambda_var(X1)+lambda_var(X2))/2
    return l_st**2/l_var

def Activation(X,dG):
    lambda_ = lambda_var(X)
    X.mean
    F = X.mean()**2/(4*lambda_) +dG
    return F

def Activation2(X1,X2,dG):
    st = lambda_st(X1,X2)
    lbda = lambda_var(X1)
    X1mean = st + lbda/st*dG
    X2mean = st - lbda/st*dG   
    F1     = X1mean**2/(4*lambda_var(X1))
    F2     = X2mean**2/(4*lambda_var(X2))
    return F1,F2

def Probability(X,F_act):
    lambda_ = lambda_var(X) 
    P     = 1/(np.sqrt(4*np.pi*lambda_*0.0259))*np.exp(-F_act/0.0259)
    return P
    
def Coupling(donor,acceptor,Ree):
    N = {"Tyr":7,"Trp":9}
    Vda  =  2.7/np.sqrt(N[donor]*N[acceptor])*np.exp(-0.72*Ree)
    print("Prefactor for V:"+str(2.7/np.sqrt(N[donor]*N[acceptor])))
    

    return Vda

def Coupling_dutton(donor,acceptor,f,Ree):
    N = {"Tyr":7,"Trp":9}
    bvac = 2.8/2
    batm = 0.9/2
    Vda  =  2.7/np.sqrt(N[donor]*N[acceptor])*np.exp(-1*(Ree-3.6)*(batm*f+bvac*(1-f)))
    return Vda

def Rate(donor,acceptor,X1,X2,Ree,dG,dir_):
    hbar = 6.58e-16
    F1,F2=Activation2(X1,X2,dG)
    if dir_=="forward":
        P = Probability(X1,F1)
    else:
        P = Probability(X2,F2)
    print("Probability of crossing barrier:"+str(P))
    V = Coupling(donor,acceptor,Ree)
    k = 2*np.pi/hbar*V**2*P
    print("Prefactor for Golden Rule:"+str(2*np.pi/hbar))
    return k

def Rate_dutton(X1,X2,Ree,dG,dir_):
    lambda_ = reorg(X1,X2)
    F1,F2=Activation2(X1,X2,dG)
    rho = 0.76
    if dir_=="forward" :
        P = Probability(X1,F1)
    else:
        P = Probability(X2,F2)
    
    K = 10**(13 - (1.2 - 0.8*rho)*(Ree-3.6))*P
   
    return K

def find_nearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def Rate_dist_dutton(Ree,lambda_,F1,F2,dir_):
    rho = 0.76
    if dir_=="forward":
        P     = 1/(np.sqrt(4*np.pi*lambda_*0.0259))*np.exp(-F1/0.0259)
    else:
        P     = 1/(np.sqrt(4*np.pi*lambda_*0.0259))*np.exp(-F2/0.0259)
    
    K = 10**(13 - (1.2 - 0.8*rho)*(Ree-3.6))*P
   
    return K

def Reorganization(filename,donor,acceptor):
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Neutral*'+str(acceptor)+'*'+str(donor)+'*'):
            
            X_red1=[]
            X_red1.append(os.path.join(filename,file))
            
           
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Charged*'+str(donor)+'*'):
            X_ox1=[]
            
            X_ox1.append(os.path.join(filename,file)) 
            
   
    
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Neutral*'+str(donor)+'*'+str(acceptor)+'*'):
            
            X_red2=[]
            X_red2.append(os.path.join(filename,file))
            
            
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Charged*'+str(acceptor)+'*'):
            X_ox2=[]
            X_ox2.append(os.path.join(filename,file))
            
      
    
    try:
        t1,X1 =   Energy_gap (X_red1[0],X_ox1[0])       
        t2,X2 =   Energy_gap (X_red2[0],X_ox2[0])
        X2 = -X2
       # dG = (-(X1.mean()/X1.std())**2/2*np.log(1/(2*np.pi*X2.std())) + (X2.mean()/X2.std())**2/(2)*np.log(1/(2*np.pi*X2.std())))*0.026
        dG = -((X1.mean()**2)/(2*X1.std()**2) + (X2.mean()**2)/(2*X2.std()**2))*0.026
        lam =   reorg(X1,X2)
        stokes =   lambda_st(X1,X2)
        F1,F2  =   Activation2(X1,X2,dG)
        return lam,dG,stokes
    except:
        return str(donor)+"->"+str(acceptor)+"Transition not available"

def egap(filename,acceptor,donor):
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Neutral*'+str(acceptor)+'*'+str(donor)+'*'):
            
            X_red1=[]
            X_red1.append(os.path.join(filename,file))
            print("Xred1"+str(X_red1))
           
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Charged*'+str(donor)+'*'):
            X_ox1=[]
            
            X_ox1.append(os.path.join(filename,file)) 
            print("X_ox1"+str(X_ox1))
   
    
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Neutral*'+str(donor)+'*'+str(acceptor)+'*'):
            
            X_red2=[]
            X_red2.append(os.path.join(filename,file))
            print("Xred2:"+str(X_red2))
            
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Charged*'+str(acceptor)+'*'):
            X_ox2=[]
            X_ox2.append(os.path.join(filename,file))
            print("X_ox2:"+str(X_ox2))
    try:
        t1,X1 =   Energy_gap (X_red1[0],X_ox1[0])   
        
        t2,X2 =   Energy_gap (X_ox2[0],X_red2[0])
        
        return X1,t1,X2,t2
    except:
        return str(donor)+"->"+str(acceptor)+"Transition not available"


def plot_X_stats_charged(filename,residue,tN,ti=0):
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Charged*'+str(residue)+'*'):
            X_ox1=[]
            
            X_ox1.append(os.path.join(filename,file)) 
            print("X_ox1"+str(X_ox1))
    t,X = functions.append_files(X_ox1[0])
    
    P1,Edges = np.histogram(X[ti:tN],bins=25,density=True)
    mid_points = []
    for i in range(len(Edges)-1):
        avg = (Edges[i]+Edges[i+1])/2
        mid_points.append(avg)
        
    return P1,mid_points

def plot_X_stats_neutral(filename,residue,host,tN,ti=0):
    for file in os.listdir(filename):
        if fnmatch.fnmatch(str(file), '*_Neutral*'+str(residue)+'*'+str(host)+'*'):
            
            X_red=[]
            X_red.append(os.path.join(filename,file))
            print("Xred"+str(X_red))
    t,X = functions.append_files(X_red[0])
    P1,Edges = np.histogram(X[ti:tN],bins=25,density=True)
    mid_points = []
    for i in range(len(Edges)-1):
        avg = (Edges[i]+Edges[i+1])/2
        mid_points.append(avg)
        
    return P1,mid_points

def plot_X_stats(filename,donor,acceptor,N):
    X1,t1,X2,t2 = egap(filename,acceptor,donor)
    P1,Edges1 = np.histogram(X1[:N],bins=25,density=True)
    P2,Edges2 = np.histogram(X2[:N],bins=25,density=True)
    mid_points1 = []
    mid_points2 =  []
    for i in range(len(Edges1)-1):
        avg1 = (Edges1[i]+Edges1[i+1])/2
        mid_points1.append(avg1)
        avg2 = (Edges2[i]+Edges2[i+1])/2
        mid_points2.append(avg2)
    dg = ((X1.mean()/X1.std())**2 - (X1.mean()/X1.std())**2)*0.026
    stats = [X1[0:int(len(X1)/2)].mean(),X1[0:int(len(X1)/2)].std(),X2[0:int(len(X2)/2)].mean(),X2[0:int(len(X2)/2)].std()]
    return P1,mid_points1,P2,mid_points2