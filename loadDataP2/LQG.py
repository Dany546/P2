"""
Created on Wed Oct 13 07:32:37 2021

Link to our images: https://drive.google.com/open?id=1D7sLIqdaSxhrZuw4FtkXCrE5e351-Erc&authuser=rimez.dany%40gmail.com&usp=drive_fs

@author: Rimez Dany and Loffet Alexandre
"""   
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from collections import deque
import numpy.linalg as lin
from PIL import Image
import numpy as np
import imageio
import time
import os
 
last = None
class LQG:
    # defining core variables
    a1 = 0.025 + 0.045 + 0.33**2
    a2 = 0.3*0.16 ; a3 = 0.045
    arm, forearm = (30 , 33)
    tau = 0.06
    b = 0.5*np.eye(2) + 0.025*np.eye(2)[:,::-1]  
    
    def __init__(self,eval_point=np.array([np.pi/4,np.pi/2,0,0,0,0]),time_range=0.6,dt=0.01,w=None,r=None):
        """
        Creates an LQG object

        Parameters
        ----------
        eval_point : equilibrium point for linearisation. The default is [45°,90°,0,0,0,0].
        time_range : time range for simulation. The default is 0.6
        dt : time step for simulation. The default is 0.01.
        w : Weights for the Q matrix. The default is [100,10,10] (for positions, speed and torques respectively).
        r : Weights for the R matrix. The default is [1e-4,1e-4]. 
        """ 
        self.time_range = time_range
        self.dt = dt ; self.ns = 12
        self.N = int(self.time_range/self.dt) + 10
        self.eq = eval_point  
        self.states = np.expand_dims(np.append(eval_point,eval_point),axis=1)
        self.states = np.tile(self.states,(1,10))  
        self.L = np.tile(np.zeros((1,2,self.ns)),(10,1,1))  
        self.w = w
        self.R = np.diag(r) if np.any(r!=None) else np.eye(2)*1e-4
        self.compute_dynamics_simple() 
                    
    def compute_dynamics_simple(self):  
        # checking eq_point
        ev_point = self.eq
        if len(ev_point)<6:
            print("there's something wrong here, evaluation point should be of dimension 6")
        else:    
            ns = self.ns  
            # initialising linear dynamics matrixes
            A_loc = np.zeros((ns,ns))
            # naming important variables
            theta1, theta2, dtheta1, dtheta2, tau1, tau2 = ev_point[:] 
            
            M = np.ones((2,2))*(self.a3)  
            M[1,1] = self.a1   
            A_loc[2:4,2:4] = np.linalg.inv(M)@(-self.b)  
            A_loc[2:4,4:6] = np.linalg.inv(M)   
            A_loc[:2,2:4] = np.eye(2)  # vu que dtheta = dtheta (:
            A_loc[4:6,4:6] = np.eye(2)*(-1/self.tau) 
            
            B_loc = np.zeros((ns,2))
            B_loc[4:6,:] = np.eye(2)*(1/self.tau) 
               
            self.Ad = (A_loc * self.dt) + np.eye(ns)
            self.Bd = (B_loc* self.dt)  
        
    def compute_target_angles_from_position(self,positions):  
        """ 
        Parameters
        ----------
        positions : positions in cartesian coordinates

        Returns
        -------
        theta1 : shoulder angles corresponding to the given positions
        theta2 : elbow angles corresponding to the given positions  
        """
        # with cosine law, we find for right arm: 
            
        arg = lin.norm(positions,axis=0)
        assert not np.any(arg > self.arm + self.forearm), "your arm is too short to reach the target !" 
        theta2 = np.pi - np.arccos( ( (self.arm**2) - (arg**2) + self.forearm**2) / (2*(self.arm)*(self.forearm)) )   
        phi = np.arccos(positions[0]/arg)  
        theta1 = phi - np.arccos( ( (arg**2) + (self.arm**2) - (self.forearm**2)) / (2*(self.arm)*arg) )  
        return theta1, theta2 
         
    def simulation(self,targets,simul_type="LQG",visu_pert=None,mech_pert=None,delay=0,plot=True,filename=None,hpar=None,notangles=False,angout=False): 
        self.visu = visu_pert[0] if np.any(visu_pert!=None) else None
        self.mech = mech_pert[0] if np.any(mech_pert!=None) else None 
        self.hpar = hpar if np.any(hpar!=None) else None 
        if delay != 0:
            self.dstate = np.expand_dims(np.tile(self.eq,2*(delay+1)),axis=1)
            self.dstate = np.tile(self.dstate,(1,10))  
            self.L = np.tile(np.zeros((1,2,self.ns*(delay+1))),(10,1,1))  
        if simul_type == "LQG": 
            self.lq = False
            self.LQ_simulation(targets,visu_pert=visu_pert,mech_pert=mech_pert,delay=delay,notangles=notangles) 
            self.states = self.dstate if delay!=0 else self.states
            self.LQG_simulation(visu_pert=visu_pert,mech_pert=mech_pert,delay=delay)
            self.states = self.dstate if delay!=0 else self.states
            print("Simulation done")
            if plot:
              self.plot_results(simul_type="LQG",filename=filename) 
              print("Done")
            if angout:
                return self.states, self.get_cart_coords(self.states[:2])[:,:,1].T, self.get_cart_coords(self.y[:2])[:,:,1].T, self.y[:2]
            return self.states, self.get_cart_coords(self.states[:2])[:,:,1].T, self.get_cart_coords(self.y[:2])[:,:,1].T
        elif simul_type == "LQ":
            self.lq = True 
            self.LQ_simulation(targets,visu_pert=visu_pert,mech_pert=mech_pert,delay=delay)
            self.states = self.dstate if delay!=0 else self.states
            print("Simulation done")
            if plot:
              self.plot_results(simul_type="LQ",filename=filename) 
              print("Done")
            return self.states, self.get_cart_coords(self.states[:2])[:,:,1].T , self.get_cart_coords(self.y[:2])[:,:,1].T
     
    def LQ_simulation(self,targets,visu_pert=None,mech_pert=None,delay=0,notangles=False): 
        ns = self.ns ; N = self.N  
        A = self.Ad ; B = self.Bd  
        
        if self.lq: 
            self.y = self.states[:self.ns,:10] 
        if delay!=0:  
                A_row,A_col = np.shape(A)
                B_row,B_col = np.shape(B)
                B_delay = np.zeros((B_row+delay*12,B_col))
                B_delay[:B_row,:B_col] = B[:,:]
                A_delay = np.eye(A_row+delay*12,A_col+delay*12,k=-12)
                A_delay[:A_row,:A_col] = A  
                 
                A, B = A_delay, B_delay 

        for mmm,target in enumerate(targets):  
            ns = (delay+1)*self.ns
            X = np.zeros((self.ns,N)) 
            Xi = np.random.normal(loc=0, scale=1e-4, size=(N, ns))
            for j in range(delay+1):   
                Xi[:,j*self.ns+self.ns//2:(j+1)*(self.ns)] = np.random.normal(loc=0, scale=1e-4, size=(N, self.ns//2))
                Xi[:,j*self.ns:(j)*(self.ns)+self.ns//2] = np.random.normal(loc=0, scale=1e-4, size=(N, self.ns//2)) 
            Xi[self.ns//2:,:] = 0
            
            # forward 
            X[:6] = self.eq[:,None] 
            target_angles = self.compute_target_angles_from_position(target) if notangles else target
            X[-6,:] = target_angles[0] ; X[-5,:] = target_angles[1]
               
            if np.all(visu_pert != None):
                pert, loc, duration = visu_pert[0][mmm],visu_pert[1][mmm],visu_pert[2][mmm] 
                loc //= self.dt
                loc -= mmm*self.N 
                duration //= self.dt 
                my_pert = self.compute_target_angles_from_position(pert) 
                X[-6,int(loc):int(loc+duration)] = my_pert[0] 
                X[-5,int(loc):int(loc+duration)] = my_pert[1]
              
            X = np.tile(X,(delay+1,1))
            L = np.zeros((N, 2, ns))
            S = np.zeros((N, ns, ns)) 
            R = self.R
            Q = np.zeros((N, ns, ns))    
            Q[-11] *= 2
            w1, w2, w3 = 100, 10, 1
            weights = [*self.w,*self.w] if np.any(self.w!=None) else [w1,w1,w2,w2,w3,w3,w1,w1,w2,w2,w3,w3]
            weights2 = [*self.w] if np.any(self.w!=None) else [w1,w1,w2,w2,w3,w3]
            QN = np.diag(weights) - np.diag(weights2,6) - np.diag(weights2,-6) 
            S[-11,:self.ns,:self.ns] = QN  
               
            for ii in range(N-11,0,-1):  
                P = B.T@S[ii]@B
                pin = np.linalg.inv(R+P)
                F = (B.T)@S[ii]@A
                L[ii-1] = pin@F
                newa = (A-B@L[ii-1])
                S[ii-1] = A.T@S[ii]@newa

            loc, dur = ((mech_pert[1][mmm]//self.dt) - mmm*self.N + 10 , (mech_pert[2][mmm]//self.dt)) if np.any(mech_pert!=None) else (0,0)
            mech = np.zeros_like(X[:,0])
            if np.any(mech_pert!=None):
                mech[4:6] = mech_pert[0][mmm] 
            
            if self.lq:  
                self.y[:,:10]  = X[:self.ns,:10] 
                s = self.ns//2 
                Y = np.zeros_like(X[-self.ns:])
                Y[:,:10] = X[:self.ns,:10]
                j = 0
                targ=np.zeros_like(X[s:self.ns,0])
                last=np.zeros_like(X[s:self.ns,0])
                for i in range(N-11):       
                    last[:] = X[s:self.ns,i]     
                    targ[:] = X[s:self.ns,i+1]
                    MM = mech if (loc<=i and i<=loc+dur) else 0 
                    X[:,i] = X[:,i] + MM if np.any(mech_pert!=None) else X[:,i]  
                    X[:,i+1] = (A-B@L[i])@X[:,i]  + Xi[i]  
                    if np.any(last!=targ) or j!=0: # if any switch caused by visual perturbation 
                        X[j*self.ns+s:(j+1)*(self.ns),i+1] = last
                        j += 1
                        if j==delay+1:
                            j = 0   
                    Y[:,i+1] = X[-self.ns:,i+1] # like H@X in LQG
                
                self.y = np.concatenate((self.y,Y[:,:-10],X[:self.ns,:10]),axis=1)
                
            if delay==0:
                self.states = np.concatenate((self.states,X), axis=1)   
            else:     
                self.dstate = np.concatenate((self.dstate,X), axis=1)  
            
            self.L = np.concatenate((self.L, L), axis=0)    
        print("LQ step done")
    
    def LQG_simulation(self,visu_pert=None,mech_pert=None,delay=0): 
        ns = self.ns ; N = self.N 
        A = self.Ad ; B = self.Bd ; X = self.states.T
        L = self.L
        # extraction of targets  
        targets = np.array([state[-self.ns//2:] for mmh, state in enumerate(X) if (mmh-10)%(N)==0])  
        self.y = self.states[:self.ns,:10] 
        nd = (delay+1)*ns 
        Z = np.zeros_like(X)
        Z[:] = X[:]
        if delay!=0: 
            # definition matrix A and B for delay
            A_row,A_col = np.shape(A)
            B_row,B_col = np.shape(B)
            B_delay = np.zeros((B_row+delay*12,B_col))
            B_delay[:B_row,:B_col] = B[:,:]
            A_delay = np.eye(A_row+delay*12,A_col+delay*12,k=-12) 
            A_delay[:A_row,:A_col] = A  
              
            A, B = A_delay, B_delay
        
        mmm0 = 0
        for index,target in enumerate(targets):  
            mmm = index  
            # Initialization of the command and observable
            Y = np.zeros((N-10, ns)) 
            U = np.zeros((N,2)) 

            # initialisation of new states
            start = index*N + 10           # there's no move for the first ten steps
            Xhat = Z[start:start+N]  

            # Initialization of the covariance matrix of the state  
            Sigma = np.zeros((N, nd, nd))
            Sigma[0,:nd,:nd] = np.eye(nd) 
             
            loc,duration = (-1,-1)
            if np.all(visu_pert != None):
                pert, loc, duration = visu_pert[0][mmm],visu_pert[1][mmm],visu_pert[2][mmm] 
                loc //= self.dt
                loc -= mmm*self.N 
                duration //= self.dt 
                my_pert = self.compute_target_angles_from_position(pert)  
                
            # Some more initialization 
            K = np.zeros((N, nd, ns))
            H = np.zeros((ns, nd))   
            H[:,-ns:] = np.eye(ns)
            useHH = False
            if np.any(self.hpar!=None):
                useHH = True
                H = np.zeros((N,ns, nd)) 
                mean_loc, scale, wich, diag = self.hpar[0], self.hpar[1], self.hpar[2], self.hpar[3] 
                hi = np.random.normal(loc=0, scale=scale, size=(N, ns, ns))  
                mask = np.zeros_like(hi) 
                if 0 in wich:
                    mask[:,:2,:2] = np.eye(2) if diag else np.ones((2,2))
                if 1 in wich:
                    mask[:,:2,10:12] = np.eye(2) if diag else np.ones((2,2))  
                if 2 in wich:
                    mask[:,10:12,10:12] = np.eye(2) if diag else np.ones((2,2))
                hi *= mask
                H[:,:,-ns:] = np.eye(ns) + np.abs(hi) + mean_loc  
                HH = H/(np.eye(ns)[None,:,:]+hi).sum(2)[:,:,None]
                # H[:,:2,10:12] = H[:,10:12,:2].T # imposing symmetry of H, not always used 
                # H[:,11,10] = H[:,10,11]       # imposing symmetry of H, not always used  
                # H[:,1,0] = H[:,0,1]           # imposing symmetry of H, not always used  
            
            s = self.ns//2   
            Xi = np.random.normal(loc=0, scale=1e-4, size=(N, nd))
            for j in range(delay+1):   
                Xi[:,j*self.ns+s:(j+1)*(self.ns)] = np.random.normal(loc=0, scale=1e-4, size=(N, s)) 
                Xi[:,j*self.ns:(j)*(self.ns)+s] = np.random.normal(loc=0, scale=1e-4, size=(N, s)) 
            Omega = np.random.normal(loc=0, scale=1e-4, size=(N, ns))  
            oXi = 0.1 * (B @ B.T)                           # from the LQG TP of Modelling of biological systems course
            oOmega = 0.1 * np.max(np.max(oXi)) * np.eye(ns) # from the LQG TP of Modelling of biological systems course
            
            loc, dur = ((mech_pert[1][mmm]/self.dt) - start , (mech_pert[2][mmm]/self.dt)) if np.any(mech_pert!=None) else (0,0)
            mech = np.zeros_like(Xhat[0])
            if np.any(mech_pert!=None):
                mech[4:6] = mech_pert[0][mmm]   
            
            sign = np.array([-1,0,1,0,0,0,0,0,0,0,0,0,0,0])[index] 
            corr = np.array([10000,100000000000000000000000,10000,1,1,1,1,1,1,1,1,1])[index]
            j=0 ; lastvx = 0 ; lastvy = 0
            targ=np.zeros_like(Z[0,s:ns])
            last=np.zeros_like(Z[0,s:ns])
            switch = -1 ; switched = False
            for i in range(N-11):    
                H = HH[i] if useHH else H 
                last[:] = X[start+i,s:ns] 
                K[i] = A@Sigma[i]@H.T@np.linalg.inv(H@Sigma[i]@H.T+ oOmega)
                Sigma[i+1] = oXi + (A-K[i]@H)@Sigma[i]@A.T  
                Y[i] = H@Z[start+i] + Omega[i] 
                mech = Z[start+i]*0
                vy = ((self.forearm*Z[start+i][2:4].sum()*np.cos(Z[start+i][:2].sum())+Z[start+i][2]*self.arm*np.cos(Z[start+i][0])))
                actu1 = self.arm*sign*0.13*vy*np.sin(Z[start+i][:2].sum())
                actu2 = -self.forearm*sign*0.13*vy*np.cos(Z[start+i][:2].sum())
                accx = (lastvx-actu1)/self.dt
                accy = (lastvy-actu2)/self.dt
                if sign == -1:
                    mech[4] = actu1/corr
                    mech[5] = actu2/corr
                else:    
                    mech[5] = actu1/corr
                    mech[4] = actu2/corr
                Z[start+i] = Z[start+i] + mech if (loc<=i and i<=loc+dur) or True else Z[start+i] 
                nextZ = ((A-B@L[start+i])@Xhat[i]) + Xi[i]
                mech = Z[start+i]*0
                vy = ((self.forearm*nextZ[2:4].sum()*np.cos(nextZ[:2].sum())+nextZ[2]*self.arm*np.cos(nextZ[0])))
                actu1 = self.arm*sign*0.13*vy*np.sin(nextZ[:2].sum())
                actu2 = -self.forearm*sign*0.13*vy*np.cos(nextZ[:2].sum()) 
                if sign == -1 or True:
                    mech[4] = actu1/corr
                    mech[5] = actu2/corr
                else:    
                    mech[5] = actu1/corr
                    mech[4] = actu2/corr
                Z[start+i][4:6] -= mech[4:6]*(self.tau) 
                Xhat[i+1] = A@Xhat[i] + B@U[i] + K[i]@(Y[i]-H@Xhat[i])  
                Z[start+i+1] = ((A-B@L[start+i])@Xhat[i]) + Xi[i]
                targ[:] = X[start+i+1,s:ns]    
                if np.any(visu_pert!=None): # if any switch caused by visual perturbation ==> this does'nt work as said in the report
                    for j in range(delay+1):  
                        Xhat[i+1,j*self.ns+s:(j+1)*(self.ns)] = Z[start+i+1,j*self.ns+s:(j+1)*(self.ns)] = last
                      
               
            Y[-1] = H@Z[start+N-11] + Omega[N-11]  
            # update of the states (don't forget the 10 first steps of no moves)
            if delay==0:
                self.states[:,start:start+N-10] = Xhat[:-10].T
                self.states[:,start+N-10:start+N] = Z[0].T[:,None]
            else:     
                self.dstate[:,start:start+N] = Xhat.T
                self.dstate[:,start+N-10:start+N] = Z[0].T[:,None] 
          
            self.y = np.concatenate((self.y,Y.T,X[:10,:self.ns].T),axis=1)
     
    def get_cart_coords(self,angles_coord=None):
        if np.any(angles_coord!=None):
            angles_coord=angles_coord
        else:
            angles_coord=self.y[:2]
        cart_coords = np.zeros((angles_coord.shape[1],2,2))
        for i,angles in enumerate(angles_coord.T):
            O1, O2 = angles    
            cart_coords[i][:,0] = self.arm * np.array([np.cos(O1), np.sin(O1)])
            cart_coords[i][:,1] = cart_coords[i][:,0] + self.forearm * np.array([np.cos(O1+O2), np.sin(O1+O2)])  
        return cart_coords
        
    def plot_results(self,simul_type="LQG",filename=None):  
        # adapted from https://matplotlib.org/stable/gallery/animation/double_pendulum.html
        print("Making you a beautiful gif :)")
         
        N = self.N 
        speed_corr = (N-10)/60
        corr = int(np.rint(speed_corr)) 
        
        coords = None
        coords = self.get_cart_coords(self.y[:2,::corr]) 
        target = self.get_cart_coords(self.states[-6:-4,::corr])[:,:,1] 
        vis = self.get_cart_coords(np.array(self.compute_target_angles_from_position(self.visu.T)))[:,:,1]  if np.any(self.visu!=None) else np.zeros_like(target)
        
        scale = 1.7
        fig = plt.figure(figsize=(5*scale, 4*scale))
        L = self.arm + self.forearm
        ax = fig.add_subplot(autoscale_on=False, xlim=(-L/2, 0.6*L), ylim=(-0.05*L, 1*L))
        ax.set_aspect('equal')
        ax.grid() 
        
        if target.size>2:
            t1, = ax.plot(target[:,0],target[:,1],'o',color='yellow',markersize=3*scale)   
            t2, = ax.plot([],[],'o',color='red',fillstyle="none", markersize=10*scale)
        else:    
            ax.plot(target[:,0],target[:,1],'o',color='red',markeredgewidth=2*scale,fillstyle="none",markersize=11*scale)
        trace, = ax.plot([], [], '*', lw=1, color="orange", markersize=4*scale)
        line, = ax.plot([], [], 'o-', lw=2*scale, markersize=8*scale)
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        history_x, history_y = deque(maxlen=N//5), deque(maxlen=N//5)
         
        
        global last
        last = None
        def animate(i):
            global last  
            thisx = np.append(0,coords[i,0,:])
            thisy = np.append(0,coords[i,1,:])
        
            if i == 0:
                history_x.clear()
                history_y.clear()
        
            if (i-10)%int(5) == 0 or (i-10)%N==0: # on dessine une étoile sur ~5
                history_x.appendleft(thisx[2])
                history_y.appendleft(thisy[2])  
            
            if np.any(last!=target[i]): 
                tg = None
                if np.any(self.visu!=None):
                    tg = np.logical_and(target[:,0]!=target[i,0],target[:,1]!=target[i,1])
                else:    
                    tg = np.logical_and(np.logical_and(target[:,0]!=target[i,0],target[:,1]!=target[i,1]),np.any(np.logical_and(target[:,0]!=vis[:,0,None],target[:,1]!=vis[:,1,None]),axis=0))
                all_but_one = np.unique(target[tg],axis=0)  
                t1.set_data(all_but_one[:,0],all_but_one[:,1])
                t2.set_data(target[i,0],target[i,1])
                last = target[i] 
              
            line.set_data(thisx, thisy)
            trace.set_data(history_x, history_y)
            
            time_text.set_text(time_template % (i*self.dt))
            return line, trace, time_text
        
        F = int(20 * (1 + speed_corr - int(speed_corr)))
        # default in 20 fps ==> we change to F 
        global anim
        myAnimation = animation.FuncAnimation(fig, animate, len(coords), interval=1000*N/(speed_corr*F/self.dt), blit=True)  
        anim = myAnimation
        if filename!=None:
          with open(filename+'.gif',"w") as file:
            myAnimation.save(filename=filename+'.gif', writer='imagemagick', fps=F)
            file.close()
        # os.system("mygif.gif") 
        
    def draw_circle(self): 
        angles = np.linspace(0,2*np.pi,self.N) 
        c1 = self.arm * np.cos(np.pi/4) + self.forearm * np.cos(3*np.pi/4)
        c2 = self.arm * np.sin(np.pi/4) + self.forearm * np.sin(3*np.pi/4)
        radius = 12
        pos = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)])  
        self.states[0,:], self.states[1,:] = self.compute_target_angles_from_position(pos)
        self.plot_results() 

if __name__ == "__main__" :        
            
    from itertools import product
    
    def traverse(d): # inspired from https://stackoverflow.com/questions/41523543/how-can-i-create-a-list-of-possible-combination-in-a-dict
        K,V = zip(*d.items())
        for v in product(*(v if isinstance(v,list) else traverse(v) for v in V)):
            yield dict(zip(K,v))
            
    if not os.path.isdir("images"):
        os.mkdir("images")    
    if not os.path.isdir("images\\gifs"):
        os.mkdir("images\\gifs")        
    #%% basic definitions & firts (not saved) plot  
    point = np.array([np.pi/4,np.pi/2,0,0,0,0]) 
    angles = np.linspace(0,2*np.pi,9) 
    c1 = 30 * np.cos(point[0]) + 33 * np.cos(point[0]+point[1])
    c2 = 30 * np.sin(point[0]) + 33 * np.sin(point[0]+point[1])
    radius = 18  
    targets = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T  
               
    dt = 0.001; point=point   
    shift = 0.3 ; time = 0.6 
    time_0 = time 
    time += 10*dt
    visu_pert = (targets[1:],np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1) 
    
    mech_pert1 = (np.tile(np.array([0.2,0]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    mech_pert2 = (np.tile(np.array([0,0.1]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    mech_pert3 = (np.tile(np.array([0.3,0.1]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    time = time_0
     
    params = {"init":{"eval_point":point,"time_range":time,"dt":dt}
              , "sim":{"visu_pert":[None],"mech_pert":[None,mech_pert1],"delay":[0]}}
    lqg = LQG(**params["init"])
    states, cart_states, out = lqg.simulation(targets[:-1],simul_type="LQG",visu_pert=visu_pert, delay=0,plot=False)
    plt.plot(np.linalg.norm(out.T[10:610]-targets[0,None],axis=1))
    plt.hlines(0,xmin=0,xmax=10+time/dt)
    plt.show()
    
    #%% report steps - this cell is really slow
    point = np.array([np.pi/4,np.pi/2,0,0,0,0]) 
    angles = np.linspace(0,2*np.pi,9) 
    c1 = 30 * np.cos(point[0]) + 33 * np.cos(point[0]+point[1])
    c2 = 30 * np.sin(point[0]) + 33 * np.sin(point[0]+point[1])
    radius = 18  
    targets = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T 
    
    for dt in [0.001,0.005,0.01]:
        point=point   
        shift = 0.3 ; time = 0.6 
        visu_pert = (targets[1:],np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
        time_0 = time
        time += 10*dt
        mech_pert1 = (np.tile(np.array([0.1,0]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
        time = time_0
        ddelay = np.rint(np.array([0,30,50]) / (dt/0.001)).astype('int64')  # conversion from ms to time steps  
        params = {"init":{"eval_point":[point],"time_range":[time],"dt":[dt],"w":[[100,100,10,10,10,10]],"r":[[1e-4,1e-4]]}
                      , "sim":{"simul_type":["LQG"],"visu_pert":[None,visu_pert],"mech_pert":[None,mech_pert1],"delay":ddelay.tolist()}}
             
        count = 1
        for k, v in params["sim"].items():
            count *= len(v)   
        countiter = count    
        for k, v in params["init"].items():
            countiter *= len(v) 
          
        plot=False # set it to True if you want to generate the gifs
        output = {}
        groups = {} # just serves later to verify that samples are sorted along with their init parameters first (assumption used for the boxplots)
        for ij,param in enumerate([*traverse(params)]): 
            kwargs = param['sim']
            init = param['init']
            words = [*param["init"].values(),*param["sim"].values()][1:]  
            p = [*param["init"].keys(),*param["sim"].keys()][1:]  
            filename = ""
            for pp, word in zip(p,words):
              if pp=="simul_type":
                  filename = word.lower() + "_" + filename
              elif pp=="time_range":
                  filename += "_time=" + str(word)
              elif pp=="w":
                  filename += "_w=" + str(word[::2])
              elif pp=="r": 
                  filename += "_r=" + str(word[::2]) 
              elif pp=="delay": 
                  filename += "_d=" + str(word) 
              elif pp=="dt": 
                  filename += "_dt=" + str(word) 
              else:
                  filename = filename + "_" + pp if np.any(word!=None) else filename
            lqg = LQG(**init) 
            states, cart_states, out = None, None, None
            states, cart_states, out = lqg.simulation(targets[:-1], plot=plot, filename="images\\gifs\\"+filename, **kwargs) 
            key, item = str(ij), {"args":param, "out":out, "out2":lqg.y[2:4], "out3":lqg.y[4:6], "cart_states":cart_states}
            words = [*item["args"]["init"].values(),*item["args"]["sim"].values()][1:]  
            p = [*item["args"]["init"].keys(),*item["args"]["sim"].keys()][1:] 
            filename = ""
            for pp, word in zip(p,words):
              if pp=="simul_type":
                  filename = word.lower() + "_" + filename
              elif pp=="time_range":
                  filename += "_time=" + str(word)
              elif pp=="w":
                  filename += "_w=" + str(word[::2])
              elif pp=="r": 
                  filename += "_r=" + str(word[::2]) 
              elif pp=="delay": 
                  filename += "_d=" + str(word) 
              elif pp=="dt": 
                  filename += "_dt=" + str(word) 
              else:
                  filename = filename + "_" + pp if np.any(word!=None) else filename
            out = item["out"] 
            out2,out3 = item["out2"][0],item["out2"][1] 
            out4,out5 = item["out3"][0],item["out3"][1]    
            outs = [out,out2,out3,out4,out5]
            ind = int(10 + item["args"]["init"]["time_range"] / item["args"]["init"]["dt"])
            subfig , subx = plt.subplots(3,1,gridspec_kw={'height_ratios': [2,1,1]})
            subx[1].axis("off") 
            subx[2].axis("off")
            s1 = subfig.add_subplot(4,2,5)
            s2 = subfig.add_subplot(4,2,6,sharey=s1)
            s3 = subfig.add_subplot(4,2,7,sharex=s1)
            s4 = subfig.add_subplot(4,2,8,sharex=s3,sharey=s2)
            plt.setp(s1.get_xticklabels(), visible=False)
            plt.setp(s2.get_xticklabels(), visible=False)
            plt.setp(s2.get_yticklabels(), visible=False)
            plt.setp(s4.get_yticklabels(), visible=False)
            subx[0].set_ylabel("reaching error in mm") 
            subx[0].set_xlabel("time in ms") 
            s1.set_ylabel("angular speeds") 
            s3.set_ylabel("torques") 
            s3.set_xlabel("time in ms") 
            s4.set_xlabel("time in ms") 
            subaxes = [subx[0],s1,s2,s3,s4] 
            subfig.set_size_inches([12,8]) 
            time = np.arange(0,item["args"]["init"]["time_range"]*1000,item["args"]["init"]["dt"]*1000)
            for shift in range(len(targets[:-1])):  
                i=0
                for subax, oout in zip(subaxes,outs): 
                    if i==0: 
                        subax.plot(time,np.linalg.norm(oout.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10,linewidth=0.5,label="target {}".format(shift+1)) 
                        subax.scatter(time[-1],np.linalg.norm(oout.T[ind*(shift+1)-1]-targets[shift])*10) 
                        if shift==0:
                            subax.hlines(0,xmin=0,xmax=time[-1],linewidth=0.5,linestyle="dashed",color='gray')
                    else:
                        subax.plot(time,oout.T[ind*shift+10:ind*(shift+1)],linewidth=0.5)  
                    i+=1
            subfig.legend(loc="upper right")  
            subfig.suptitle(   " ".join(filename.split("_"))   ) 
            if not os.path.isdir("images\\report"):
                os.mkdir("images\\report")
            subfig.savefig("images\\report\\"+filename+".png")
            plt.show()
            print(ij+1, "of "+str(countiter))
    print("Finished (:")        
    plt.close('all')
            
    #%% W weights
    point = np.array([np.pi/4,np.pi/2,0,0,0,0]) 
    angles = np.linspace(0,2*np.pi,9) 
    c1 = 30 * np.cos(point[0]) + 33 * np.cos(point[0]+point[1])
    c2 = 30 * np.sin(point[0]) + 33 * np.sin(point[0]+point[1])
    radius = 18  
    targets = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T 
    time = 0.6
    
    ww = []
    for w1 in [25,50,100,200,500]:
        for w2 in [0,10,100]:
            for w3 in [0,10,100]:
                ww.append([w1,w1,w2,w2,w3,w3])
             
    for dt in [0.01,0.005,0.001]:
        for w in ww: 
            params = {"init":{"eval_point":[point],"time_range":[time],"dt":[dt],"w":[w],"r":[[1e-1,1e-1],[1e-2,1e-2],[1e-4,1e-4],[1,1]]}
                      , "sim":{"simul_type":["LQG"],"visu_pert":[None],"mech_pert":[None],"delay":[0]}}
             
            count = 1
            for k, v in params["sim"].items():
                count *= len(v)
            countiter = count    
            for k, v in params["init"].items():
                countiter *= len(v)    
             
            plot=False
            output = {}
            groups = {} # just serves later to verify that samples are sorted along with their init parameters first (assumption used for the boxplots)
            for ij,param in enumerate([*traverse(params)]): 
                kwargs = param['sim']
                init = param['init']
                words = [*param["init"].values(),*param["sim"].values()][1:]  
                p = [*param["init"].keys(),*param["sim"].keys()][1:]  
                filename = "lqg"
                for pp, word in zip(p,words):
                  if pp=="time_range":
                      filename += "_time=" + str(word)
                  elif pp.split("_")[-1]=="pert":
                      filename = filename + "_" + pp  if np.any(word!=None) else filename  
                  elif word!=True and word!=False:
                      filename += "_" + pp + "=" + str(word) 
                  elif pp=="delay" and word!=0:
                      filename += "_" + pp + "=" + str(word) 
                  else:
                      filename = filename + "_" + pp  if word else filename
                lqg = LQG(**init) 
                states, cart_states, out = None, None, None
                states, cart_states, out = lqg.simulation(targets[:-1], plot=plot, filename="images\\gifs\\"+filename+'.gif', **kwargs)
                output[str(ij)] = {"args":param, "out":out, "cart_states":cart_states, "states":states} #, "obj":lqg}
                if groups == {}:
                    groups[str(1)] = {"args":param["init"], "list":[ij]}
                else: 
                    found_a_group = False
                    for group_name, group in groups.items():
                        group_args = group["args"] 
                        for (_, val1), (_, val2) in zip(group_args.items(),param["init"].items()):
                            if np.any(val1!=val2): # diff groups => diff params
                                break
                        else:
                            found_a_group = True
                            group['list'] = np.append(group['list'],ij)
                        if found_a_group:
                            break
                    else:
                        groups[str( int([*groups.keys()][-1]) +1)] = {"args":param["init"], "list":[ij]} 
                print(ij+1, "of "+str(countiter))  
            print("Finished ;)")  
            
            fig, (ax0, ax1, ax2) = plt.subplots(3,1,gridspec_kw={'height_ratios': [1,1,1]})
            fig.set_size_inches([24,15]) 
            ax1.set_title("label convention: sum of 100 if visual pert, 10 if mechanical pert, 0 else and delay in seconds") 
            ax1.set_ylabel("reaching error in mm") 
            labels = []
            box = [] 
            box2 = [] 
            words = [*output["0"]["args"]["init"].values(),*output["0"]["args"]["sim"].values()][1:]  
            p = [*output["0"]["args"]["init"].keys(),*output["0"]["args"]["sim"].keys()][1:] 
            filename = ""
            folder = "stab"
            for pp, word in zip(p,words):
              if pp=="simul_type":
                  filename = word.lower() + "_" + filename
              elif pp=="time_range":
                  filename += "_time=" + str(word)
              elif pp=="w":
                  filename += "_w=" + str(word[::2])
              elif pp=="r":
                  filename += "_r=" + str(word[::2]) 
                  folder += "_r=" + str(word[::2]) 
              elif pp=="dt":
                  filename += "_dt=" + str(word) 
                  folder += "_dt=" + str(word) 
              else:
                  filename = filename  
            fig.suptitle(  " ".join(filename.split("_"))  )
            for key, item in output.items(): 
                outs = []
                outs2 = []
                out = item["out"] 
                ind = int(10 + item["args"]["init"]["time_range"] / item["args"]["init"]["dt"])
                subfig, subax = plt.subplots(1,1) 
                subfig.set_size_inches([15,8])
                subax.hlines(0,xmin=0,xmax=ind-10+ind//8)   
                for shift in range(len(targets[:-1])): 
                    outs.append([np.linalg.norm(out[:,ind*(shift+1)-1]-targets[shift])*10])
                    outs2.append([(np.linalg.norm(out.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10).mean()])
                    if (int(key)+1)%count==0:
                      ax0.plot(np.linalg.norm(out.T[ind*shift+10:ind*(shift+1)+ind//8]-targets[shift,None],axis=1)*10,linewidth=0.5) 
                      ax0.scatter(ind-11,np.linalg.norm(out.T[ind*(shift+1)-1]-targets[shift])*10)  
                    subax.plot(np.linalg.norm(out.T[ind*shift+10:ind*(shift+1)+ind//8]-targets[shift,None],axis=1)*10,linewidth=0.5) 
                    subax.scatter(ind-11,np.linalg.norm(out.T[ind*(shift+1)-1]-targets[shift])*10)
                if not os.path.isdir("images\\"+folder):
                    os.mkdir("images\\"+folder)
                subfig.savefig("images\\"+folder+"\\"+filename+".png")
                box.append(np.array(outs).flatten())
                box2.append(np.array(outs).flatten())
                label = 0
                items = np.sort([*item['args']["sim"].keys()])
                for nval, ikey in enumerate(items):
                    val = item['args']["sim"][ikey]
                    label += (10**nval)*1 if np.any(val!=None) and ikey!="delay" else (10**nval)*val*item["args"]["init"]["dt"] if ikey=="delay" else 0  
                labels.append(str(np.around(label,decimals=4)))  
                if (int(key)+1)%(20)==0 or (int(key)+1)%count==0: 
                    axes = [ax1,ax2] # [ax2,ax4] if (int(key)+1)%count==0 else [ax1,ax3]
                    boxes = [box,box2]
                    for ax, box in zip(axes,boxes):
                        ax.set_title("label convention: sum of 100 if visual pert, 10 if mechanical pert, 0 else and delay in seconds") 
                        ax.set_ylabel("reaching error in mm") 
                        ax.boxplot(box)  
                        ax.set_xticklabels(labels)
                    labels = []
                    box = []   
                    box2 = []   
                if (int(key)+1)%count==0:
                    ax0.hlines(0,xmin=0,xmax=ind+ind//8)    
                    words = [*item["args"]["init"].values(),*item["args"]["sim"].values()][1:]  
                    p = [*item["args"]["init"].keys(),*item["args"]["sim"].keys()][1:]  
                    filename = ""
                    folder = "stab"  
                    for pp, word in zip(p,words):
                      if pp=="simul_type":
                          filename = word.lower() + "_" + filename
                      elif pp=="time_range":
                          filename += "_time=" + str(word)
                      elif pp=="w":
                          filename += "_w=" + str(word[::2])
                      elif pp=="r":
                          filename += "_r=" + str(word[::2]) 
                          folder += "_r=" + str(word[::2]) 
                      elif pp=="dt":
                          filename += "_dt=" + str(word) 
                          folder += "_dt=" + str(word) 
                      # elif pp=="delay":
                      #     filename = filename + "_delay=" + str(word) if word!=0 else filename
                      # elif pp in ("visu_pert", "mech_pert"):
                      #     filename = filename + pp.split('_')[0] + "=" + str(np.round(word[0][0],3)) if np.any(word!=None) else filename
                      else:
                          filename = filename
                    fig.suptitle(   " ".join(filename.split("_"))   )
                  #  fig.savefig("images\\"+filename+".png") 
                    plt.show()
                    plt.close('all')
                    if key!=[*output.keys()][-1]: 
                        fig, (ax0, ax1, ax2) = plt.subplots(3,1,gridspec_kw={'height_ratios': [1,1,1]})
                        fig.set_size_inches([24,15])  
                    labels = []
                    box = []  
                    box2 = []  
    #%% h with non-diagonal terms
    point = np.array([np.pi/4,np.pi/2,0,0,0,0]) 
    angles = np.linspace(0,2*np.pi,9) 
    c1 = 30 * np.cos(point[0]) + 33 * np.cos(point[0]+point[1])
    c2 = 30 * np.sin(point[0]) + 33 * np.sin(point[0]+point[1])
    radius = 18 
    # targets = np.tile([c1,c2],(2*angles.size,1))
    # targets[np.arange(2*angles.size)%2!=0] = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T 
    targets = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T 
    time = 0.6
    
    targs = [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)] 
    hps = []
    for mean in [0,1e-4,1e-3,1e-2]:
        for std in [1e-4,1e-3,1e-2]:
            for targ_or_not in targs:
                for is_diag in [True,False]:
                    hps.append([mean,std,targ_or_not,is_diag])
    dt = 0.001   
    params = {"init":{"eval_point":[point],"time_range":[time],"dt":[dt],"w":[[100,100,10,10,10,10]],"r":[[1e-4,1e-4]]}
              , "sim":{"simul_type":["LQG"],"visu_pert":[None],"mech_pert":[None],"delay":[0],"hpar":hps}}
     
    count = 1
    for k, v in params["sim"].items():
        count *= len(v)    
    countiter = count    
    for k, v in params["init"].items():
        countiter *= len(v)
        
    plot=False
    output = {}
    groups = {} # just serves later to verify that samples are sorted along with their init parameters first (assumption used for the boxplots)
    for ij,param in enumerate([*traverse(params)]): 
        kwargs = param['sim']
        init = param['init']
        words = [*param["init"].values(),*param["sim"].values()][1:]  
        p = [*param["init"].keys(),*param["sim"].keys()][1:]  
        filename = "lqg"
        for pp, word in zip(p,words):
          if pp=="time_range":
              filename += "_time=" + str(word)
          elif pp.split("_")[-1]=="pert":
              filename = filename + "_" + pp  if np.any(word!=None) else filename  
          elif word!=True and word!=False:
              filename += "_" + pp + "=" + str(word) 
          elif pp=="delay" and word!=0:
              filename += "_" + pp + "=" + str(word) 
          else:
              filename = filename + "_" + pp  if word else filename
        lqg = LQG(**init) 
        states, cart_states, out = None, None, None
        states, cart_states, out = lqg.simulation(targets[:-1], plot=plot, filename="images\\gifs\\"+filename+'.gif', **kwargs)
        output[str(ij)] = {"args":param, "out":out, "cart_states":cart_states, "states":states} #, "obj":lqg}
        if groups == {}:
            groups[str(1)] = {"args":param["init"], "list":[ij]}
        else: 
            found_a_group = False
            for group_name, group in groups.items():
                group_args = group["args"] 
                for (_, val1), (_, val2) in zip(group_args.items(),param["init"].items()):
                    if np.any(val1!=val2): # diff groups => diff params
                        break
                else:
                    found_a_group = True
                    group['list'] = np.append(group['list'],ij)
                if found_a_group:
                    break
            else:
                groups[str( int([*groups.keys()][-1]) +1)] = {"args":param["init"], "list":[ij]} 
        print(ij+1, "of "+str(countiter))   
    print("Finished ;)")  
    
    boxmean = {"mean":{},"var":{},"t":{},"d":{}}
    for m,v,t,d in hps:
        boxmean["mean"][str(m)] = []
        boxmean["var"][str(v)] = []   
        boxmean["t"][str(t)] = []   
        boxmean["d"][str(d)] = [] 
              
    box = [] 
    boxall = [] 
    time = None
    for key, item in output.items():   
        words = [*item["args"]["init"].values(),*item["args"]["sim"].values()][1:]  
        p = [*item["args"]["init"].keys(),*item["args"]["sim"].keys()][1:] 
        filename = ""
        for pp, word in zip(p,words):
          if pp=="simul_type":
              filename = word.lower() + "_" + filename
          elif pp=="time_range":
              filename += "_time=" + str(word)
          elif pp=="w":
              filename += "_w=" + str(word[::2])
          elif pp=="r": 
              filename += "_r=" + str(word[::2]) 
          elif pp=="delay": 
              filename += "_d=" + str(word) 
          elif pp=="dt": 
              filename += "_dt=" + str(word) 
          elif pp=="hpar": 
              filename += "_hpar=" + str(word) 
          else:
              filename = filename + "_" + pp + "=" + str(word) if np.any(word!=None) else filename
        out = item["out"]      
        ind = int(10 + item["args"]["init"]["time_range"] / item["args"]["init"]["dt"])
        subfig , subx = plt.subplots(1,1) 
        subfig.set_size_inches([12,8]) 
        time = np.arange(0,item["args"]["init"]["time_range"]*1000,item["args"]["init"]["dt"]*1000)
        for shift in range(len(targets[:-1])):  
            i=0
            box.append(np.linalg.norm(out.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10)
            subx.plot(time,np.linalg.norm(out.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10,linewidth=0.5,label="target {}".format(shift+1))
            
        boxall.append(box)
        for m in ["0","0.0001","0.001","0.01"]:
            if m==str(item["args"]["sim"]["hpar"][0]):
                boxmean["mean"][str(m)].append(box)   
                
        for v in ["0.0001","0.001","0.01"]:
            if v==str(item["args"]["sim"]["hpar"][1]):
                boxmean["var"][str(v)].append(box)   
        for t in [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)] :        
            if len(set({t,item["args"]["sim"]["hpar"][2]}))==len(set({t})):
                boxmean["t"][str(t)].append(box)   
        for d in ["True","False"]:        
            if d==str(item["args"]["sim"]["hpar"][3]):
                boxmean["d"][str(d)].append(box)  
                
        box = []      
        subfig.legend(loc="upper right")  
        subfig.suptitle(   " ".join(filename.split("_"))   ) 
        if not os.path.isdir("images\\reportH"):
            os.mkdir("images\\reportH")
        subfig.savefig("images\\reportH\\"+filename+".png")
        plt.show()
        plt.close('all') 
    
    for condition in ["mean","var","t","d"]:
        moyennes = [] # mean accross conditions 
        for key in boxmean[condition].keys():
            moyennes.append(boxmean[condition][key])
        
        for k,moy in enumerate(moyennes):  
            mean_m = np.mean(moy,axis=0) # mean accross samples of similar conditions
            variances = np.std(moy,axis=0) 
            for ia, moyenne, var in zip(np.arange(8),mean_m,variances):
                fig, ax = plt.subplots(1,1)
                fig.set_size_inches([16,8])    
                ax.fill_between(x=time,y1=moyenne+var,y2=moyenne-var,alpha=0.5,color="lightgreen")   
                ax.plot(time,moyenne)  
                ax.set_title("target {}".format(ia+1))  
                if not os.path.isdir("images\\reportH_{}".format(condition)):
                    os.mkdir("images\\reportH_{}".format(condition))
                fig.savefig("images\\reportH_{}\\{}_target{}.png".format(condition,[*boxmean[condition].keys()][k],ia)) 
                plt.show()
        plt.close('all') 
        
    for condition in ["mean","var","t","d"]:
        moyennes = [] # mean accross conditions 
        for key in boxmean[condition].keys():
            moyennes.append(boxmean[condition][key])
        
        for k,moy in enumerate(moyennes):  
            moyenne = np.mean(moy,axis=(0,1)) # mean accross samples of similar conditions
            variances = np.std(moy,axis=(0,1))  
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches([16,8])    
            ax.fill_between(x=time,y1=moyenne+var,y2=moyenne-var,alpha=0.5,color="lightgreen")   
            ax.plot(time,moyenne)  
            ax.set_title("Mean accross targets")  
            fig.savefig("images\\reportH_{}\\{}_all.png".format(condition,[*boxmean[condition].keys()][k],ia)) 
            plt.show()
        plt.close('all') 
    
    moyennes = np.mean(boxall,axis=1) # mean accross targets
    variances = np.std(boxall,axis=1) 
    for k, moyenne, var in zip(np.arange(count),moyennes,variances):
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches([16,8])    
        ax.fill_between(x=time,y1=moyenne+var,y2=moyenne-var,alpha=0.5,color="lightgreen")   
        ax.plot(time,moyenne)  
        cond = output[str(k)]["args"]["sim"]["hpar"]
        ax.set_title("Condition {}".format(cond))  
        if not os.path.isdir("images\\reportH_cond"):
            os.mkdir("images\\reportH_cond")
        fig.savefig("images\\reportH_cond\\{}.png".format(cond)) 
        plt.show()
    plt.close('all') 
    
    moyennes = np.mean(boxall,axis=0) # mean accross all conditions
    variances = np.std(boxall,axis=0) 
    for ia, moyenne, var in zip(np.arange(8),moyennes,variances):
        fig, ax = plt.subplots(1,1)
        fig.set_size_inches([16,8])    
        ax.fill_between(x=time,y1=moyenne+var,y2=moyenne-var,alpha=0.5,color="lightgreen")   
        ax.plot(time,moyenne)   
        ax.set_title("Target {}".format(ia))  
        if not os.path.isdir("images\\reportH_targ"):
            os.mkdir("images\\reportH_targ")
        fig.savefig("images\\reportH_targ\\{}.png".format(ia)) 
        plt.show()
    plt.close('all') 
    #%% h with non-diagonal terms
    point = np.array([np.pi/4,np.pi/2,0,0,0,0]) 
    angles = np.linspace(0,2*np.pi,9) 
    c1 = 30 * np.cos(point[0]) + 33 * np.cos(point[0]+point[1])
    c2 = 30 * np.sin(point[0]) + 33 * np.sin(point[0]+point[1])
    radius = 18  
    targets = np.array([c1 + radius*np.cos(angles), c2 + radius*np.sin(angles)]).T 
    time = 0.6
    
    targs = [(0,1,2),(1,),(0,1),(0,2)]  
    hps = [None]
    means = [0,1e-2,2*1e-2]
    stds = [1e-1,1e-2,1e-3]
    for mean in means:
        for std in stds:
            for targ_or_not in targs:
                for is_diag in [True,False]:
                    hps.append([mean,std,targ_or_not,is_diag])
    dt = 0.001   
    shift = 0.3 ; time = 0.6 
    time_0 = time
    time += 10*dt
    mech_pert0 = (np.tile(np.array([0.1,0]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    mech_pert1 = (np.tile(np.array([0.2,0]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    mech_pert2 = (np.tile(np.array([0.3,0]),(targets.shape[0]-1,1)),np.arange(0,time*targets.shape[0],time)+shift,np.ones(targets.shape[0])*0.1)
    time = time_0
    params = {"init":{"eval_point":[point],"time_range":[time],"dt":[dt],"w":[[100,100,10,10,10,10]],"r":[[1e-4,1e-4]]}
              , "sim":{"simul_type":["LQG"],"visu_pert":[None],"mech_pert":[mech_pert0,mech_pert1,mech_pert2],"delay":[0],"hpar":hps}}
     
    count = 1
    for k, v in params["sim"].items():
        count *= len(v) 
    countiter = count    
    for k, v in params["init"].items():
        countiter *= len(v)    
     
    plot=False
    output = {}
    groups = {} # just serves later to verify that samples are sorted along with their init parameters first (assumption used for the boxplots)
    for ij,param in enumerate([*traverse(params)]): 
        kwargs = param['sim']
        init = param['init']
        words = [*param["init"].values(),*param["sim"].values()][1:]  
        p = [*param["init"].keys(),*param["sim"].keys()][1:]  
        filename = "lqg"
        for pp, word in zip(p,words):
          if pp=="time_range":
              filename += "_time=" + str(word)
          elif pp.split("_")[-1]=="pert":
              filename = filename + "_" + pp  if np.any(word!=None) else filename  
          elif word!=True and word!=False:
              filename += "_" + pp + "=" + str(word) 
          elif pp=="delay" and word!=0:
              filename += "_" + pp + "=" + str(word) 
          else:
              filename = filename + "_" + pp  if word else filename
        lqg = LQG(**init) 
        states, cart_states, out = None, None, None
        states, cart_states, out = lqg.simulation(targets[:-1], plot=plot, filename="images\\gifs\\"+filename+'.gif', **kwargs)
        output[str(ij)] = {"args":param, "out":out, "cart_states":cart_states, "states":states} #, "obj":lqg}
        if groups == {}:
            groups[str(1)] = {"args":param["init"], "list":[ij]}
        else: 
            found_a_group = False
            for group_name, group in groups.items():
                group_args = group["args"] 
                for (_, val1), (_, val2) in zip(group_args.items(),param["init"].items()):
                    if np.any(val1!=val2): # diff groups => diff params
                        break
                else:
                    found_a_group = True
                    group['list'] = np.append(group['list'],ij)
                if found_a_group:
                    break
            else:
                groups[str( int([*groups.keys()][-1]) +1)] = {"args":param["init"], "list":[ij]} 
        print(ij+1, "of "+str(countiter))  
    print("Finished ;)")   
    
    boxmean ={"None":{},"hpar":{"mean":{},"var":{},"t":{},"d":{}}}
    for anhp in hps:
        if np.any(anhp!=None):
            m,v,t,d = anhp 
            boxmean["hpar"]["mean"][str(m)] = {str(mech_pert0[0][0]):[],str(mech_pert1[0][0]):[],str(mech_pert2[0][0]):[]} 
            boxmean["hpar"]["var"][str(v)] = {str(mech_pert0[0][0]):[],str(mech_pert1[0][0]):[],str(mech_pert2[0][0]):[]} 
            boxmean["hpar"]["t"][str(t)] = {str(mech_pert0[0][0]):[],str(mech_pert1[0][0]):[],str(mech_pert2[0][0]):[]} 
            boxmean["hpar"]["d"][str(d)] = {str(mech_pert0[0][0]):[],str(mech_pert1[0][0]):[],str(mech_pert2[0][0]):[]} 
       
    boxmean["None"]= {str(mech_pert0[0][0]):[],str(mech_pert1[0][0]):[],str(mech_pert2[0][0]):[]} 
              
    box = [] 
    boxall = [] 
    time = None
    for key, item in output.items():   
        words = [*item["args"]["init"].values(),*item["args"]["sim"].values()][1:]  
        p = [*item["args"]["init"].keys(),*item["args"]["sim"].keys()][1:] 
        filename = ""
        for pp, word in zip(p,words):
          if pp=="simul_type":
              filename = word.lower() + "_" + filename
          elif pp=="time_range":
              filename += "_time=" + str(word)
          elif pp=="w":
              filename += "_w=" + str(word[::2])
          elif pp=="r": 
              filename += "_r=" + str(word[::2]) 
          elif pp=="delay": 
              filename += "_d=" + str(word) 
          elif pp=="dt": 
              filename += "_dt=" + str(word) 
          elif pp=="hpar": 
              filename += "_hpar=" + str(word)
          elif pp=="mech_pert": 
              filename += "_pert=" + str(word[0][0]) 
          else:
              filename = filename  
        outs = item["out"]      
        ind = int(10 + item["args"]["init"]["time_range"] / item["args"]["init"]["dt"])
        subfig , subx = plt.subplots(1,1) 
        subfig.set_size_inches([12,8]) 
        time = np.arange(0,item["args"]["init"]["time_range"]*1000,item["args"]["init"]["dt"]*1000)
        for shift in range(len(targets[:-1])):  
            i=0
            box.append(np.linalg.norm(outs.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10)
            subx.plot(time,np.linalg.norm(outs.T[ind*shift+10:ind*(shift+1)]-targets[shift,None],axis=1)*10,linewidth=0.5,label="target {}".format(shift+1))
            
        boxall.append(box)
        if np.any(item["args"]["sim"]["hpar"]!=None):
            for m in [str(mmean) for mmean in means]:
                if m==str(item["args"]["sim"]["hpar"][0]):
                    for kk in boxmean["hpar"]["mean"][str(m)].keys():
                        if np.all(kk == str(item["args"]["sim"]["mech_pert"][0][0])):
                            boxmean["hpar"]["mean"][str(m)][kk].append(box) 
                    
            for v in [str(sstd) for sstd in stds]:
                if v==str(item["args"]["sim"]["hpar"][1]):
                    for kk in boxmean["hpar"]["var"][str(v)].keys():
                        if np.all(kk == str(item["args"]["sim"]["mech_pert"][0][0])):
                            boxmean["hpar"]["var"][str(v)][kk].append(box)
                            
            for t in [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)] :        
                if len(set({t,item["args"]["sim"]["hpar"][2]}))==len(set({t})):
                    for kk in boxmean["hpar"]["t"][str(t)].keys():
                        if np.all(kk == str(item["args"]["sim"]["mech_pert"][0][0])):
                            boxmean["hpar"]["t"][str(t)][kk].append(box) 
                            
            for d in ["True","False"]:        
                if d==str(item["args"]["sim"]["hpar"][3]):
                    for kk in boxmean["hpar"]["d"][str(d)].keys():
                        if np.all(kk == str(item["args"]["sim"]["mech_pert"][0][0])):
                            boxmean["hpar"]["d"][str(d)][kk].append(box) 
        else:            
            for kk in boxmean["None"].keys():
                if np.all(kk == str(item["args"]["sim"]["mech_pert"][0][0])):
                    boxmean["None"][kk].append(box)                   
                
        box = []      
        subfig.legend(loc="upper right")  
        subfig.suptitle(   " ".join(filename.split("_"))   ) 
        if not os.path.isdir("images\\reportH_pert"):
            os.mkdir("images\\reportH_pert")
        subfig.savefig("images\\reportH_pert\\"+filename+".png")
        plt.show()
        plt.close('all') 
    
    for condition in ["mean","var","t","d"]:
        moyennes1 = []   
        moyennes2 = []   
        moyennes0 = []  
        controle1 = []   
        controle2 = []   
        controle0 = []  
        for key in boxmean["hpar"][condition].keys():
            moyennes0.append(boxmean["hpar"][condition][key][str(mech_pert0[0][0])])
            moyennes1.append(boxmean["hpar"][condition][key][str(mech_pert1[0][0])])
            moyennes2.append(boxmean["hpar"][condition][key][str(mech_pert2[0][0])])
            controle0.append(boxmean["None"][str(mech_pert0[0][0])])
            controle1.append(boxmean["None"][str(mech_pert1[0][0])])
            controle2.append(boxmean["None"][str(mech_pert2[0][0])]) 
        
        for k,moy,moy1,moy2,c0,c1,c2 in zip(np.arange(len(moyennes1)),moyennes0,moyennes1,moyennes2,controle0,controle1,controle2):  
            mean_m = np.mean(moy,axis=0) # mean accross samples of similar conditions
            variances = np.std(moy,axis=0) 
            mean_m1 = np.mean(moy1,axis=0) # mean accross samples of similar conditions
            variances1 = np.std(moy1,axis=0) 
            mean_m2 = np.mean(moy2,axis=0) # mean accross samples of similar conditions
            variances2 = np.std(moy2,axis=0) 
            mean_c = np.mean(c0,axis=0) # mean accross samples of similar conditions
            vc = np.std(c0,axis=0) 
            mean_c1 = np.mean(c1,axis=0) # mean accross samples of similar conditions
            vc1 = np.std(c1,axis=0) 
            mean_c2 = np.mean(c2,axis=0) # mean accross samples of similar conditions
            vc2 = np.std(c2,axis=0) 
            for ia,moyenne,var,moyenne1,var1,moyenne2,var2,mc1,v1,mc2,v2,mc3,v3 in zip(np.arange(8),mean_m,variances,mean_m1,variances1,mean_m2,variances2,mean_c,vc,mean_c1,vc1,mean_c2,vc2):
                fig, ax = plt.subplots(1,1)
                fig.set_size_inches([16,8])       
                ax.plot(time,moyenne,color="red",label="pert 1")      
                ax.plot(time,moyenne1,color="red",label="pert 2")     
                ax.plot(time,moyenne2,color="red",label="pert 3")       
                ax.plot(time,mc1,color="green",label="controle with pert 1")  
                ax.plot(time,mc2,color="green",label="controle with pert 2")       
                ax.plot(time,mc3,color="green",label="controle with pert 3")  
                ax.set_title("target {}".format(ia+1))  
                if not os.path.isdir("images\\reportH_pert_{}".format(condition)):
                    os.mkdir("images\\reportH_pert_{}".format(condition))
                fig.legend()    
                fig.savefig("images\\reportH_pert_{}\\{}_target{}.png".format(condition,[*boxmean["hpar"][condition].keys()][k],ia)) 
                plt.show() 
        plt.close('all') 
        
    for condition in ["mean","var","t","d"]:
        moyennes1 = []   
        moyennes2 = []   
        moyennes0 = []  
        controle1 = []   
        controle2 = []   
        controle0 = []  
        for key in boxmean["hpar"][condition].keys():
            moyennes0.append(boxmean["hpar"][condition][key][str(mech_pert0[0][0])])
            moyennes1.append(boxmean["hpar"][condition][key][str(mech_pert1[0][0])])
            moyennes2.append(boxmean["hpar"][condition][key][str(mech_pert2[0][0])])
            controle0.append(boxmean["None"][str(mech_pert0[0][0])])
            controle1.append(boxmean["None"][str(mech_pert1[0][0])])
            controle2.append(boxmean["None"][str(mech_pert2[0][0])]) 
        
        for k,moy,moy1,moy2,c0,c1,c2 in zip(np.arange(len(moyennes1)),moyennes0,moyennes1,moyennes2,controle0,controle1,controle2):  
            mean_m = np.mean(moy,axis=(0,1)) # mean accross samples of similar conditions
            variances = np.std(moy,axis=(0,1)) 
            mean_m1 = np.mean(moy1,axis=(0,1)) # mean accross samples of similar conditions
            variances1 = np.std(moy1,axis=(0,1)) 
            mean_m2 = np.mean(moy2,axis=(0,1)) # mean accross samples of similar conditions
            variances2 = np.std(moy2,axis=(0,1)) 
            mean_c = np.mean(c0,axis=(0,1)) # mean accross samples of similar conditions
            vc = np.std(c0,axis=(0,1)) 
            mean_c1 = np.mean(c1,axis=(0,1)) # mean accross samples of similar conditions
            vc1 = np.std(c1,axis=(0,1)) 
            mean_c2 = np.mean(c2,axis=(0,1)) # mean accross samples of similar conditions
            vc2 = np.std(c2,axis=(0,1)) 
            moyenne,var,moyenne1,var1,moyenne2,var2,mc1,v1,mc2,v2,mc3,v3 = mean_m,variances,mean_m1,variances1,mean_m2,variances2,mean_c,vc,mean_c1,vc1,mean_c2,vc2 
            fig, ax = plt.subplots(1,1)
            fig.set_size_inches([16,8])       
            ax.plot(time,moyenne,color="red",label="pert 1")      
            ax.plot(time,moyenne1,color="red",label="pert 2")     
            ax.plot(time,moyenne2,color="red",label="pert 3")       
            ax.plot(time,mc1,color="green",label="controle with pert 1")  
            ax.plot(time,mc2,color="green",label="controle with pert 2")       
            ax.plot(time,mc3,color="green",label="controle with pert 3")  
            ax.set_title("Condition {} = {}".format(condition,[*boxmean["hpar"][condition].keys()][k]))  
            if not os.path.isdir("images\\reportH_pert_{}".format(condition)):
                os.mkdir("images\\reportH_pert_{}".format(condition))
            fig.legend()    
            fig.savefig("images\\reportH_pert_{}\\{}_all.png".format(condition,[*boxmean["hpar"][condition].keys()][k])) 
            plt.show() 
        plt.close('all') 
         
    for mymy in [str(mmean) for mmean in means]:
        condition0 = "mean = " + mymy
        for keyn, key in enumerate([str(sstd) for sstd in stds]): 
            moyennes1 = []   
            moyennes2 = []   
            moyennes0 = []  
            controle1 = []   
            controle2 = []   
            controle0 = [] 
            condition = condition0 + ", std = {}".format(key)
            
            seg = boxmean["hpar"]["mean"][mymy]
            pert0, pert1, pert2 = tuple(boxmean["hpar"]["mean"][mymy].keys()) 
            moyennes0.append(seg[pert0][::(len(targs)+2)])
            moyennes1.append(seg[pert1][::(len(targs)+2)])
            moyennes2.append(seg[pert2][::(len(targs)+2)])
            controle0.append(boxmean["None"][pert0][::(len(targs)+2)])
            controle1.append(boxmean["None"][pert1][::(len(targs)+2)])
            controle2.append(boxmean["None"][pert2][::(len(targs)+2)])  
        
            for k,moy,moy1,moy2,c0,c1,c2 in zip(np.arange(len(moyennes1)),moyennes0,moyennes1,moyennes2,controle0,controle1,controle2):  
                mean_m = np.mean(moy,axis=(0,1)) # mean accross samples of similar conditions
                variances = np.std(moy,axis=(0,1)) 
                mean_m1 = np.mean(moy1,axis=(0,1)) # mean accross samples of similar conditions
                variances1 = np.std(moy1,axis=(0,1)) 
                mean_m2 = np.mean(moy2,axis=(0,1)) # mean accross samples of similar conditions
                variances2 = np.std(moy2,axis=(0,1)) 
                mean_c = np.mean(c0,axis=(0,1)) # mean accross samples of similar conditions
                vc = np.std(c0,axis=(0,1)) 
                mean_c1 = np.mean(c1,axis=(0,1)) # mean accross samples of similar conditions
                vc1 = np.std(c1,axis=(0,1)) 
                mean_c2 = np.mean(c2,axis=(0,1)) # mean accross samples of similar conditions
                vc2 = np.std(c2,axis=(0,1))  
                moyenne,var,moyenne1,var1,moyenne2,var2,mc1,v1,mc2,v2,mc3,v3 = mean_m,variances,mean_m1,variances1,mean_m2,variances2,mean_c,vc,mean_c1,vc1,mean_c2,vc2 
                fig, ax = plt.subplots(1,1)
                fig.set_size_inches([16,8])       
                ax.plot(time,moyenne,color="red",label="pert 1")      
                ax.plot(time,moyenne1,color="red",label="pert 2")     
                ax.plot(time,moyenne2,color="red",label="pert 3")       
                ax.plot(time,mc1,color="green",label="controle with pert 1")  
                ax.plot(time,mc2,color="green",label="controle with pert 2")       
                ax.plot(time,mc3,color="green",label="controle with pert 3")  
                ax.set_title("Condition {}".format(condition))  
                if not os.path.isdir("images\\reportH_pert_{}".format( "".join(tuple(condition0.split(" "))) ) ):
                    os.mkdir("images\\reportH_pert_{}".format( "".join(tuple(condition0.split(" "))) ))
                fig.legend()    
                fig.savefig("images\\reportH_pert_{}\\{}_all.png".format("".join(tuple(condition0.split(" "))),key)) 
                plt.show() 
        plt.close('all')   
         
        for keyn, key in enumerate([str(sstd) for sstd in stds]): 
            moyennes1 = []   
            moyennes2 = []   
            moyennes0 = []  
            controle1 = []   
            controle2 = []   
            controle0 = [] 
            condition = condition0 + ", std = {}".format(key)
            
            seg = boxmean["hpar"]["mean"][mymy]
            pert0, pert1, pert2 = tuple(boxmean["hpar"]["mean"][mymy].keys()) 
            moyennes0.append(seg[pert0][::(len(targs)+2)])
            moyennes1.append(seg[pert1][::(len(targs)+2)])
            moyennes2.append(seg[pert2][::(len(targs)+2)])
            controle0.append(boxmean["None"][pert0][::(len(targs)+2)])
            controle1.append(boxmean["None"][pert1][::(len(targs)+2)])
            controle2.append(boxmean["None"][pert2][::(len(targs)+2)])  
        
            for k,moy,moy1,moy2,c0,c1,c2 in zip(np.arange(len(moyennes1)),moyennes0,moyennes1,moyennes2,controle0,controle1,controle2):  
                mean_m = np.mean(moy,axis=0)  
                variances = np.std(moy,axis=0) 
                mean_m1 = np.mean(moy1,axis=0)  
                variances1 = np.std(moy1,axis=0) 
                mean_m2 = np.mean(moy2,axis=0)  
                variances2 = np.std(moy2,axis=0) 
                mean_c = np.mean(c0,axis=0)  
                vc = np.std(c0,axis=0) 
                mean_c1 = np.mean(c1,axis=0)  
                vc1 = np.std(c1,axis=0) 
                mean_c2 = np.mean(c2,axis=0)  
                vc2 = np.std(c2,axis=0) 
                for ia,moyenne,var,moyenne1,var1,moyenne2,var2,mc1,v1,mc2,v2,mc3,v3 in zip(np.arange(8),mean_m,variances,mean_m1,variances1,mean_m2,variances2,mean_c,vc,mean_c1,vc1,mean_c2,vc2):
                    fig, ax = plt.subplots(1,1)
                    fig.set_size_inches([16,8])       
                    ax.plot(time,moyenne,label="pert 1",linewidth=0.7)      
                    ax.plot(time,moyenne1,label="pert 2",linewidth=0.7)     
                    ax.plot(time,moyenne2,label="pert 3",linewidth=0.7)       
                    ax.plot(time,mc1,color="gray",label="control with pert 1")  
                    ax.plot(time,mc2,color="gray",label="control with pert 2")       
                    ax.plot(time,mc3,color="gray",label="control with pert 3")  
                    ax.set_title("Target {}".format(ia)) 
                    if not os.path.isdir("images\\reportH_pert_{}".format("".join(tuple(condition0.split(" "))))):
                        os.mkdir("images\\reportH_pert_{}".format("".join(tuple(condition0.split(" ")))))
                    fig.legend()    
                    fig.savefig("images\\reportH_pert_{}\\{}_target{}.png".format("".join(tuple(condition0.split(" "))),key,ia))
                    plt.show()
    plt.close('all')       
     
    #%% to empty nearly empty folders into images folder ==> play with len(files)<=x to tune it (: 
    import os 
    import glob
    import shutil
    
    root = os.path.join(os.getcwd(),"images\\")
    with os.scandir(root) as it:
        
        for entry in it:
            path = os.path.join(root,entry.name)
            if os.path.isdir(path):
                files = os.listdir(path)   
                if len(files)<=-1:
                    shutil.move(os.path.join(path, files[0]) , root + files[0] )
                    shutil.rmtree(path) 
                if  entry.name.startswith("reportH") and entry.name!="reportH":
                    shutil.move(path, os.path.join("images\\reportH",entry.name))
                elif  entry.name.startswith("stab") and entry.name!="stability":
                    shutil.move(path, os.path.join("images\\stability",entry.name)) 
        it.close()
        
        
        
        
        
        
         
        
        
        
        
        
        
        
        
        