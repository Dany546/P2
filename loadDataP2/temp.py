# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import matplotlib.pyplot as plt
import scipy.signal as sg
import scipy.ndimage as im
import numpy as np
import os 
import pywt  
from scipy.ndimage import gaussian_filter as gauss   
from scipy.signal import detrend

minimums = None ; aligned_cells = None
#%% Load data
from loadDataP2 import Load
  
os.chdir("G:\\Mon Drive\\Modelling_in_neurosciences\\Projet2\\loadDataP2")
ld = Load() 
neurons = ld.dictNeurons
muscles = ld.dictMuscles
print(neurons["target1"]["trial3"]["shoang"])
print(muscles["Deltoid"])

#%% Data exploration playgroud
print(muscles["extracted"][:,0])

#%% Data filtering

filtered = im.gaussian_filter(neurons["target1"]["trial3"]["shoang"].flatten(),sigma=1.5)
plt.plot(neurons["target1"]["trial3"]["shoang"].flatten(),lw=6)
plt.plot(filtered,',-',lw=1)
plt.show()


#%% muscles
time = np.arange(1200)
  
des = muscles["descriptions"][0] # idk wtf is this
ext = muscles["extracted"]
for patient in range(0,len(ext)//60): # travels through patients
    fig222, ax222 = plt.subplots(1,1)
    mvmts = [[],[],[]]
    for trial_type in [1,2,3]: # through trial types 
        trials = np.arange(patient*60,(patient+1)*60)[ext[patient*60:(patient+1)*60,2]==trial_type] # gets indices of trials belonging to the same type 
        # some variables initialisation
        stt = [] ; sD = 0 ; sP = 0
        summedP = 0
        summedD = 0
        count = 0
        for trial in trials: # travels through trials of type trial_type
         
        # gets corresponding variables 
            x, y = muscles["HandX"][trial],muscles["HandY"][trial]
            vx, vy = muscles["HandXVel"][trial],muscles["HandYVel"][trial]
            P, D = muscles["Pectoralis"][trial],muscles["Deltoid"][trial]
            Fx, Fy = muscles["HandXForce"][trial],muscles["HandYForce"][trial]
            theta = np.arctan2(np.diff(x),np.diff(y) ) * np.sqrt(np.diff(x)**2 + np.diff(y)**2)
            theta = np.append(theta[0],theta)
            
            if abs(theta.std())>0.075 or True:
                count += 1
                # store to plot afterwards
                stt.append(theta) 
                sP, sD = sP+P-P.mean(), sD+D-D.mean()
                # if we want to use fft 
                
                # plots variables & CWT
                colormap = "jet"
                scales = 10**np.linspace(0.5,1.55,200)
                dt = 0.001      
                wavelet = 'cmor1.5-0.5'
                coefP, freqsP = pywt.cwt(P,scales, wavelet,sampling_period=dt) 
                coefD, freqsD = pywt.cwt(D,scales, wavelet,sampling_period=dt) 
                summedP += np.abs(coefP)
                summedD += np.abs(coefD)
                coefP = np.abs(coefP) ; coefD = np.abs(coefD) 
                mvmts[trial_type-1].append([x,y])
                # fig = plt.figure()
                # ax1 = fig.add_subplot(321) ; ax2 = fig.add_subplot(322,sharey=ax1)  
                # ax3 = fig.add_subplot(323) ; ax4 = fig.add_subplot(324,sharey=ax3)
                # ax5 = fig.add_subplot(325) ; ax6 = fig.add_subplot(326,sharey=ax5)
                # ax1.plot(x)                ; ax2.plot(y)
                # ax3.plot(vx)               ; ax4.plot(vy)
                # ax5.plot(Fx)               ; ax6.plot(Fy)
                # plt.show() 
                
                # f = plt.figure()
                # ax0 = f.add_subplot(311)                                        ; ax0.plot(theta,label=trial)
                # ax7 = f.add_subplot(323)                                        ; ax8 = f.add_subplot(324)
                # ax9 = f.add_subplot(325)                                        ; ax10 = f.add_subplot(326, sharey=ax9)
                # ax7.plot(P)                                                       ; ax8.plot(D)
                # ax9.pcolormesh(time, freqsP, coefP, shading='auto',cmap=colormap) ; ax10.pcolormesh(time, freqsD, coefD, shading='auto',cmap=colormap) 
                # ax9.set_yscale('log')                                             ; ax10.set_yscale('log')  
                # f.legend()
            else:
                print("trial %s not shown" %trial)   
        
        if count!= 0 and True:
            pass
            """
            # trying to use the mean accross trial of the same type like in bioninstru
            mt = np.mean(stt,axis=0) ; stt = np.std(stt,axis=0)
            meanP = summedP/count
            meanD = summedD/count
            
            
            colormap = "jet" 
            dt = 0.001       
            coefP = np.abs(meanP)
            coefD = np.abs(meanD) 
            
            f = plt.figure()
            ax0 = f.add_subplot(311)                                        ; ax0.plot(mt)
            ax7 = f.add_subplot(312)                                        ; ax0.fill_between(x=time,y1=mt+stt,y2=mt-stt,alpha=0.2,color='lightblue')
            ax9 = f.add_subplot(313)                                         
            ax7.pcolormesh(time, freqsP, coefP, shading='auto',cmap=colormap) ; ax9.pcolormesh(time, freqsD, coefD, shading='auto',cmap=colormap)    
            ax7.set_yscale('log')                                             ; ax9.set_yscale('log')  
           
            
            colormap = "jet"
            scales = 10**np.linspace(0.5,1.55,200)
            dt = 0.001      
            wavelet = 'cmor1.5-0.5'
            coefP, freqsP = pywt.cwt(sP/count,scales, wavelet,sampling_period=dt) 
            coefD, freqsD = pywt.cwt(sD/count,scales, wavelet,sampling_period=dt) 
            coefP = np.log1p(np.abs(coefP))
            coefD = np.log1p(np.abs(coefD))
            
            f = plt.figure()
            ax1 = f.add_subplot(221) ; ax2 = f.add_subplot(222)  
            ax3 = f.add_subplot(223) ; ax4 = f.add_subplot(224,sharey=ax3) 
            ax1.plot(sP/count)               ; ax2.plot(sD/count)
            ax3.pcolormesh(time, freqsP, coefP, shading='auto',cmap=colormap)            ; ax4.pcolormesh(time, freqsP, coefD, shading='auto',cmap=colormap)     
            ax3.set_yscale('log')                                                        ; ax4.set_yscale('log')  
             
            print('minimal freq: %s'%freqsD.min(),'\n')"""
        else: 
            print("what the heck")
        # break # stops the loop after 1st trial
    x,y = 0
    ax222.plot(x,y)
    plt.show() 
    # break # stops the loop after patient 1
  
#%% Filtering don't look
# minimums = {}
# for gg in range(1,len(neurons)+1):
#     end = 7 # 7   
#     time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten()[[0,-1]] 
#     minimums['target%s'%gg] = {}
#     time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten().size
#     for j in range(1,7): 
#         time_range = np.minimum(time_range,neurons['target%s'%gg]['trial%s'%j]["time"].flatten().size) 
#     interval = time_range  
    
#     sigma = 1   
#     for i in range(1,end):
#         plt.close("all")
#         time = neurons['target%s'%gg]['trial%s'%i]["time"].flatten()   
#         lendiff = neurons['target%s'%gg]['trial%s'%i]["time"].flatten().size - interval
#         time_range = slice(50,-50-lendiff) 
#         time = time[time_range]
#         handypos = abs(neurons['target%s'%gg]['trial%s'%i]["handypos"])
#         handypos  =  handypos.flatten()[time_range]
#         handypos -=  handypos[:100].mean()
#         derivY = np.diff(handypos) 
#         FFTderiv = np.fft.fftshift(np.fft.fft(derivY))
#         freq = np.fft.fftshift(np.fft.fftfreq(len(FFTderiv))) #shift?
#         copFFT = np.zeros_like(FFTderiv)
#         copFFT[:] = FFTderiv[:]
#         copFFT[freq<0.01]*= 0 
#         derivFilt = np.fft.ifft(np.fft.ifftshift(copFFT)).real
#         derivFilt = derivY -derivFilt
#         handx = handypos[:]
        
#         grady = gauss(derivY,sigma=sigma) ; y = handypos
           
#         handypos = abs(neurons['target%s'%gg]['trial%s'%i]["handxpos"])
#         time = neurons['target%s'%gg]['trial%s'%i]["time"].flatten()[time_range]
#         handypos  =  handypos.flatten()[time_range]
#         handypos -=  handypos[:100].mean()
#         derivY = np.diff(handypos) 
#         FFTderiv = np.fft.fftshift(np.fft.fft(derivY))
#         freq = np.fft.fftshift(np.fft.fftfreq(len(FFTderiv))) #shift?
#         copFFT = np.zeros_like(FFTderiv)
#         copFFT[:] = FFTderiv[:]
#         copFFT[freq<0.01]*= 0 
#         derivFilt = np.fft.ifft(np.fft.ifftshift(copFFT)).real
#         derivFilt = derivY -derivFilt
        
#         gradx = gauss(derivY,sigma=sigma)
         
#         norm = np.sqrt(gradx**2 +grady**2)  
#         gradir = abs(np.diff(np.arctan2(handypos,handypos))*norm*180/np.pi)
        
#         colormap = "jet"
#         scales = 10**np.linspace(1.4,1.5,500)
#         dt = time[1]-time[0]  # 100 Hz sampling    
#         fig2, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,gridspec_kw={"height_ratios":[1,1,2,1]})
#         ax1.plot(time,handypos.flatten()) 
#         ax1.plot(time,handx.flatten()) 
#         ax2.plot(time[1:],norm)  
#         coef, freqs = pywt.cwt(norm,scales,'gaus1',sampling_period=dt/1000) 
#         finalcoef2 = coef
#         ax3.pcolormesh(time[1:], freqs, coef, shading='auto',cmap=colormap) 
#         ax3.set_yscale('log') 
#         coef, freqs = pywt.cwt(gradir,scales,'gaus1',sampling_period=dt/1000) 
#         ax4.pcolormesh(time[1:], freqs, coef, shading='auto',cmap=colormap) 
#         ax4.set_yscale('log') 
#         plt.show() 
#         plt.plot(finalcoef2.mean(0))
#         finalcoef2 = coef  
#         plt.plot(finalcoef2.mean(0)/60)
#         plt.show()
             
#         means = finalcoef2.mean(0) /5
#         mins = np.argmin(means) -1  
        
#         plt.figure()
#         plt.plot(time,handypos.flatten()) 
#         plt.plot(time,handx.flatten()) 
#         plt.plot(time[1:],-means) 
#         plt.plot(time[1:],norm*0.5*(handypos.max()+handx.max())) 
#         plt.plot(time[1:],gradir/3) 
#         plt.scatter(time[mins],handypos.flatten()[mins])
#         plt.scatter(time[mins],handx.flatten()[mins])  
#         plt.scatter(time[mins],-means[mins])   
#         plt.ylim([-15,70])
#         plt.show()     
#         minimums['target%s'%gg]['trial%s'%i] = {}
#         minimums['target%s'%gg]['trial%s'%i]['gauss1'] = mins 
        
#         fig2, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,gridspec_kw={"height_ratios":[1,1,2,1]})
#         ax1.plot(time,handypos.flatten()) 
#         ax1.plot(time,handx.flatten()) 
#         ax2.plot(time[1:],norm)  
#         coef, freqs = pywt.cwt(norm,scales,'gaus2',sampling_period=dt/1000) 
#         finalcoef2 = coef
#         ax3.pcolormesh(time[1:], freqs, coef, shading='auto',cmap=colormap) 
#         ax3.set_yscale('log') 
#         coef, freqs = pywt.cwt(gradir,scales,'gaus2',sampling_period=dt/1000) 
#         ax4.pcolormesh(time[1:], freqs, coef, shading='auto',cmap=colormap) 
#         ax4.set_yscale('log') 
#         plt.show() 
#         plt.plot(finalcoef2.mean(0))
#         # finalcoef2 = coef  
#         plt.plot(finalcoef2.mean(0)/60)
#         plt.show()
             
#         means = finalcoef2.mean(0) #/5
#         mins = np.argmin(means) -1  
#         thresh = -0.5*means.min()
#         mins = sg.find_peaks(-means,height=thresh,width=10,distance=100)[0] +1 
#         while len(mins)!=2:
#             if len(mins)<2:
#                 thresh -= 0.5*thresh
#             else:     
#                 thresh += 0.5*thresh 
#             mins = sg.find_peaks(-means,height=thresh,width=10,distance=100)[0] +1  
        
#         plt.figure()
#         plt.plot(time,handypos.flatten()) 
#         plt.plot(time,handx.flatten()) 
#         plt.plot(time[1:],-means) 
#         plt.plot(time[1:],norm*0.5*(handypos.max()+handx.max())) 
#         plt.plot(time[1:],gradir/3) 
#         plt.scatter(time[mins],handypos.flatten()[mins])
#         plt.scatter(time[mins],handx.flatten()[mins])  
#         plt.scatter(time[mins],-means[mins])   
#         plt.ylim([-15,70])
#         plt.show()    
#         # print(mins)
#         # print(finalcoef2[:,mins].std()) 
#         # print(freqs[where].size)  
#         minimums['target%s'%gg]['trial%s'%i]['gauss2_start'] = mins[0]
#         minimums['target%s'%gg]['trial%s'%i]['gauss2_end'] = mins[1]
        
#     # print(freqs.min(),freqs.max()) 
     
# for gg in range(1,len(neurons)+1):    
#     end = 7 # 7   
#     time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten()[[0,-1]]  
#     for i in range(1,end):
#         time_range[0] = np.minimum(time_range[0],neurons['target%s'%gg]['trial%s'%i]["time"].flatten()[0])
#         time_range[1] = np.minimum(time_range[1],neurons['target%s'%gg]['trial%s'%i]["time"].flatten()[-1]) 
#     interval = time_range  
#     fig, (ax1,ax2) = plt.subplots(2,1) 
#     ref = 0.5*(minimums['target%s'%gg]['trial%s'%1]['gauss1']+minimums['target%s'%gg]['trial%s'%1]['gauss2_start'])  
#     for i in range(1,end): 
#         quid = np.argmax(time==interval[:,None],axis=1)
#         time_range = slice(quid[0]+100,quid[1]-50) 
#         minimum = 0.5*(minimums['target%s'%gg]['trial%s'%i]['gauss1']+minimums['target%s'%gg]['trial%s'%i]['gauss2_start'])
#         shift = int(minimum-ref)  
#         print(shift)
#         if shift>0:
#             handypos = neurons['target%s'%gg]['trial%s'%i]["handypos"].flatten()[time_range][shift:]
#             handxpos = neurons['target%s'%gg]['trial%s'%i]["handxpos"].flatten()[time_range][shift:]
#         else:    
#             handypos = neurons['target%s'%gg]['trial%s'%i]["handypos"].flatten()[time_range][:shift]
#             handxpos = neurons['target%s'%gg]['trial%s'%i]["handxpos"].flatten()[time_range][:shift]
#         ax1.plot(handypos)
#         ax2.plot(handxpos)
#     plt.show()

#%% Alignement 
# maybe tester la même avec les angles 'shoang' pour voir si c'est pas mieux

minimums = {} 
for gg in range(1,len(neurons)+1):
    end = 7 # 7   
    minimums['target%s'%gg] = {} # initialisation
    
    # finds the length of the shortest trial
    time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten().size-100
    for j in range(1,7): 
        time_range = np.minimum(time_range,neurons['target%s'%gg]['trial%s'%j]["time"].flatten().size) 
    interval = time_range 
    
    sigma = 2 # for gaussian smoothing  
    for i in range(1,end):
        plt.close("all")
        
        time = neurons['target%s'%gg]['trial%s'%i]["time"].flatten()   
        lendiff = neurons['target%s'%gg]['trial%s'%i]["time"].flatten().size - interval
        time_range = slice(100,-100-lendiff) # crop to the length of the shortest -50 at both edges
         
        # gets varibles
        time = time[time_range]
        handypos = neurons['target%s'%gg]['trial%s'%i]["handypos"]  # smooth a bit to denoise
        elb = neurons['target%s'%gg]['trial%s'%i]["elbang"] # smooth a bit to denoise
        sho = neurons['target%s'%gg]['trial%s'%i]["shoang"] # smooth a bit to denoise
        
        # remove the initial position ==> to remove ?
        handypos -=  handypos[:50].mean() 
        handypos  =  handypos.flatten()[time_range]
        elb -=  elb[:50].mean() 
        elb  =  elb.flatten()[time_range]
        sho -=  sho[:50].mean() 
        sho  =  sho.flatten()[time_range]
        
        
        # derives & smoothes 2 times to get 2nd derivative 
        derivY = np.diff(handypos)   
        grady = gauss(derivY,sigma=2*sigma) ; y = handypos
        varang = np.diff(elb)**2 + np.diff(sho)**2
        
        # same with x
        handxpos = neurons['target%s'%gg]['trial%s'%i]["handxpos"] 
        time = neurons['target%s'%gg]['trial%s'%i]["time"].flatten()[time_range]
        handxpos -=  handxpos[:50].mean()
        handxpos  =  handxpos.flatten()[time_range]
        derivX = np.diff(handxpos)  
        
        gradx = gauss(derivY,sigma=2*sigma) 
         
        # compute 2nd derivative's norm ==> is it better to compute the norm & then derivate ? or derivate & then compute the norm ?
        norm = gauss(abs(np.diff(np.sqrt(gradx**2 + grady**2 + varang))),sigma=2*sigma) # *abs(np.arctan2(gradx,grady) )    
          
        # find peaks + plots & store the 1st one
        vect = np.sqrt(handxpos**2+handypos**2+elb**2+sho**2)
        mini = np.argmin(abs(vect-vect.std()))
        for i in range(14):
            mini = mini - 0.5*(vect[mini]-vect[:mini-50].mean())
            mini = int(mini)
        
        maxi = sg.find_peaks(norm,height=max(norm)*0.6)[0][0]
        plt.plot(time,vect)
        plt.plot(time[mini],vect[mini],'ro')
        plt.show()
               
        minimums['target%s'%gg]['trial%s'%i] = maxi
          
# for gg in range(1,len(neurons)+1):  
#     end = 7 # 7   
#     time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten().size
#     for j in range(1,7): 
#         time_range = np.minimum(time_range,neurons['target%s'%gg]['trial%s'%j]["time"].flatten().size) 
#     interval = time_range    
#     fig, (ax1,ax2) = plt.subplots(2,1) 
#     time = neurons['target%s'%gg]['trial%s'%1]["time"].flatten()
#     ref = (minimums['target%s'%gg]['trial%s'%1])  
#     for i in range(1,end):   
#         handypos = neurons['target%s'%gg]['trial%s'%i]["handypos"].flatten() 
#         handxpos = neurons['target%s'%gg]['trial%s'%i]["handxpos"].flatten() 
#         ax1.plot(handypos)
#         ax2.plot(handxpos)
#     plt.show()
    
#%% Cells 
aligned_cells = {}
for i in range(1,9): 
    gg = i 
    
    # finds the length of the shortest trial
    time_range = neurons['target%s'%gg]['trial%s'%1]["time"].flatten().size-100
    for j in range(1,7): 
        time_range = np.minimum(time_range,neurons['target%s'%gg]['trial%s'%j]["time"].flatten().size-100) 
    interval = time_range 
    
    # variable initialisation
    aligned_cells["target%s"%i] = {}
    aligned_cells["target%s"%i]['cells'] = 0
    aligned_cells["target%s"%i]['handy'] = 0
    aligned_cells["target%s"%i]['handx'] = 0
    aligned_cells["target%s"%i]['ang'] = 0
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,sharex=True) 
    ref = (minimums['target%s'%gg]['trial%s'%1])   # début du mvm pr essai 1 (pris comme référence)
    for j in range(1,7): 
        
        # là tu prends le temps, la diff de taille avec le + court et l'indice du début du mvmt
        time = neurons["target%s"%i]["trial%s"%j]["time"].flatten()
        lendiff = time.size - interval 
        minimum = (minimums['target%s'%gg]['trial%s'%j])
        
        shift = int(minimum-ref)  # décalage par rapport à la ref
        
        # on prend les variables en alignant tt sur la ref + on coupe 100 au début et 50 à la fin
        if shift>=0:  
            time_range = slice(100+shift,-lendiff-50)
            acti = np.append(np.zeros(shift),neurons["target%s"%i]["trial%s"%j]["cells"].flatten()[time_range])
            y = neurons["target%s"%i]["trial%s"%j]["handypos"].flatten()[time_range]
            y = np.append(np.repeat(y[0],shift),y) 
            x = neurons["target%s"%i]["trial%s"%j]["handxpos"].flatten()[time_range]
            x = np.append(np.repeat(x[0],shift),x)  
            ang = neurons["target%s"%i]["trial%s"%j]["shoang"].flatten()[time_range]
            ang = np.append(np.repeat(ang[0],shift),ang) 
        else:        
            time_range = slice(100,shift-lendiff-50)
            acti = np.append(neurons["target%s"%i]["trial%s"%j]["cells"].flatten()[time_range],np.zeros(-shift))
            y = neurons["target%s"%i]["trial%s"%j]["handypos"].flatten()[time_range]
            y = np.append(y,np.repeat(y[-1],-shift)) 
            x = neurons["target%s"%i]["trial%s"%j]["handxpos"].flatten()[time_range]
            x = np.append(x,np.repeat(x[-1],-shift)) 
            ang = neurons["target%s"%i]["trial%s"%j]["shoang"].flatten()[time_range]
            ang = np.append(ang,np.repeat(ang[-1],-shift)) 
          
        # on stocke tt ds le dico    
        aligned_cells["target%s"%i]['cells'] += acti 
        aligned_cells["target%s"%i]['handy'] += y 
        aligned_cells["target%s"%i]['handx'] += x 
        aligned_cells["target%s"%i]['ang'] += ang 
            
        # on fait de jolis plots
        theta = np.arctan2(y,x)*180/np.pi
        #dtheta = np.diff(theta) 
        dtheta = np.arctan2(np.diff(y),np.diff(x)) 
        ax1.plot(y)
        ax2.plot(x)  
        ax3.plot(theta) 
        # ax4.plot(time[1:],abs(dtheta)/abs(dtheta).max(),alpha=0.5,color='lightblue') 
        ax4.plot(acti,color="orange")
        
        # aligned_cells["target%s"%i] = aligned_cells["target%s"%i] + acti 
        
    plt.show() 


# fig = plt.figure()
# plt.plot(aligned_cells["target%s"%1])    
    
#%% Firing rate
global aligned_cells

time = neurons["target1"]["trial3"]["time"].flatten()
for target,dicta in aligned_cells.items(): # parcourt le dico target par target: à l'iter tu auras target='targeti', dicta=aligned_cells['targeti']
  
    acti = dicta['cells']    # prend l'acti de la target actuelle
    
    rect, ra = plt.subplots(1,1)
    gau, ga = plt.subplots(1,1)  
    
    for width in [10,15,20]: # ça donne mieux avec 15 nan ?
    
        size = acti.size; width_time = width*(time[1]-time[0])*0.001
        
        # firering rate comme en tp
        
        # fenetre rectangulaire
        ra.plot(acti) 
        kernel = np.concatenate((np.zeros(size//2 - width//2),np.ones(width)/width_time,np.zeros(size//2 - width//2)))
        fire = im.convolve(acti,kernel)*width_time 
        ra.plot(fire)  
        
        # kernel gaussien
        corr = 0.5
        sigma = width_time*corr 
        ga.plot(acti)
         
        kernel = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-((time-time.mean())*0.001)**2/(sigma**2))
        fire = im.convolve(acti,kernel)*width_time 
        ga.plot(fire) 
    plt.show() 