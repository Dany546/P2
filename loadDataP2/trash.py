# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 09:07:16 2021

@author: rimez
"""

        # if count!= 0 and True:
        #     pass
        #     """
        #     # trying to use the mean accross trial of the same type like in bioninstru
        #     mt = np.mean(stt,axis=0) ; stt = np.std(stt,axis=0)
        #     meanP = summedP/count
        #     meanD = summedD/count
            
            
        #     colormap = "jet" 
        #     dt = 0.001       
        #     coefP = np.abs(meanP)
        #     coefD = np.abs(meanD) 
            
        #     f = plt.figure()
        #     ax0 = f.add_subplot(311)                                        ; ax0.plot(mt)
        #     ax7 = f.add_subplot(312)                                        ; ax0.fill_between(x=time,y1=mt+stt,y2=mt-stt,alpha=0.2,color='lightblue')
        #     ax9 = f.add_subplot(313)                                         
        #     ax7.pcolormesh(time, freqsP, coefP, shading='auto',cmap=colormap) ; ax9.pcolormesh(time, freqsD, coefD, shading='auto',cmap=colormap)    
        #     ax7.set_yscale('log')                                             ; ax9.set_yscale('log')  
           
            
        #     colormap = "jet"
        #     scales = 10**np.linspace(0.5,1.55,200)
        #     dt = 0.001      
        #     wavelet = 'cmor1.5-0.5'
        #     coefP, freqsP = pywt.cwt(sP/count,scales, wavelet,sampling_period=dt) 
        #     coefD, freqsD = pywt.cwt(sD/count,scales, wavelet,sampling_period=dt) 
        #     coefP = np.log1p(np.abs(coefP))
        #     coefD = np.log1p(np.abs(coefD))
            
        #     f = plt.figure()
        #     ax1 = f.add_subplot(221) ; ax2 = f.add_subplot(222)  
        #     ax3 = f.add_subplot(223) ; ax4 = f.add_subplot(224,sharey=ax3) 
        #     ax1.plot(sP/count)               ; ax2.plot(sD/count)
        #     ax3.pcolormesh(time, freqsP, coefP, shading='auto',cmap=colormap)            ; ax4.pcolormesh(time, freqsP, coefD, shading='auto',cmap=colormap)     
        #     ax3.set_yscale('log')                                                        ; ax4.set_yscale('log')  
             
        #     print('minimal freq: %s'%freqsD.min(),'\n')"""
        # else: 
        #     print("what the heck")
#%%% Filtering don't look
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