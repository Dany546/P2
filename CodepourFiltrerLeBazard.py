import numpy as np
import scipy.io as spio
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg

"""
Loading the different files and storing them in a dictionnary. See the project 
instruction file available on moodle for a description of the data dictionaries. 

The folder dataNeuron and the folder dataMuscle must be in the same folder as 
this script.
dataP2 : ___ loadDataP2.py
        |___ dataNeuron : *.mat
        |___ dataMuscle : *.mat
               
"""

def loadMuscle():
    """
    Loading the different files and storing them in a dictionnary
    """
    
    
    HandX = spio.loadmat("dataMuscle/HandX.mat")
    HandXVel = spio.loadmat("dataMuscle/HandXVel.mat")
    HandXForce = spio.loadmat("dataMuscle/HandXForce.mat")
    
    HandY = spio.loadmat("dataMuscle/HandY.mat")
    HandYVel = spio.loadmat("dataMuscle/HandYVel.mat")
    HandYForce = spio.loadmat("dataMuscle/HandYForce.mat")
    
    Pectoralis = spio.loadmat("dataMuscle/Pectoralis.mat")
    Deltoid = spio.loadmat("dataMuscle/Deltoid.mat")
    
    extracted = spio.loadmat("dataMuscle/extracted.mat")
    descriptions = spio.loadmat("dataMuscle/descriptions.mat")
    
    
    """ 
    Creation of the first dictionnary - strMuscles
    """
    
    dictMuscles = {"HandX": HandX["HandX"],
                  "HandXVel": HandXVel["HandXVel"],
                  "HandXForce": HandXForce["HandXForce"],
                  "HandY": HandY["HandY"],
                  "HandYVel": HandYVel["HandYVel"],
                  "HandYForce": HandYForce["HandYForce"],
                  "Pectoralis": Pectoralis["Pectoralis"],
                  "Deltoid": Deltoid["Deltoid"],
                  "extracted": extracted["extracted"],
                  "descriptions": descriptions["descriptions"]}

    return dictMuscles

dictMuscles = loadMuscle()

def loadNeuron():
    
    namesSignals= [
        ('time'    ),
        ('shoang'  ),
        ('elbang'  ),
        ('handxpos'),
        ('handypos'),
        ('cells'   )]
    
    dictNeurons = {}
    for targetNum in range(1,9):
            
        target = {}
            
        for trialNum in range(1,7):
            trial = {}
            for nam in namesSignals:
                key = nam
                value = spio.loadmat('dataNeuron/target'+str(targetNum)+'trial' + str(trialNum) + 'signals'+nam+'.mat')
                trial[key]=value['a']
                
            target['trial'+str(trialNum)] = trial
    
        dictNeurons['target'+str(targetNum)] = target
        
    return dictNeurons
    
dictNeurons = loadNeuron()

#%%
#

"""
Dans le DictNeurons on a 8 targets, 6 essais pour chaque target, et 'cells' firing cells, 'elbang' elbow angle,
 'shoang' shoulder angle, 'handxpos' et 'handypos' (plot 2D?), et 'time' en ms?

"""

"""
Muscles ans le même sens que les essais
extracted 
1-60 ordre enregistrement
2e colonne
last colon is the reaction time

muscles
deltoid=> dos
pectoralis
"""

#Save all the figures from neuroMonkey for one feature ('cells','elbang' ,'shoang','handxpos', 'handypos') to put in string
def SaveImagesNeurons(string=None,dictN=dictNeurons):
    for i in range(1,9):
        for j in range(1,7):
            axY=None
            if string == "elbang":
                name = "Elbow angle"
                axY = "Angle [degree]"
            elif string=="shoang":
                name = "Shoulder angle"
                axY = "Angle [degree]"
            elif string=="cells":
                name = "Cells activity"
                axY=name
                axY= "Cells firing"
            elif string=="handxpos":
                name="Hand abscisse postion"
            elif string=="handypos":
                name="Hand ordonate position"
                axY=name
                
            title = name + ' Target %s'%i + ' Trial %s'%j
            dict_content = dictN['target%s'%i]['trial%s'%j][string]
            time = dictN['target%s'%i]['trial%s'%j]['time']
            if (string == 'shoang' or string=='elbang'):
                dict_content *= 180/np.pi #To put angles in degree
            plt.figure()
            plt.title(title)
            plt.ylabel(axY)
            plt.xlabel("Time [ms]")
            plt.plot(time,dict_content)
            plt.savefig('/Users/alexandreloffet/Documents/Unif/Master/M2Q1/Math model in neuroscience/Projet 2/loadDataP2/Graphs/Monkey/'+string+'/' + title)    

#Pour refaire toutes les images de base de dictNeurons
"""
listFeature = ['cells','elbang' ,'shoang','handxpos', 'handypos']
for i in listFeature:
    SaveImagesNeurons(i)
"""
def MovementApe(dictN=dictNeurons,HandX="handxpos",HandY="handypos"):
    for i in range(1,9):
        for j in range(1,7):                
            title = "Hand movement"  + ' Target %s'%i + ' Trial %s'%j
            handxpos= dictN['target%s'%i]['trial%s'%j][HandX]
            handypos=dictN['target%s'%i]['trial%s'%j][HandY]
            plt.figure()
            plt.title(title)
            plt.ylabel("Y position")
            plt.xlabel("X position")
            plt.plot(handxpos,handypos)
            plt.savefig('/Users/alexandreloffet/Documents/Unif/Master/M2Q1/Math model in neuroscience/Projet 2/loadDataP2/Graphs/Monkey/Movements/' + title)    


#MovementApe()
#%%
#Ici il faudrait faire qqch pour visualiser les données du muscle EMG + ...
    
    
    
    
    










#%%




















































#%%     TESTS à supprimer in fine
plt.close('all')
from scipy.ndimage import gaussian_filter as gauss            
for i in range(1,7):
    handypos=dictNeurons['target1']['trial%s'%i]["handypos"]
    time = dictNeurons['target1']['trial%s'%i]["time"]
    derivY = np.diff(handypos.flatten())
    FFTderiv = np.fft.fftshift(np.fft.fft(derivY))
    freq = np.fft.fftshift(np.fft.fftfreq(len(FFTderiv))) #shift?
    copFFT = np.zeros_like(FFTderiv)
    copFFT[:] = FFTderiv[:]
    copFFT[freq<0.01]*= 0
    """
    plt.figure()
    plt.plot(freq[freq>0],FFTderiv[freq>0])
    plt.plot(freq[freq>0],copFFT[freq>0])
    """
    derivFilt = np.fft.ifft(np.fft.ifftshift(copFFT)).real
    derivFilt=derivY -derivFilt
    
    GaussFilt = gauss(derivY,sigma=2)
    
    PeakofDev = sg.find_peaks(GaussFilt,height=0.1,width=20,distance=900) 
    print(PeakofDev[1]['peak_heights'])
    plt.figure()
    plt.title("Derivate of the Y hand position")
    plt.ylabel("Y position")
    plt.xlabel("Time [ms]")
    plt.plot(time[:-1],derivY)
    plt.figure()
    plt.plot(time[:-1],derivFilt)
    plt.figure()
    plt.plot(time[:-1],GaussFilt)
    plt.plot(time.flatten()[PeakofDev[0][:]],PeakofDev[1]['peak_heights'],'-or')
    
            
              
        
            
"""  
for i in range(1,7):
    T1 = dictNeurons['target1']['trial%s'%i]['shoang'] 
    
    T1*=(180/np.pi) #angles en degre
    plt.figure()
    plt.title("shoulder angle")
    plt.plot(T1)
    plt.plot()
    #plt.savefig(name + '.png', bbox_inches='tight')
"""








"""
#Faire une fonction pour chopper ts les ess
for i in range(1,7):
    T1 = dictNeurons['target1']['trial%s'%i]['shoang'] 
    T1*=(180/np.pi) #angles en degre
    plt.figure()
    plt.title("shoulder angle")
    plt.plot(T1)
    plt.plot()
""" 
"""
for i in range(1,7):
    T1 = dictNeurons['target1']['trial%s'%i]['elbang']
    T1*=(180/np.pi) #angles en degre
    plt.figure()
    plt.title("elbow angle")
    plt.plot(T1)
    plt.plot()
"""

"""
for i in range(1,7):
    T1 = dictNeurons['target1']['trial%s'%i]['cells'] 
    plt.title("cells")
    plt.plot(T1)
    plt.show()
"""

"""
#di = dictNeurons['target1']
i=1
T1 = dictNeurons['target%s'%i]['trial6']['shoang'] 
T1*=(180/np.pi) #angles en degre
plt.figure()
plt.plot(T1)
plt.plot()


HandX =  dictMuscles['HandX']

plt.figure()
plt.plot(HandX[0])
plt.plot()
"""
"""
HandX = dictNeurons['target1']['trial6']['handxpos'] 
HandY = dictNeurons['target1']['trial6']['handypos'] 
plt.figure()
plt.plot(HandX,HandY)
plt.show()
"""






