import numpy as np
import scipy.io as spio
import csv
import numpy as np
import matplotlib.pyplot as plt
import os

"""
Loading the different files and storing them in a dictionnary. See the project 
instruction file available on moodle for a description of the data dictionaries. 

The folder dataNeuron and the folder dataMuscle must be in the same folder as 
this script.
dataP2 : ___ loadDataP2.py
        |___ dataNeuron : *.mat
        |___ dataMuscle : *.mat
               
"""
class Load:
    
    def __init__(self): 
        old = os.getcwd()
        os.chdir("G:\\Mon Drive\\Modelling_in_neurosciences\\Projet2\\loadDataP2")
        self.dictMuscles = self.loadMuscle()
        self.dictNeurons = self.loadNeuron() 
        os.chdir(old)
    
    def loadMuscle(self):
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
    
    
    
    def loadNeuron(self):
        
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
        
    
