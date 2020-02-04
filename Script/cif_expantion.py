import sys 
import numpy as np
import pandas as pd
from Script.ExpBasedStr import ExpBasedStr
#from Script.MentalLine import MentalLine
from Script.DataOut import DataOut
np.set_printoptions(threshold=sys.maxsize)

class cif_expantion(object):
        def __init__(self,cifName): #No need to implement
            self.cifName=cifName
            self.fileName='exp_cifs/'+self.cifName+'/'+self.cifName+'.cif'
            
#delete the headings
            
        def _symOp_(self):
            sym=pd.read_csv(self.fileName+'sym',sep=', ')
            position=np.loadtxt(self.fileName+'pos')
            self.index_H_bonded=[0,0,0,0,1,0,1,0,0,0]
            self.nonzero_index=[i for i, e in enumerate(self.index_H_bonded) if e != 0]
            self.atomName=np.loadtxt(self.fileName+'label',dtype=str)
            self.cellParameter=np.loadtxt(self.fileName+'cellPara')
            self.occupancy=np.loadtxt(self.fileName+'occupancy')
            
            position=self._Occupancy_(position)
            self.atomName=self._Occupancy_(self.atomName)
        
            x=position[:,0]
            y=position[:,1]
            z=position[:,2]
            
            atomNames=self.atomName
            index_H_bonded=self.index_H_bonded
            
            self._output_(position,atomNames,'ligand')
            
            for i in range(0,np.shape(sym)[0]):
                xT=eval(sym.iloc[i]['x'])
                yT=eval(sym.iloc[i]['y'])
                zT=eval(sym.iloc[i]['z'])    
                temp=np.stack((xT,yT,zT), axis=-1)
                position=np.append(position[:],temp,axis=0)
                atomNames=np.append(atomNames[:],self.atomName,axis=0)
                index_H_bonded=np.append(index_H_bonded,self.index_H_bonded,axis=0)
            for i in range(0,len(index_H_bonded)):
                if index_H_bonded[i]:
                    index_H_bonded[i]+=i
                    
            #self.position=np.expand_dims(position, axis=0)
            self._uniquePos_(position)
            self.position=np.delete(position,self.dup_label,0)[:,:]
            
            self.atomNames=np.delete(atomNames,self.dup_label,0)
            self.index_H_bonded=index_H_bonded
            
            self._output_(self.position,self.atomNames,'Rigid')
            
            return self.position,self.atomNames,self.cellParameter
        
        def _Occupancy_(self,inputProperty):
            deletList=np.where(self.occupancy<0.5)
            outputProperty=np.delete(inputProperty,[deletList],0)
            return outputProperty
        
        def _uniquePos_(self,position):
            self.dup_label=[]
            positionPBC=np.round(position-0.5)+position
            for i in range(0,np.shape(position)[0]):
                for j in range (0,i):
                    if np.array_equal(positionPBC[i,:],positionPBC[j,:]):
                        self.dup_label.append(i)
            self.dup_label=np.unique(self.dup_label)
            
            return self.dup_label
        
        def _CCSD_GenGaussian_(self):
            
            posCartn=self._frac2Cartesian_(self.position)
            num=49
            # decide cartn input or fract input
            var2=ExpBasedStr(self.cifName,self.position,self.atomNames,self.index_H_bonded,self.occupancy,num,self.cellParameter)
            pattern=var2._read_()
            
            #position_gaussian=var2._strGen_()
            #self._output_(position_gaussian,self.atomNames,'Guassian')
            return pattern
            
            
        def _output_(self,position,atomNames,outputFilename):
            var=DataOut()
            var._out2cif_('cifBasedOut/'+self.cifName,atomNames, self.cellParameter, position, 'none',outputFilename)
            
            return
        
        def _frac2Cartesian_(self,positionIn):
            M=np.zeros([6])
            M[0:3]=self.cellParameter[0:3]
            M[3:6]=self.cellParameter[3:6]*np.pi/180

            omega=M[0]*M[1]*M[2]*np.sqrt(1-np.cos(M[3])**2-np.cos(M[4])**2-np.cos(M[5])**2+2*np.cos(M[3])*np.cos(M[4])*np.cos(M[5]))
            self.transM=[[M[0],M[1]*np.cos(M[5]),M[2]*np.cos(M[4])],
                     [0,M[1]*np.sin(M[5]),M[2]*(np.cos(M[3])-np.cos(M[4])*np.cos(M[5]))/np.sin(M[5])],
                     [0,0,omega/(M[0]*M[1]*np.sin(M[5]))]]
            positionOut=np.matmul(self.transM,(self.position).T).T
            
            return positionOut
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
#        def _uniquePos_(self,position):
#            self.dup_label=[]
#            for i in range(0,np.shape(position)[0]):
#                for j in range (0,i):
#                    if position[i,0]==position[j,0] or (np.abs(position[i,0])+np.abs(position[j,0]))==1 or np.abs((np.abs(position[i,0])-np.abs(position[j,0])))==1:
#                        if position[i,1]==position[j,1] or (np.abs(position[i,1])+np.abs(position[j,1]))==1 or np.abs((np.abs(position[i,1])-np.abs(position[j,1])))==1:
#                            if position[i,2]==position[j,2] or (np.abs(position[i,2])+np.abs(position[j,2]))==1 or np.abs((np.abs(position[i,2])-np.abs(position[j,2])))==1:
#                                self.dup_label.append(i)
#            self.dup_label=np.unique(self.dup_label)
#            
#            return self.dup_label