import sys 
import numpy as np
import os
from Script.GaussianVib import GaussianVib
from scipy import spatial
np.set_printoptions(threshold=sys.maxsize)

class ExpBasedStr(object):
        def __init__(self,fileName,position,atomNames,index_H_bonded,occupancy,num,cellParameter): #No need to implement
            self.fileName=fileName
            self.num=num
            self.position=position
            self.atomName=atomNames
            self.occupancy=occupancy
            self.cellParameter=cellParameter
            
            self.atomType=np.loadtxt('exp_cifs/'+self.fileName+'/'+self.fileName+'.cifatomType',dtype=str)
            self.index_H_bonded=self._findHList_(self.position,self.atomType)
            
            

        def cat(self,outfilename, *infilenames):
            with open(outfilename, 'w') as outfile:
                for infilename in infilenames:
                    with open(infilename) as infile:
                        for line in infile:
                            if line.strip():
                                outfile.write(line)
            
        def _strGen_(self):
            var1=GaussianVib()
            position=self.position
            self.pos_Gaussian_fract=var1._atom_type_add_Gaussian_(position,self.pattern,self.atomName,self.num,self.cellParameter)
            return self.pos_Gaussian_fract
            
        
        def _read_(self):
            #self.positionR=np.loadtxt('cifFile'+'/'+self.fileName+'/'+self.fileName+'.cifpos')
            pattern=np.loadtxt('exp_cifs/'+self.fileName+'/'+self.fileName+'.cifpattern',dtype=str)
            self.pattern=self._Occupancy_(pattern)
            #self.label=np.loadtxt('cifFile'+'/'+self.fileName+'/'+self.fileName+'.ciflabel',dtype=str)
            #self.prePos=np.loadtxt('cifFile'+'/'+self.fileName+'/'+self.fileName+'.cifprePos',dtype=str)
            #self.afterPos=np.loadtxt('cifFile'+'/'+self.fileName+'/'+self.fileName+'.cifafterPos',dtype=str)
            #self.head=np.loadtxt('cifFile'+'/'+self.fileName+'/'+self.fileName+'.cifhead')
            return self.pattern
            
        def _write_(self):
            for i in range(0,self.num):
                posFilename = 'outGaussianVib_cifBased/'+self.fileName+'/'+str(i)+'.data'
                os.makedirs(os.path.dirname(posFilename), exist_ok=True)
                all_info=np.column_stack((self.prePos,self.pos_Gaussian_fract[i],self.afterPos))
                
                for k in range(0,np.shape(self.pos_Gaussian_fract)[1]):
                    print (*all_info[k,:], file = open(posFilename,"a"))
                
                headFileName='cifFile'+'/'+self.fileName+'/'+self.fileName+'.cifhead'
                outFileName='outGaussianVib_cifBased/'+self.fileName+'/'+str(i)+'.cif'
                self.cat(outFileName, headFileName, posFilename)
                  
            return
        
        
        def _Occupancy_(self,inputProperty):
            deletList=np.where(self.occupancy<0.5)
            outputProperty=np.delete(inputProperty,[deletList],0)
            return outputProperty
        
        def _findHList_(self,pos,atomType):
        
            index_H=np.where(self.atomType=='H')
            temp=np.arange(np.shape(pos)[0])
            index_other=np.setdiff1d(temp,np.array(index_H[0]))
            varH=spatial.KDTree(pos[index_H])
            varOther=spatial.KDTree(pos[index_other])
            result=varH.sparse_distance_matrix(varOther,2)
            index_nonzero=result.nonzero()
            
            #########based on H, but output in heavy atom#########
            index_associatedH=(np.shape(pos)[0])*np.ones(np.shape(pos)[0]+1,dtype=int)
            distanceArray=2*np.ones(len(index_H[0]),dtype=float)
            index_HAssociated2=(np.shape(pos)[0])*np.ones(len(index_H[0]),dtype=int)
            if len(index_nonzero[0]) < len(index_H[0]):
                print('may have serious H neignbour error!!')
            for i in range(0,len(index_nonzero[0])):
                if result[index_nonzero[0][i],index_nonzero[1][i]] < distanceArray[index_nonzero[0][i]]:
                    distanceArray[index_nonzero[0][i]]=result[index_nonzero[0][i],index_nonzero[1][i]]                 
                    index_HAssociated2[index_nonzero[0][i]]=index_other[index_nonzero[1][i]]
            for i in range(0,len(index_H[0])):
                index_associatedH[index_HAssociated2[i]]=index_H[0][i]
            
            return index_associatedH,index_H[0]