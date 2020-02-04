import sys
import numpy as np
from Script.DataOut import DataOut

np.set_printoptions(threshold=sys.maxsize)

class GaussianVib(object):
        def __init__(self): #No need to implement
             pass

        def _atom_type_add_Gaussian_(self,positionR,pattern,label,imgOutNum,cellParameter):
            #modifiedPositionCartn=np.ones((imgOutNum,np.shape(positionR)[0],np.shape(positionR)[1]))
            modifiedPositionFract=np.ones((imgOutNum,np.shape(positionR)[0],np.shape(positionR)[1]))
            for i in range(0,imgOutNum):
                for j in range(0,len(pattern)):
                    index_atom=np.where(label[:]==pattern[j][0])
                    for k in range(0,3):
                        vib=np.random.normal(0,np.sqrt(float(pattern[j][1]))/cellParameter[k],(np.size(index_atom),1)).T
                        #vib=np.random.normal(0,0/cellParameter[k],(np.size(index_atom),1)).T
                        modifiedPositionFract[i,index_atom,k]=positionR[index_atom,k]+vib
                    #np.sqrt(float(pattern[j][1]))
                
                #modifiedPositionFract[i]=np.matmul(np.linalg.inv(lattice),modifiedPositionCartn[i].T).T
            return modifiedPositionFract

        def _fullRank_add_Gaussian_(self,positionR,pattern,lattice,index_associatedH,index_H,imgOutNum):
            temp=np.arange(np.shape(positionR)[0])
            index_other=np.setdiff1d(temp,np.array(index_H))
            modifiedPositionFract=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            modifiedPositionCartn=np.ones((imgOutNum,np.shape(positionR)[0],3),dtype=float)
            for i in range(0,imgOutNum):
                for k in range(0,np.shape(index_other)[0]):
                    vib=np.random.normal(0,pattern[index_other[k]])
                    modifiedPositionCartn[i,index_other[k]]=positionR[index_other[k]]+vib
                    if index_associatedH[index_other[k]] < np.shape(positionR)[0]:
                           modifiedPositionCartn[i,index_associatedH[index_other[k]]]=positionR[index_associatedH[index_other[k]]]+vib
                            
                modifiedPositionFract[i]=np.matmul(np.linalg.inv(lattice),modifiedPositionCartn[i].T).T
            return modifiedPositionFract,modifiedPositionCartn
