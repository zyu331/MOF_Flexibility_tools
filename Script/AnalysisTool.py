import sys 
import numpy as np
from sklearn.cluster import DBSCAN
import pymatgen
np.set_printoptions(threshold=sys.maxsize)
import seaborn as sns
import matplotlib.pyplot as plt

class mofAnalysis(object):
        def __init__(self,position,atomType,charge,lattice,atomTypeStr): #No need to implement
             self.position=position
             self.atomType=atomType
             self.atomTypeNum=max(atomType)
             self.charge=charge
             self.chargeBasedType=charge
             self.lattice=lattice
             self.atomTypeStr=atomTypeStr
    
        def _chargeBasedTypeGen_(self,index_atom,acc=0.01):
    
            clustering = DBSCAN(eps=acc, min_samples=3).fit(self.charge[index_atom].reshape(-1, 1))
            e = clustering.labels_
            self.chargeBasedType[index_atom]=np.array(e,dtype=str)
            return self.chargeBasedType
       
        def _atom_type_get_parameter_(self,index_atom):
# calculate mean and std  #################################specified#########################
            atom3DMeanOneType=[]
            atom3DstdOneType=[]
            for j in range(0,self.atom_typeNum_chargeBased+1):
                index=[x for x in index_atom if self.chargeBasedType[x] == j]
                atom3DMeanT=np.mean(self.position[:,index,:],axis=0)
                atom3DstdT=np.std(self.position[:,index,:],axis=0)
                atom3DstdOneType.append(np.mean(atom3DstdT,axis=0))
                atom3DMeanOneType.append(np.mean(atom3DMeanT,axis=0))          
            self.atom3Dstd.append(np.array(atom3DstdOneType))
            self.atom3DMean.append(np.array(atom3DMeanOneType))
                
        def _chargeBasedStdGen_(self):
 #initial for chargeBased argrithm
            self.atom3Dstd=[]
            self.atom3DMean=[]

            for j in range(0,self.atomTypeNum):
                index_atom=np.where(self.atomType==j+1)
                self._chargeBasedTypeGen_(index_atom)
                self.atom_typeNum_chargeBased=int(max(self.chargeBasedType[index_atom]))
                # for each cahrge type given atom type
                self._atom_type_get_parameter_(index_atom[0])
            
            return self.atom3Dstd,self.atom3DMean
       
        def _fullRankStdGen(self):
            fullStd=[]
            fullMean=[]
            for i in range(np.shape(self.position)[1]):
                atom3DMeanT=np.mean(self.position[:,i,:],axis=0)
                atom3DstdT=np.std(self.position[:,i,:],axis=0)
                fullStd.append(atom3DstdT)
                fullMean.append(atom3DMeanT)
            
            return fullStd, fullMean
        
        def _fullRankPairwiseAnalysis_(self):
            for i in range(0,len(self.postion)):
                np
            
            return
        
        def _getNeighborList_(self,posRF):
            #expansion
            expandedPos=posRF
            expandedIndex=[]
            for i in range(0,3):
                pbcAtoms=np.where(posRF[:,i]<0.1)
                expandedPosT=posRF[pbcAtoms[0]]
                expandedPosT[:,i]=expandedPosT[:,i]+1
                expandedPos=np.concatenate((expandedPos, expandedPosT), axis=0)
                expandedIndex=np.concatenate((expandedIndex[:], pbcAtoms[0]), axis=0)
                
                pbcAtoms=np.where(posRF[:,i]>0.9)
                expandedPosT=posRF[pbcAtoms[0]]
                expandedPosT[:,i]=expandedPosT[:,i]-1
                expandedPos=np.concatenate((expandedPos, expandedPosT), axis=0)
                expandedIndex=np.concatenate((expandedIndex[:], pbcAtoms[0]), axis=0)
                
            #find neighbor
         
            return expandedPos
            
        def _bondCorrelation_(self,bondInfo,pos):
                      
            pairR=np.ones((np.shape(bondInfo)[0],np.shape(pos)[0]))
            for i in range(0,np.shape(pos)[0]):
                for j in range(0,np.shape(bondInfo)[0]):
                    pairR[j,i]=np.linalg.norm(pos[i][int(bondInfo[j,2]-1),:]-pos[i][int(bondInfo[j,3]-1),:])
            
            deletePBC = np.ma.masked_array(pairR, mask = (pairR >3))
            typeIndex=[]
            std=np.ones(int(max(bondInfo[:,1])))
            mean=np.ones(int(max(bondInfo[:,1])))
            for i in range(0,int(max(bondInfo[:,1]))):
                typeIndexT=np.where(bondInfo[:,1]==i+1)[0]
                std[i]=np.std(deletePBC[typeIndexT])
                mean[i]=np.mean(deletePBC[typeIndexT])
                typeIndex.append(typeIndexT)
            
            for i in range(0,30):
                plt.clf()
                sns.distplot(pairR[i,:])
                plt.xlim(1, 3)
                plt.savefig("figure/"+str(i)+".jpg")
                
            return pairR
            
               
        #type_chargeBased_4output=np.zeros((np.shape(positionR)[0],2))
        #type_4cifOut=np.ones((np.shape(positionR)[0],1))
        #type_4cifOut=np.array(type_chargeBased,dtype=str
            # give back value
            #atomTypeArray=j*np.ones((len(index_atom[0]),1))
            #type_chargeBased_4output[index_atom,:]=np.column_stack((atomTypeArray,type_chargeClustered))
            #df=pd.DataFrame({'col':type_chargeClustered})
            #type_chargeClustered=atomNames[j][0:2]+df.astype(str)
            #type_4cifOut[index_atom,:]=np.array(atomNames[j][0:2],dtype=str)
            #type_chargeBased[index_atom]=np.array(type_chargeClustered,dtype=str)
                
            # call output
            #d._out2xyz_(type_chargeBased,position,filenames[i],'snapshots') 