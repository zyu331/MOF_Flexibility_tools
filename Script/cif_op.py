import numpy as np
import pymatgen
import copy
from Script.DataOut import DataOut
from scipy import spatial
import pandas as pd

class cif_op:
    def __init__(self, snapshotfiledir,mofName,fileNames):
        self.snapshotfile_dir=snapshotfiledir
        self.mofName=mofName
        self.fileNames=fileNames
        self.ref=[]
        
    def _batch_op_(self):
        self.posBatch=[]
        self.posFracBatch=[]
        for i in range(0,len(self.fileNames)):
            posT,posFracT,charge,atomType,lattice,index=self._read_snapshots_LAMMPS_(self.fileNames[i])
            if i==0:
                self.ref=posFracT
            self.posBatch.append(posT)
            self.posFracBatch.append(posFracT)
        
        self._read_bondInfo_()
        self.atomTypeStr=[self.atomName[x-1][0:2] for x in atomType]
        return self.posBatch,self.posFracBatch,charge,atomType,self.atomTypeStr,self.atomName,lattice,index,self.bondInfo
    
    def _read_snapshots_LAMMPS_(self,fileName):
        File=self.snapshotfile_dir+self.mofName+'_clean_min_charges'+'/'+fileName
        with open(File, 'r') as the_file:
            all_data = [line.strip() for line in the_file.readlines()]
            temp=np.array([element.split(' ', 6) for element in all_data[9:]],dtype=float)
            data=temp.view(np.ndarray)[np.lexsort((temp[:, 0], ))]
        
            self.cell_info=np.array([line.split() for line in np.array(all_data[5:8])],dtype=float)
            self.box,self.charge,self.position,self.atomType,index=self.cell_info[:,1],data[:,5],data[:,2:5],np.array(data[:,1],dtype=int),np.array(data[:,0],dtype=int)
            self.positionFrac=self.position[:].copy()
            self._frac2Cartesian_()
            self._fetch_atom_name_()
            
        return (self.position,self.positionFrac,self.charge,self.atomType,self.lattice,index)
   
    def _read_bondInfo_(self):
        File='cifs/'+self.mofName+'_clean_min_charges'+'.bond'
        self.bondInfo=np.loadtxt(File)
        
   
    def _frac2Cartesian_(self):
        xy,xz,yz=self.cell_info[0,2],self.cell_info[1,2],self.cell_info[2,2]          
        lx,ly,lz=(self.cell_info[0,1]-max((0,xy,xz,xy+xz)))-(self.cell_info[0,0]-min((0,xy,xz,xy+xz))),(self.cell_info[1,1]-max(0,yz))-(self.cell_info[1,0]-min(0,yz)),self.cell_info[2,1]-self.cell_info[2,0]        
        a,b,c=lx,np.sqrt(ly**2+xy**2),np.sqrt(lz**2+xz**2+yz**2)
        self.box=[a,b,c]
        alpha,beta,gamma=(xy*xz+ly*yz)/(b*c),xz/c,xy/b
        omega=a*b*c*np.sqrt(1-alpha**2-beta**2-gamma**2+2*alpha*beta*gamma)            
# add PBC
        transM=np.array(([a,b*gamma,c*beta],[0,b*np.sin(np.arccos(gamma)), c*(alpha-beta*gamma)/np.sin(np.arccos(gamma))],[0,0,omega/(a*b*np.sin(np.arccos(gamma)))]))
        self.lattice=transM
        
        if np.shape(self.ref)[0]!=0:
            positionN=self._pbc_(self.position)
            self.position=positionN
            
        self.position=np.matmul(transM,(self.position).T).T
#for cif output
        self.cell_para_4cif=[a,b,c,np.rad2deg(np.arccos(alpha)),np.rad2deg(np.arccos(beta)),np.rad2deg(np.arccos(gamma))]
             
        
    def _pbc_(self,positionOld):
        self.box=np.array(self.box)
        size=np.shape(positionOld)
        positionN=self.position
        for j in range(0,size[0]):
            for k in range(0,size[1]):
                change=positionOld[j][k]-self.ref[j][k]
                positionN[j,k]=positionOld[j][k]-round(change)
                # posCartitial based: positionN[j,k]=positionOld[j][k]-round(change/self.box[k]*2)*self.box[k]/2
        return positionN
            
    def _fetch_atom_name_(self):

            self.filename2=self.mofName.rstrip('/')
            lammps_dir='/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/LammpsData/'
            with open(lammps_dir+'data.'+self.filename2+'_clean_min_charges', 'r') as the_file:
                all_data = [line.strip() for line in the_file.readlines()]
                index_l,index_h=all_data.index('Masses'),all_data.index('Bond Coeffs')                             
                self.atomName=[]
                data=[element.split(' ', 7) for element in all_data[index_l+2:index_h-1]]   
                for i in range(0,len(data)):    
                    nameT=data[i][-1]
                    self.atomName.append(nameT)
                    
    def _findHList_(self,pos,atomType):
        
        H_index=self.atomName.index('H_')
        index_H=np.where(atomType==H_index+1)
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
                    
    def _out2cif_(self,fileFolder,positionIn):
        DataOut._out2cif_(self,fileFolder,self.atomTypeStr, self.cell_para_4cif,positionIn, self.charge, self.filename2)
            
        return
    
    def _outxyz_(self,positionIn):
        DataOut._out2xyz_(self,self.atomTypeStr,positionIn,self.filename2,'tag')
        
        return
     
    def _symFind_(self,pos_thisf):
        
        #for i in range(0,np.shape(pos_thisf)[1]):
        #    self.charge[i]=self.charge[i]+10*i
        positionUnitCellOrdered=[]
        for i in range(0,len(self.fileNames)):
            ori_struc=pymatgen.core.structure.IStructure(self.lattice,self.atomType,pos_thisf[i],self.charge,validate_proximity=False,to_unit_cell=True,coords_are_cartesian=False)
            sys_struc=pymatgen.symmetry.analyzer.SpacegroupAnalyzer(ori_struc)
            lattice_matrix,chargeNew,struc=sys_struc.find_primitive()
            coord_struc=struc.frac_coords
           #sss=struc.charge
            
            allMatrix=np.matmul(lattice_matrix,coord_struc.T).T
            allMatrix=self._pbc_(allMatrix)
            positionUnitCellOrdered.append(copy.deepcopy(allMatrix))
            
            #coord_sturc2=sys_struc.get_conventional_standard_structure()
            
            #find index
            if i==0:
                indexAll=[]
                aT=np.round(np.matmul(self.lattice,pos_thisf[0].T).T,decimals=2)
                bT=np.round(allMatrix,decimals=2)
                for i in range(0,np.shape(pos_thisf)[1]):
                    indexAtom=np.where((aT==bT[i,:]).all(axis=1))
                    indexAll.append(indexAtom)
                  
        superCellNum=8            
        atomNuminUnitCell=int(np.shape(self.position)[0]/superCellNum)
            #output=np.savetxt(,ppos, fmt="%.3f")
        return positionUnitCellOrdered,superCellNum,atomNuminUnitCell,indexAll
    
            
        
        
    
    