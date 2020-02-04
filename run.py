from Script.cif_op import cif_op
from Script.AnalysisTool import mofAnalysis
from Script.GaussianVib import GaussianVib
from Script.cif_expantion import cif_expantion
from Script.DataOut import DataOut
import Script.cluster_chargeBased as cluster_chargeBased
import numpy as np
import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt

def __init__(vibIn=None,vibAdd=None,outFile=None,mofName=None):
    sys.path.append('/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/Structures_Rigid&Snapshots')
    mofNames=os.listdir("/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/Structures_Rigid&Snapshots/.")
    sys.path.append('/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/LammpsData')
    lammps_dir='/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/LammpsData/'
    snapshotfile_dir='/home/zyu/Dropbox (GaTech)/simulation/Flexibility-Project/Shared_Zhenzi/Structures_Rigid&Snapshots/'
    
    #########################      position[-1]=rigid        ##################################

#def __VibInput__(self,vibIn,mofNum,snapshotfile_dir,mofNames):
    if vibIn=='snapshots' or vibIn=='cmp':
        fileNames=[]
        dataRange=np.arange(0,50000,1000)
        for i in dataRange:
            fileNames.append('snapshot_'+str(i))
        fileNames.append('RigidStructure')
            #print(mofNames[i])
        var1=cif_op(snapshotfile_dir,mofName,fileNames)
        pos,posFrac,charge,atomType,atomTypeStr,atomName,lattice,index,bondInfo=var1._batch_op_()
        index_associatedH,index_H=var1._findHList_(pos[-1],atomType)
        var2=mofAnalysis(np.array(pos[:-1],dtype=float),atomType,charge,lattice,atomTypeStr)
        fullRankStd,fullRankMean=np.array(var2._fullRankStdGen(),dtype=float)
    if vibIn=='charge':
        atom3Dstd,atom3DMean=var2._chargeBasedStdGen_()
    if vibIn=='cif' or vibIn=='cmp':
        var_cif=cif_expantion(mofName)
        posRigid,atomName,lattice=var_cif._symOp_()
    if vibIn=='cmp':
        chargeBasedstd,chargeBasedMean=var2._chargeBasedStdGen_()
        cifBasedVib=var_cif._CCSD_GenGaussian_()
        snapAtomName=np.array([atomTypeStr[np.argwhere(atomType==i+1)[0,0]][0] for i in range(0,max(atomType))])
        DataOut._outParityPlot4cmp_(snapAtomName,cifBasedVib,chargeBasedstd,mofName)
        raise ValueError('ok')

        
        
#def __VibAdd__(self,vibAdd):  
    if vibAdd=='atom':    
        var4=GaussianVib()
        pos_Gaussian_fract,pos_Gaussian_cartn=var4._fullRank_add_Gaussian_(pos[-1],fullRankStd,lattice,index_associatedH,index_H,50)
    elif vibAdd=='correlated':
        var2._bondCorrelation_(bondInfo,pos_Gaussian_cartn)    
    elif vibAdd=='cif':
        var_cif._CCSD_GenGaussian_()
    else:
        raise ValueError('vibOut not specified')
        
 
#def __OutPut__(sekf,outFile):
    if outFile=='cif':
        var1._out2cif_('outGaussian',pos_Gaussian_fract)
        var1._out2cif_('outSnapshot',posFrac)
    elif outFile=='xyz':
        var1._outxyz_(pos_Gaussian_cartn)
    else:
        raise ValueError('no Output Specified')
  
    
    return 


      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
           
        
#def __Others__(self):        
    # for unit cell
    switch=False
    if switch:
        corMofVar=cif_op(snapshotfile_dir,mofNames[i],fileNames)
        positionUnitCellOrdered,superCellNum,atomNuminUnitCell,indexAll=var1._symFind_(posFrac)
        var3=mofAnalysis(np.array(positionUnitCellOrdered[:-1],dtype=float),atomType,charge,lattice,atomTypeStr)
        a_std,a_mean=var3._chargeBasedStdGen_()  
    #test
    switch=False
    if switch:
        var12=mofAnalysis(np.array(pos_Gaussian_cartn,dtype=float),atomType,charge,lattice,atomTypeStr)
        atom3Dstd,atom3DMean=var2._chargeBasedStdGen_()
        fullRankStd2,fullRankMean2=np.array(var12._fullRankStdGen(),dtype=float)
        temp=(fullRankStd2/fullRankStd).ravel()
        
        ax=sns.distplot(temp,axlabel="std_Artifitial/std_Snapshot",norm_hist=True)
        plt.savefig('std_test.png',dpi=500)
