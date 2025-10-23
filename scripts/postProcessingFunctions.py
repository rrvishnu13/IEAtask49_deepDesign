import pCrunch as pc
from pathlib import Path
import shutil
import pickle
import os

def genExtList(resDirectory, extList):
    
    fileList = []

    dirs = [entry for entry in resDirectory.iterdir() if entry.is_dir()] #also looks into the first subdirectory
    dirs.extend([resDirectory]) #also look ito the first level directory

    if isinstance(extList, str):
        extList = [extList]

    for directory in dirs:
        for ext in extList:
            fileList.extend(list((directory.glob(f'*{ext}'))))
    
    fileList = [str(path) for path in fileList]

    return fileList



def replaceValueInFile(fileName, text2replace, newValList, newFileNameList, replaceSafeFlag = False):
    '''
        Can put each new value as a seperate file to make different varients of the same file

        replaceSafeFlag - if you are replacing a file, turn it on so that it can store a copy of the original file
    
    '''
        

    with open(fileName, 'r') as file :
        filedata = file.readlines()

    linList = []
    for ind, line in enumerate(filedata) :
        if text2replace in line:
            linList.append(ind)

    print(f'Found {len(linList)} instances of {text2replace} in the {fileName}')

    if len(linList) == 0:
        return None


    for newVal, newFileName in zip(newValList, newFileNameList):

        if fileName == newFileName:
            if replaceSafeFlag == True:
                #make a copy of the filedata
                shutil.copy(file, file + '_original')


        fildeDataCopy = filedata.copy()

        for ind in linList:
            fildeDataCopy[ind] = fildeDataCopy[ind].replace(text2replace, f'{newVal}')

        with open(newFileName, 'w') as file:
            file.writelines(fildeDataCopy)
        

def processOutputs(resDirectory, ext, nCore, pqFile = None, extreme_channels = [], fatigue_channels = {}, magnitude_channels ={}, trans = 0):

    '''
    Notes:
        - Make sure the outb files have differnet names, else it will be overwritten 
        - The output files can be on the first level under resDirectory or in the subdirectory
    '''

    outbFileList = genExtList(resDirectory, ext)
        

    la =pc.LoadsAnalysis(
                            outputs             =   outbFileList,               # The primary input is a list of output files
                            fatigue_channels    =   fatigue_channels,  
                            trim_data           =   (trans, ),                  # If 'trim_data' is passed, all input files will
                            extreme_channels    =   extreme_channels,
                            magnitude_channels  =   magnitude_channels
                        )

    if len(fatigue_channels) == 0:
        return_damage = False

    else:
        return_damage = True


    la.process_outputs(cores=min(nCore, len(outbFileList)), return_damage = return_damage)

    #summary stats
    res = {}
    res['summary_stats'] = la.summary_stats
    la.summary_stats.to_parquet(pqFile) 

    #extreme events
    if extreme_channels !=[]:
        with open(f'{pqFile}_extremeEvents.pickle', 'wb') as f:
            pickle.dump(la.extreme_events, f)
        res['extreme_events'] = la.extreme_events
    
    #damage
    if return_damage == True:
        dir_name, file_name = os.path.split(pqFile)
        base_name, extension = os.path.splitext(file_name)
        la.damage.to_parquet(os.path.join(dir_name, base_name + '_damage' + extension))
        res['damage'] = la.damage

    
    return res



if __name__ == "__main__" :
    #fatigue simulations

    # folders = [Path(r'W:\task49\fast_v4\fatigue\P200_C135_L1430_clump40_fat_-50')]
    # folders = [Path(r'W:\task49\fast_v4\fatigue\P200_C142_L1430_clump40_fat_0_v3')]
    folders = [f for f in Path(r'M:\test').iterdir() if f.is_dir()]

    #outb file processing
    ext = 'outb'
    nCore = 50

    extreme_channels = [] #['AnchTen4', 'LINE4N0FX', 'LINE4N0FY', 'LINE4N0FZ'] 
    fatigue_channels = {}
    magnitude_channels = {}
    trans = 600

    chainMbl = 16873826.79 #chain 142 mm dia + 5 mm corrosion
    # chainMbl = 15452629.04 #chain 135 mm dia + 5 mm corrosion
    # chainMbl = 14457922 #chain 130 mm dia + 5 mm corrosion
    chainIDList = [1, 6, 7, 12, 13, 18] #list the IDs of lines which are chain - both ends of these are considered for fatigue check
    chainFatParam = pc.FatigueParams(lifetime=25, load2stress=1/chainMbl, slope=3, S_intercept=316**(1/3))
    # chainFatParam = pc.FatigueParams(load2stress=1/chainMbl, slope=3, S_intercept=316**(1/3))

    #add fairlead tensions
    fatigue_channels = {f'FAIRTEN{chainID}': chainFatParam for chainID in chainIDList}
    fatigue_channels.update({f'ANCHTEN{chainID}': chainFatParam for chainID in chainIDList})

    anchLineIDList = [1, 7, 13]
    magnitude_channels = {f'ANCH{anchLineID}HORZ': [f'LINE{anchLineID}N0FX', f'LINE{anchLineID}N0FY'] for anchLineID in anchLineIDList}

    for resDirectory in folders :
        pqFile = resDirectory.joinpath(os.path.basename(resDirectory) + '.parquet')
        sum = processOutputs(resDirectory, ext, nCore, pqFile, extreme_channels, fatigue_channels, magnitude_channels, trans)



