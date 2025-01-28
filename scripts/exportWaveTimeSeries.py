from openfast_toolbox.io import FASTOutputFile
import pandas as pd
from pathlib import Path
import multiprocessing as mp
import os
from functools import partial


def genExtFileList(resDirectory, ext):

    resDirectory = Path(resDirectory)
    
    fileList = []

    dirs = [entry for entry in resDirectory.iterdir() if entry.is_dir()] #also looks into the first subdirectory
    dirs.extend([resDirectory]) #also look ito the first level directory

    for directory in dirs:
        fileList.extend(list((directory.glob(f'*{ext}'))))
    
    fileList = [str(path) for path in fileList]

    return fileList



def expWaveElevation(file_path, resDir):

    df = FASTOutputFile(file_path).toDataFrame()
    df = df[['Time_[s]','Wave1Elev_[m]']]
    df.to_parquet(os.path.join(resDir, 'waveElev_' + os.path.basename(file_path).replace('.outb','.parquet')))

if __name__ == '__main__':

    resFol = r"W:\task49\Final\P200_C135_L1430_clump40_extreme"
    fileList = genExtFileList(resFol, '.outb')

    fileList = genExtFileList(resFol, '.outb')
    nCore = min(40, mp.cpu_count()-2, len(fileList))


    with mp.Pool(nCore) as pool:
        dictList = pool.map(partial(expWaveElevation, resDir = r'W:\task49\Final\P200_C135_L1430_clump40_extreme\exportWaveElevation'), fileList)










