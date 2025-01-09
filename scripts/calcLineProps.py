import numpy as np
import pandas as pd
import os




def calcEquivalentDia(lineWeightInAir, lineWeightInWater):

    '''
    Calculate the equivalent diameter and cross section area of a cylinder having the same upthrust per unit length as
    the line in water.
    
    Input - the weight in air and water (N/m)
    '''

    rhoSea = 1025 #kg/m3

    eqArea = (lineWeightInAir - lineWeightInWater) / rhoSea / 9.81 #m2

    eqDia = np.sqrt(4 * eqArea / np.pi) 
    return eqDia, eqArea



def calcEquivalentCoefficients(d1, d2, hydroCoeff):

    '''
    Find the equivalent hydro coefficent for d2 such that the hydro force is same as that given by combination of 
    d1 and hydroCoeff
    '''


    correctedCoeff = {}

    for term, value in hydroCoeff.items():
        
        if d2 !=1: #to account for marine growth correction on nominal diameter
            term = term.replace('nomDia', 'eqDia')

        if 'Drag' in term:
           correctedCoeff[term] = d1/d2 * value

        elif 'AddedMass' in term:
            correctedCoeff[term] = (d1/d2)**2 * value
    
    return correctedCoeff



def marineGrowthCorrection(dNom, dt_mg, lType,  hydroCoeff):

    '''
    Calculate the change in line properties due to marine growth as specified in DNV-OS-E301 Chapter 2 Section 2.8

    
    '''

    g               = 9.81 #m/s2
    rhoSea          = 1025 #kg/m3
    rhoMarineGrowth = 1325 #kg/m3 #DNVGL-ST-0437 2.4.11

    #DNV-OS-E301 Chapter 2 Section 2.8
    if lType == 'sl_chain':
        mu = 2

    elif lType == 'polyester':
        mu = 1
    
    #mass of marine growth per unit length - kg/m
    m_growth = np.pi/4 * (((dNom) + 2 * (dt_mg))**2 - (dNom)**2) * rhoMarineGrowth * mu

    #submerged weight of marine growth per unit length
    w_growth = m_growth * (1 - rhoSea/rhoMarineGrowth) * g #N/m

    #hydroCoeff correction factor on Nominal diameter
    hydroCorrFactor = (dNom + 2 * dt_mg) / dNom

    #Make corrections on the pristine values
    massDensity, subWeight, volEqDia, mbl, EA = nrelLineProps(dNom, lType)
    massDensity = massDensity + m_growth #kg/m
    subWeight = subWeight + w_growth #N/m
    volEqDia, _ = calcEquivalentDia( massDensity * g, subWeight)

    mgHydroCorrection = calcEquivalentCoefficients(hydroCorrFactor, 1, hydroCoeff[lType]) #hydro coefficent correction for marine growth, Note this is still on Nominal diameter
    
    return massDensity, subWeight, volEqDia, mbl, EA, mgHydroCorrection



def nrelLineProps(dNom, lType):

    '''
    dNom in mm - Line properties calculated as per NREL's equations

    Updated the NREL's equation with trimmed decimal places based on the design basis document on 21-05-2024 


    Lm is the mean load in % of MBL
    '''

    g = 9.81 #m/s2

    rho = 1025 #kg/m3

    match lType :

        case 'sl_chain':
            
            massDensity = 20e3 * dNom**2 #kg/m
            volEqDia    = 1.8 * dNom #m
            subWeight   = massDensity * g  - np.pi/4 * volEqDia**2 * rho * g #N/m
            mbl         = 9.11e2 * dNom + 1.21e9 * dNom **2 - 2.19e9 * dNom**3 #N
            EA          = 8.56e10 * dNom**2 - 3.93e7 * dNom **3 #N

        case 'polyester':

            massDensity = 679 * dNom**2 #kg/m
            # wAir = massDensity * g
            # spGravity = 1.38
            # subWeight = wAir * (1 - 1/spGravity)
            # volEqDia, _ =  calcEquivalentDia(wAir, subWeight) #mm
            volEqDia = 0.79 * dNom #m
            subWeight = massDensity * g  - np.pi/4 * volEqDia**2 * rho * g #N/m
            mbl = 308e6 * dNom**2
            EA = [14* mbl  , (11.615, 0.396)]  

    return massDensity, subWeight, volEqDia, mbl, EA


def corrosionCorrection(mbl, dcorr, dnew):
    # DNV-OS-E301 Chapter 2 Section 2 5.2.2 
    return mbl * (dcorr/dnew)**2



def calcLineProps(dNom, lType, dt_mg, dt_corr, printRes = False):

    '''
    TODO : Can add the option here to use other equations to get line properties instead of the NREL's equations
    '''

    
    #Hydro coefficents based on the Nominal diameter
    hydroCoeff = {}
    hydroCoeff['sl_chain'] = {   #studless chain link
                                'transDrag_nomDia'          : 2.4,              # DNV-OS-E301 Chapter 2 Section 1, Table 1
                                'longDrag_nomDia'           : 1.15,             # DNV-OS-E301 Chapter 2 Section 1, Table 1
                                'transAddedMass_nomDia'     : 1.8**2 * 1.0,     # BV NR493 Table 2
                                'longAddedMass_nomDia'      : 1.8**2 * 0.5,     # BV NR493 Table 2
                            }

    hydroCoeff['polyester'] = {   #fiber rope
                                'transDrag_nomDia'          : 1.6,              # DNV-OS-E301 Chapter 2 Section 1, Table 1
                                'longDrag_nomDia'           : 0,                # DNV-OS-E301 Chapter 2 Section 1, Table 1
                                'transAddedMass_nomDia'     : 1.1,              # BV NR493 Table 2
                                'longAddedMass_nomDia'      : 0.15,             # BV NR493 Table 2
                                }





    #calculate the line properties with marine growth
    #coefficents are based on nominal diameter
    massDensity, subWeight, volEqDia, mbl, EA, mgHydroCorrection =  marineGrowthCorrection(dNom, dt_mg, lType, hydroCoeff)
    
    #convert the hydro coefficents to equivalent diameter
    eqHydroCoeff = calcEquivalentCoefficients(dNom, volEqDia, mgHydroCorrection)


    if lType == 'sl_chain':
        #apply corrosion correction on MBL
        mbl = corrosionCorrection(mbl, dNom-dt_corr, dNom)


    if printRes:
    #print results
        print(f'mbl                 = {mbl:.3f} N')
        print(f'dNom                = {dNom:.3f}*1000 mm')
        print(f'Submerged weight    = {subWeight:.3f} N/m')
        print('\n')
        print('MoorDynInput:')
        print(f'volEqDia            = {volEqDia:.3f}*1000 mm')
        print(f'masDensity          = {massDensity:.3f} kg/m')
        if isinstance(EA, np.ndarray):
            print(f'Static EA           = {EA[0]:.3e} N')
            print(f'Dynamic EA          = {EA[1]:.3e} N')
        else :
            print(f'EA                  = {EA:.3e} N')
        for key, value in eqHydroCoeff.items():
            print(f'{key} = {value:.3f}')
    
    props = {}
    props['massDensity'] = massDensity
    props['subWeight'] = subWeight
    props['volEqDia'] = volEqDia
    props['mbl'] = mbl
    props['EA'] = EA
    resDict = props|eqHydroCoeff

    return resDict


def getMoorpyLineType(dNom, linTyp, dt_mg, dt_corr, printRes = False):
    '''
    Convenience function to get the line properties in the format required by moorpy

    dNom in meters
    linTyp : 'polyester' or 'sl_chain'
    dt_mg in meters - marine growth thickness
    dt_corr in meters - corrosion thickness for chain
    '''

    resDict = calcLineProps(dNom, linTyp, dt_mg, dt_corr, printRes)

    #get a name for the line type
    match linTyp:
        case 'sl_chain':
           linName = 'C'
        case 'polyester':
           linName = 'P' 
        case _:
           linName = lineType
    typestring = f"{linName}{dNom*1000:.0f}MG{dt_mg*1000:.0f}CR{dt_corr*1000:.0f}"
    
    d_vol = round(resDict['volEqDia'],4)
    mass  = round(resDict['massDensity'],4)
    w     = round(resDict['subWeight'],4)
    MBL   = round(resDict['mbl'],4)
    Cd    = round(resDict['transDrag_eqDia'],4)
    CdAx  = round(resDict['longDrag_eqDia'],4)
    Ca    = round(resDict['transAddedMass_eqDia'],4)
    CaAx  = round(resDict['longAddedMass_eqDia'],4)
    
    lineType = dict(name=typestring, d_vol=d_vol, m=mass,  w=w,
                MBL=MBL, input_d=dNom,
                Cd=Cd, CdAx=CdAx, Ca=Ca, CaAx=CaAx)
    
    #to account for the dynamic EA
    if isinstance(resDict['EA'], list):
        ##if its polyetser then it return a tuple with the constant and linear term as the dynamic EA
        lineType['EA'] = resDict['EA'][0]
        lineType['EAd'] = resDict['EA'][1][0] * MBL #this is the constant term in the dynamic EA
        lineType['EAd_Lm'] = resDict['EA'][1][1] * 100 #this is the linear term in the dynamic EA, need 
                                                       #to multiply by 100 because the NREL equation assume
                                                         #Lm in percentage but in Moorpy Lm is a fraction

    else:
        lineType['EA'] = resDict['EA']
        

    return lineType





if __name__ == '__main__':

    
    os.chdir(os.path.dirname(__file__))

    #generate a database of line properties
    diaListDict = {
                'polyester' : np.array([160, 168, 175, 183, 186, 190, 196, 200, 203, 210])/1000, 
                'sl_chain' : np.array([121, 125, 130, 135, 137, 142])/1000 
                }

    dt_mgList = np.array([0, 50, 100, 200])/1000 
    dt_corrList = np.array([0, 5, 10])/1000 
    lineTypeList = ['polyester', 'sl_chain']
    lineList = []
    Lm = 15 #for polyester line the mean load is assumed to be 15% of MBL


    for linTyp, diaList in diaListDict.items():
        for dNom in diaList:
            for dt_mg in dt_mgList:
                for dt_corr in dt_corrList:
                    if linTyp == 'polyester' and dt_corr != 0: #corrosion correction only for chain
                        continue
                    lineList.append({'lineType':linTyp, 'nomDia': dNom, 'dt_mg': dt_mg, 'dt_corr': dt_corr}|calcLineProps(dNom, linTyp, dt_mg, dt_corr, printRes = False))



    pd.DataFrame(lineList).to_csv(r"C:\GitRepos\openFastModels\fastSimulations\data\lineProps.csv", index = False)