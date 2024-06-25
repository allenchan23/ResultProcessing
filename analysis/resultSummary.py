import os
import json
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import tstd
from rocketcea.cea_obj import CEA_Obj
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant

"""
Summary of test firing results to CSV File
"""
# Standard temperature
TEMP_STANDARD = 273.15
R_UNIVERSAL = 8314

# Result Summary Class
class PerformanceSummary:
    def __init__(self,test):
        self.test = test
        self.pOx = {}
        self.mdotOx = {}
        self.pOx = {}
        self.mdotF = {}
        self.pCo = {}
        self.mdotCo = {}
        self.mdotProp = {}
        self.mdotTot = {}
        self.ER = {}
        self.CR = {}
        self.P0 = {}
        self.P1 = {}
        self.P2 = {}
        self.P3 = {}
        self.Pb = {}
        self.F = {}
        self.cstar = {}
        self.tcstar = {}
        self.nstar = {}
        self.cISP = {}
        self.tISP = {}
        self.nISP = {}


# Representative value generator class
class timeHistoryCalculator():

    def __init__(self, target_df, time_column, t_start, t_end):
        self.target_df   = target_df
        self.time_column = time_column
        self.t_start     = t_start
        self.t_end       = t_end


    def get_time_average(self, target_columnm, allow_return_tstd_err=False):
        target_values = self.target_df[(self.target_df[self.time_column] >= self.t_start) & (self.target_df[self.time_column] <= self.t_end)][target_columnm]
        if allow_return_tstd_err:
            tstd_err = tstd(target_values, ddof=1)[0]
            return target_values.mean().iloc[0], tstd_err
        return target_values.mean().iloc[0]


# Injector Class
class Injector:
    def __init__(self ,injectorJSON):
        self.diameter = injectorJSON.get("diameter")
        self.cd = injectorJSON.get("cd")
        self.cdErr = injectorJSON.get("cdErr")
        self.count = injectorJSON.get("nOrifice")
        self.area = self.calInjArea()


    def calInjArea(cls):
        return (cls.diameter/2) ** 2 * np.pi * cls.count


    def getCd(cls):
        return cls.cd


    def getCdErr(cls):
        return cls.cdErr


    def getArea(cls):
        return cls.area

class Chamber:
    def __init__(self, chamberJSON):
        self.diameter = chamberJSON.get("diameter")
        self.At = np.pi*(self.diameter/2)**2

    def getThroatArea(cls):
        return cls.At

# Propellant Class
class Prop:
    def __init__(self, propellantJSON):
        self.type = propellantJSON.get("compound")
        self.gamma = propellantJSON.get("gamma")
        self.molMass = propellantJSON.get("molMass")
        self.gasConst = R_UNIVERSAL/propellantJSON.get("molMass")
        self.sigmaStar = self.calSigStar()


    def calSigStar(cls):
        """
        Calculate mass flow coefficient
        """
        gamma = cls.gamma
        return np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(gamma-1)))


    def getMolMass(cls):
        return cls.molMass


    def getGamma(cls):
        return cls.gamma


    def getSigmaStar(cls):
        return cls.sigmaStar


    def getR(cls):
        return cls.gasConst


def ftps2mps(speed_ftps):
    """
    Converts feet per second to meters per second
    """
    return speed_ftps*0.3048


def deg2kel(temp_deg, temp_stndard = TEMP_STANDARD):
    """
    Converts degree to kelvin
    """
    return temp_deg + temp_stndard


def kel2deg(temp_kel, temp_stndard = TEMP_STANDARD):
    """
    Converts kelvin to degree
    """
    return temp_kel - temp_stndard


def MPa2Pa(pMPa):
    """
    Converts MPa to Pa
    """
    return pMPa*1000000

def getStoichometricRatio(oxidizer, fuel):
    if (oxidizer.type == "O2") and (fuel.type == "C2H4"):
        FA_stoic = (fuel.getMolMass())/(3*oxidizer.getMolMass())
        return FA_stoic
    else:
        print("ERROR NO PROPELLANT COMBINATION")
        return 0


def calMdot(prop, inj, p, pErr, T):
    """
    Calculate experimental mass flow rate and error
    """
    # Get injector properties
    cd = inj.getCd()
    cdErr = inj.getCdErr()
    A = inj.getArea()

    # Get gas properties
    R = prop.getR()
    sigStar = prop.getSigmaStar()

    # Calcualte mdot and propagation error
    mdot = cd * p * A / np.sqrt(R * T) * sigStar
    mdotErr = calcPropErr(val = mdot, var1 = p, var2 = cd, var1Err = pErr, var2Err = cdErr)
    return mdot, mdotErr


def calcER(oxidizer, fuel, mdotOx, mdotOF, mdotOx_err, mdotF_err):
    """
    Calculate experimental equivalence ratio and error
    """
    FA_stoic = getStoichometricRatio(oxidizer,fuel)
    FA_exp= mdotOF/mdotOx
    equiv = FA_exp/FA_stoic
    equiv_err = calcPropErr(val = equiv, var1 = mdotOx, var2 = mdotOF, var1Err = mdotOx_err, var2Err = mdotF_err)
    return equiv, equiv_err


def calcCP(mdotTotal, mdotCoolant, mdotTotal_err, mdotCoolant_err):
    """
    Calculate experimental cooling percentage
    """
    if mdotCoolant == 0:
        coolantPer = 0
        coolantPer_err = 0
    else:
        coolantPer = mdotCoolant/mdotTotal
        coolantPer_err = calcPropErr(val = coolantPer, var1 = mdotTotal, var2 = mdotCoolant, var1Err = mdotTotal_err, var2Err = mdotCoolant_err)

    return coolantPer, coolantPer_err

def calcExpIsp():
    return 0

def calcIdealIsp():
    return 0

def calcExpCstar(chamber,m):
    return 0

def calcIdealCstar(chamber, mdot, pc, mdotT_err = 0, pc_er = 0):
    At = chamber.getThroatArea()
    cstar = pc*At/mdot
    cstar_err = calcPropErr(val = cstar, var1 = mdot, var2 = pc, var1Err = mdotT_err, var2Err = pc_er)
    return cstar, cstar_err

def calcPropErr(val, var1, var2, var1Err, var2Err):
    """
    Calculate standard deviation error of time values
    """
    val_err = np.sqrt((val*var1Err/var1)**2+(val*var2Err/var2)**2)
    return val_err


def performanceSummary(filePath, testName, stage):
    """
    Summaruse engiene performance from test
    """
     # Load file data
    if not os.path.isfile(filePath):
        print("     WARNING: " + testName + "_P.csv does not exist. skipping conversion...")
        return
    else:
        df = pd.read_csv(filePath)
        print("Summarising: " + testName)

        # Load test conditions
        conditionPath = os.path.join(Path(__file__).parent,"conditions",testName + ".json")
        conditions = json.load(open(conditionPath, encoding="utf-8"))

        # Overall Test Conditions
        temp = conditions.get("testConditions").get("temp")
        t0 = conditions.get("testConditions").get("Analysis_tS")
        t1 = conditions.get("testConditions").get("Analysis_tE")
        cooling = conditions.get("testConditions").get("filmCooling")
        detonation = conditions.get("testConditions").get("detonation")

        # Chamber Class
        chamber = Chamber(conditions.get("chamber"))

        # Injector Classes
        injectors = conditions.get("injectors")
        oxidizerInj = Injector(injectors.get("oxidizer"))
        fuelInj = Injector(injectors.get("fuel"))
        CoolantInj = Injector(injectors.get("coolant"))

        # Propellant Classes
        propellants = conditions.get("propellants")
        oxidizer = Prop(propellants.get("oxidizer"))
        fuel = Prop(propellants.get("fuel"))
        coolant = Prop(propellants.get("coolant"))

        # Load test file data
        df = pd.read_csv(filePath)

        # Initialise time hsitory calcualtors
        burnTimeTHC = timeHistoryCalculator(
            target_df=df,
            time_column="t",
            t_start=t0,
            t_end=t1
        )

        preBurnTHC = timeHistoryCalculator(
            target_df=df,
            time_column="t",
            t_start=-3,
            t_end=-2
        )


        # initialize reps_dict
        res_dict = {}

        # Oxygen Calculations
        pOx, pOx_err = burnTimeTHC.get_time_average(target_columnm=["POx"], allow_return_tstd_err=True)
        mdotOx, mdotOx_err = calMdot(prop = oxidizer, inj = oxidizerInj, p = MPa2Pa(pOx), pErr = MPa2Pa(pOx_err), T = deg2kel(temp))
        res_dict["POx_MPa"                       ] = pOx
        res_dict["POx_err_MPa"                   ] = pOx_err
        res_dict["MdotOx_kgps"                   ] = mdotOx
        res_dict["MdotOx_err_kgps"               ] = mdotOx_err

        # Fuel Calculations
        pF, pF_err = burnTimeTHC.get_time_average(target_columnm=["PF"], allow_return_tstd_err=True)
        mdotF, mdotF_err = calMdot(prop = fuel, inj = fuelInj, p = MPa2Pa(pF), pErr = MPa2Pa(pF_err), T = deg2kel(temp))
        res_dict["PFu_MPa"                       ] = pF
        res_dict["PFu_err_MPa"                   ] = pF_err
        res_dict["MdotFu_kgps"                   ] = mdotF
        res_dict["MdotFu_err_kgps"               ] = mdotF_err

        # Coolant Calculations
        if cooling:
            pCo, pCo_err = burnTimeTHC.get_time_average(target_columnm=["PC"], allow_return_tstd_err=True)
            mdotCo, mdotCo_err = calMdot(prop = coolant, inj = CoolantInj, p = MPa2Pa(pCo), pErr = MPa2Pa(pCo_err), T = deg2kel(temp))
        else:
            pCo = 0
            pCo_err = 0
            mdotCo = 0
            mdotCo_err = 0

        res_dict["Film_Cooling"                  ] = True
        res_dict["PCo_MPa"                       ] = pCo
        res_dict["PCo_err_MPa"                   ] = pCo_err
        res_dict["MdotCo_kgps"                   ] = mdotCo
        res_dict["MdotCo_Err_kgps"               ] = mdotCo_err

        # Total values
        mdotP = mdotOx + mdotF
        mdotP_err = mdotOx_err + mdotF_err
        mdotT = mdotP + mdotCo
        mdotT_err = mdotP_err + mdotCo_err

        # Equivalence Ratio
        ER, ER_err = calcER(oxidizer = oxidizer, fuel = fuel, mdotOx = mdotOx, mdotOF = mdotF, mdotOx_err = mdotOx_err, mdotF_err = mdotF_err)
        res_dict["ER"                           ] = ER
        res_dict["ER_err"                       ] = ER_err

        # Coolant Percentage
        CP, CP_err = calcCP(mdotTotal = mdotT, mdotCoolant = mdotCo, mdotTotal_err = mdotT_err, mdotCoolant_err = mdotCo_err)
        res_dict["CP"                           ] = CP
        res_dict["ER_err"                       ] = CP_err

        # Calculate pressures
        p0, p0_err = burnTimeTHC.get_time_average(target_columnm=["P0"], allow_return_tstd_err=True)
        p1, p1_err = burnTimeTHC.get_time_average(target_columnm=["P1"], allow_return_tstd_err=True)
        p2, p2_err = burnTimeTHC.get_time_average(target_columnm=["P2"], allow_return_tstd_err=True)
        p3, p3_err = burnTimeTHC.get_time_average(target_columnm=["P3"], allow_return_tstd_err=True)
        pb, pb_err = preBurnTHC.get_time_average(target_columnm=["PV"], allow_return_tstd_err=True)
        res_dict["P0_MPa"                       ] = p0
        res_dict["P0_err_MPA"                   ] = p0_err
        res_dict["P1_MPa"                       ] = p1
        res_dict["P1_err_MPA"                   ] = p1_err
        res_dict["P2_MPa"                       ] = p2
        res_dict["P2_err_MPA"                   ] = p2_err
        res_dict["P3_MPa"                       ] = p3
        res_dict["P3_err_MPA"                   ] = p3_err
        res_dict["Pb_MPa"                       ] = pb
        res_dict["Pb_err_MPA"                   ] = pb_err

        # Back thrust
        F, F_err = burnTimeTHC.get_time_average(target_columnm=["T"], allow_return_tstd_err=True)

        # calculate peformance values
        cstarEx, cstarEx_err = calcExpCstar(chamber, mdotT, MPa2Pa(p1), mdotT_err, MPa2Pa(p1_err))
        print(cstarEx, cstarEx_err)
        # cstarIsp, cstarIsp_err = calcExpCstar(chamber, mdotT, MPa2Pa(p1), mdotT_err, MPa2Pa(p1_err))



def write2file(summary):
    return


if __name__ == "__main__":

    # Set input files
    settings = json.load(open("conditions//settings.json", encoding="utf-8"))
    stage = settings.get("stage")
    tests = settings.get("tests")

    # Check existance of files
    parentDirectory = os.path.join(Path(__file__).parent.parent.parent, "data_processed",stage)
    print(parentDirectory)

    # Iterate through all test cases
    if "All" in tests:
        for t in sorted(os.listdir(parentDirectory)):
            # Prepare file and directory names
            fileDirectory = os.path.join(parentDirectory,t)
            filePath = str(os.path.join(fileDirectory,t)) + "_P.csv"
            testName = t
            # Summarise results
            summary = performanceSummary(filePath, testName, stage)

            # Save to summary cvs
            write2file(summary)

    elif len(tests) > 0:
        for t in tests:
             # Prepare file and directory names
            fileDirectory = os.path.join(parentDirectory,t)
            filePath = str(os.path.join(fileDirectory,t)) + "_P.csv"
            testName = t
            performanceSummary(filePath, testName, stage)

