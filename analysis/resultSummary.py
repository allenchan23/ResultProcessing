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
class RepresentativeValuesGenerater():

    def __init__(self, target_df, time_column, t_start, t_end):
        self.target_df   = target_df
        self.time_column = time_column
        self.t_start     = t_start
        self.t_end       = t_end


    def get_time_average(self, target_columnm, allow_return_tstd_err=False):
        target_values = self.target_df[(self.target_df[self.time_column] >= self.t_start) & (df1[self.time_column] <= self.t_end)][target_columnm]
        if allow_return_tstd_err:
            tstd_err = tstd(target_values, ddof=1)[0]
            return target_values.mean().iloc[0], tstd_err
        return target_values.mean().iloc[0]


# Injector Class
class Injector:
    def __init__(self ,diameter, count, cd = 1, cdErr = 0):
        self.diameter = diameter
        self.cd = cd
        self.cdErr = cdErr
        self.count = count
        self.area = self.calInjArea()


    def calInjArea(cls):
        return (cls.diameter/2) ** 2 * np.pi * cls.count


    def getCd(cls):
        return cls.cd


    def getCdErr(cls):
        return cls.cdErr


    def getArea(cls):
        return cls.area


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


def calMdot(prop, inj, p, T):
    """
    Calculate mass flow rate
    """
    # Get injector properties
    cd = inj.getCd()
    A = inj.getArea()
    # Get gas properties
    gamma = prop.get("gasConst")
    sigStar = calSigStar(gamma)

    # Calcualte mdot
    mdot = cd * p * A / np.sqrt(gamma * T) * sigStar
    return mdot


def calSigStar(gamma):
    """
    Calculate mass flow coefficient
    """
    return np.sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)))


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
            # Open Pressure Time History
            fileDirectory = os.path.join(parentDirectory,t)
            filePath = str(os.path.join(fileDirectory,t)) + "_P.csv"
            testName = t
            performanceSummary(filePath, testName, stage)

    elif len(tests) > 0:
        for t in tests:
             # Open Pressure Time History
            fileDirectory = os.path.join(parentDirectory,t)
            filePath = str(os.path.join(fileDirectory,t)) + "_P.csv"
            testName = t
            performanceSummary(filePath, testName, stage)

