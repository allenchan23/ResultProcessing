import numpy as np
import pandas as pd
from pathlib import Path
import os
import json

import matplotlib.pyplot as plt
import scienceplots

"""
Conversion of raw csv data to values for data processing
"""

def processTemps(inputFilePath, outputFilePath,shotName):
    """
    Temperature file conversion
    """
    # Get file information
    specs = json.load(open(os.path.join("file_spec",(shotName + ".json")), encoding="utf-8")).get("temperature")
    fileName = os.path.join(inputFilePath,(shotName + "_T.csv"))
    sensorArray = specs.get("header")
    err = specs.get("errors")

    # Load file data
    if not os.path.isfile(fileName):
        print("     WARNING: " + shotName + "_T.csv does not exist. skipping conversion...")
        return
    else:
        df = pd.read_csv(fileName, names=sensorArray)

    # Check if data needs to be reduced
    if specs.get("compress"):
        df = df.iloc[::10,:]

    # Convert to celcius
    if specs.get("units") == "K":
        for s in sensorArray:
            if s != "t":
                df[s] = df[s] - 273.15

    # Revrse incorrect thermal couples orders
    for s in sensorArray:
        if (s != "t") and (s not in err) and (df[s] < -df.iloc[0][s]*1.25).any():
            print("     WARNING: " + shotName + ": Thermocouple " + s + " is reversed")
            # Get temepratures from other data points
            temps = []
            for o in sensorArray:
                if (o != s) and (o not in err) and (o != "t"):
                    temps.append(df.iloc[0][o])
            # Correct data
            T0_est = sum(temps)/len(temps)
            df[s] = T0_est+(df.iloc[0][s] -df[s])

    # Remove error values
    for e in err:
        # Set NaN to all values
        if len(err.get(e)) == 0:
            df[e] = np.nan
        # average error values
        elif len(err.get(e)) == 1:
            c = df.columns.get_loc(e)
            df.iloc[:,c-1]
            df[e] = (df.iloc[:,c-1]+0.85*df.iloc[:,c+1])/2
        # set NaN based on range
        else:
            tLower = err.get(e)[0]
            tUpper = err.get(e)[1]
            df.loc[(df['t'] >= tLower) & (df['t'] <= tUpper), e] *= np.nan

    # Save file
    df.to_csv(os.path.join(outputFilePath,(shotName + "_T.csv")),index = False),
    print("     - " + shotName + "_T.csv created" )

def processPress(inputFilePath, outputFilePath,shotName):
    """
    Pressure and thrust file conversion
    """
    # Get file information
    specs = json.load(open(os.path.join("file_spec",(shotName + ".json")), encoding="utf-8")).get("pressureThrust")
    fileName = os.path.join(inputFilePath,(shotName + "_P.csv"))

    # Load file data
    if not os.path.isfile(fileName):
        print("     WARNING: " + shotName + "_P.csv does not exist. skipping conversion...")
        return
    else:
        fileHeader = specs.get("header")
        df = pd.read_csv(fileName, names=fileHeader)

    # Get sensors and calibration
    sensorArray = specs.get("header")
    calibration = specs.get("calibration", None)

    # Convert Data
    for s in sensorArray:

        # Convert load cell voltage to thrust (N)
        if s == "T":
            v0Range = df[(df["t"] > -4) & (df["t"] < -2)]["T"]
            v0 =  v0Range.mean()
            cal = calibration[sensorArray.index(s)]
            df[s] = (df[s]-v0)*cal

        # Convert PT voltage data to presure (MPa)
        elif (s != "t") & (s != "trig") & (s != "T"):
            cal = calibration[sensorArray.index(s)]/10
            df[s] = df[s]*cal

    # Save file
    df.to_csv(os.path.join(outputFilePath,(shotName + "_P.csv")),index = False),
    print("     - " + shotName + "_P.csv created" )

def processVibe(inputFilePath, outputFilePath,shotName):
    """
    Vibration file conversion
    """
    # Get file information
    specs = json.load(open(os.path.join("file_spec",(shotName + ".json")), encoding="utf-8")).get("vibration")
    fileName = os.path.join(inputFilePath,(shotName + "_V.csv"))

    # Load file data
    if not os.path.isfile(fileName):
        print("     WARNING: " + shotName + "_V.csv does not exist. skipping conversion...")
        return
    else:
        fileHeader = specs.get("header")
        df = pd.read_csv(fileName, names=fileHeader)

    # Save file
    df.to_csv(os.path.join(outputFilePath,(shotName + "_V.csv")),index = False),
    print("     - " + shotName + "_V.csv created" )

if __name__ == "__main__":

    # Set conversion files
    settings = json.load(open("file_spec//runSettings.json", encoding="utf-8"))
    dir = settings.get("stage")
    tests = settings.get("testNames")
    process = settings.get("process")
     # Create data filepath
    inputParent = os.path.join(Path(__file__).parent.parent.parent, "data_raw", dir)
    outputParent = os.path.join(Path(__file__).parent.parent.parent, "data_processed", dir)

    # Iterate through files and process
    for t in tests:
        # Set test directories for input and output
        inputDir = str(os.path.join(inputParent,t))
        outputDir = str(os.path.join(outputParent,t))
        print("\n Processing:" + t +":")

        # Make directory if it doesnt already exist
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        # Process temperatures
        if "T" in process:
            processTemps(inputDir, outputDir,t)
        # Process pressures
        if "P" in process:
            processPress(inputDir, outputDir,t)
        # Process vibration
        if "V" in process:
            processVibe(inputDir, outputDir,t)

    print("Processing completed")


"""
 plt.style.use('science')
    plt.rcParams["figure.figsize"] = (6,3)
    plt.rcParams["axes.linewidth"] = 1.25
    plt.rcParams["xtick.major.width"] = 1.25
    plt.rcParams["xtick.minor.width"] = 1.25
    plt.rcParams["ytick.major.width"] = 1.25
    plt.rcParams["ytick.minor.width"] = 1.25
    fig, (ax1,ax2) = plt.subplots(1, 2)
    ax1.plot(df["t"],df["T1"])
    ax1.plot(df["t"],df["T2"])
    ax1.plot(df["t"],df["T3"])
    ax1.plot(df["t"],df["T4"])
    ax1.plot(df["t"],df["T5"])
    ax1.plot(df["t"],df["T6"])
    ax1.plot(df["t"],df["T7"])
    ax1.plot(df["t"],df["T8"])
    ax1.plot(df["t"],df["T9"])
    ax1.set_xlim(-0.1,1.5)
    ax1.set_ylim(15,240)
    ax1.legend([ "$T_1$", "$T_2$", "$T_3$","$T_4$","$T_5$","$T_6$","$T_7$","$T_8$","$T_9$"],loc="upper left")

    T1 = df.loc[df['t'] == 1.0].iloc[:,1:].T
    print(T1)
    z = [1, 4.5, 8, 11.5, 15, 25, 35, 45, 55]
    ax2.plot(z, T1[:-1],'#0C5DA5', marker = 'o',markersize='5')
    ax2.set_ylim(15,240)
    plt.savefig("tcomparison.pdf", format="pdf", bbox_inches="tight")



"""