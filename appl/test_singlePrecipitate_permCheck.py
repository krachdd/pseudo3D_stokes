#!/usr/bin/env python3
import pandas as pd
import sys
import subprocess

def read_and_compare_csv(csvFile_reference, columnOpenFoam, columnPoremaps, columnErrors, csvFileDumux, columnDumux, columnDomainName):
    # Read the CSV files
    df_reference = pd.read_csv(csvFile_reference)
    df_dumux = pd.read_csv(csvFileDumux)
    # Check if the specified columns exist in the dataframes
    if columnOpenFoam not in df_reference.columns:
        raise ValueError(f"Column {columnOpenFoam} not found in {csvFile_reference}")
    if columnPoremaps not in df_reference.columns:
        raise ValueError(f"Column {columnPoremaps} not found in {csvFile_reference}")
    if columnErrors not in df_reference.columns:
        raise ValueError(f"Column {columnErrors} not found in {csvFile_reference}")
    if columnDumux not in df_dumux.columns:
        raise ValueError(f"Column {columnDumux} not found in {csvFileDumux}")
    if columnDomainName not in df_dumux.columns:
        raise ValueError(f"Column {columnDomainName} not found in {csvFileDumux}")
    
    # Extract the columns to compare
    perm_openFOAM = df_reference[columnOpenFoam]
    perm_poremaps = df_reference[columnPoremaps]
    perm_errors = df_reference[columnErrors]
    perm_dumux = df_dumux[columnDumux]
    domainNames_dumux = df_dumux[columnDomainName]

    # calculate the absolute difference
    permDiff_abs = perm_openFOAM - perm_dumux
    permDiff_norm = abs(permDiff_abs / perm_openFOAM)*100
    deviationFromInitialErrors = permDiff_norm - perm_errors
    return permDiff_abs, permDiff_norm, deviationFromInitialErrors, domainNames_dumux

def checkErrors(deviationFromInitialErrors, domainNames_dumux, tolerance):
    exceeding_indices = []
    exceeding_values = []
    print("\n*********************************")
    print("    Test results:")
    print("*********************************")
    for index, value in enumerate(deviationFromInitialErrors):
        if abs(value) > tolerance:
            sys.stderr.write("Permeability test for " + domainNames_dumux[index] + " failed\n" + " the deviation is {:.2f}%".format(abs(value)) + " and exceeds the given tolerance of {:.2}%\n".format(tolerance))
            sys.exit(1)
        else:
            print("Permeability test for " + domainNames_dumux[index] + " passed")
    

# defining the paths to csv files and defining parameters
csvFile_reference       = 'references/permeabilities_single_precipitate-reference.csv'
columnOpenFoam          = ' k11 fully 3D OpenFOAM [m^2]'
columnPoremaps          = ' k11 fully 3D poremaps [m^2]'
columnErrors            = ' total error [%]'
csvFile_simulation      = 'test_singlePrecipitate_2d_total/test_results.csv'
columnDumux             = ' k11 [m^2]'
columnDomainName_dumux  = 'file'
tolerance = 1.0 # [%]

# start the simulation of the single precipitates with various sizes
script_to_run = "run_simulations_test.sh"
print(f'Starting the simulations with bash {script_to_run} .')
command = ["bash", script_to_run]
process = subprocess.Popen( command , stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
for line in process.stdout:
    print(line, end='')

# Wait for the process to complete and get the return code
return_code = process.wait()

# read the results and calculate the differences to the reference solution
permDiff_abs, permDiff_norm, deviationFromInitialErrors, domainNames_dumux = read_and_compare_csv(csvFile_reference, columnOpenFoam, columnPoremaps, columnErrors, csvFile_simulation, columnDumux, columnDomainName_dumux)

# check if the errors are within the defined tolerance
checkErrors(deviationFromInitialErrors, domainNames_dumux, tolerance)
