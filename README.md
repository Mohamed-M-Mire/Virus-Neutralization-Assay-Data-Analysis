# Virus-Neutralization-Data-Analysis
## Description

This application is designed to process and analyze data from luciferase-based virus neutralization assays.

The lab makes use of the Promega Luciferase Assay System. The plate reader provides data in excel form in a specified format.

Within the excel file, there are five sheets; the middle three sheets have the assay results in three different formats.

## Excel File Structure

###   Firefly luc - Contains assay run and machine information

###   Results - Assay data in 96 well plate format (difficult to work with)

###   Results By Well - Common format used for data processing

###   Results Table - Data presented in a manner conducive to creating dictionaries

### Results Charts - Empty 96 well plate (not used in this application)

# Data Normalization

Normalization Formula: x' = ((x - MIN) / (MAX - MIN)) * 100

Where:

x' = normalized value (as a percentage)

x = original value

MIN = minimum value or defined baseline (0% neutralization)

MAX = maximum value or defined maximum (100% neutralization)

# IC50 Calculations

The dose-response function is based on the GraphPad Prism 10 Curve Fitting Guide:

100 / (1 + 10**((logIC50 - x) * HillSlope))

We use the lmfit python module to derive the best fit line for our normalized data.

IC50 values are adjusted based on the adjustment factors derived from the respective positive control sera calculations.

Samples with IC50 values < 10 are adjusted to 10.

# Application Usage Guide

  1. Data Upload

    Start by uploading your data files in the 'Analyze Data' tab:
    
    Sample Names (TXT file) - with columns: sample_number and Sample
    
    Sample Dilution Factors (TXT file) - with column: dilution_factors
    
    Control Dilution Factors (TXT file) - with column: dilution_factors
    
    Excel Files containing raw data (multiple files allowed) - Example: RSV-A2-1 or A2-1

  2. Processing Options

    Select how to derive the reference value and specify the variants.
    
    NOTE: The variant input needs to match in both spelling and font case; file name is RSV-A2-1 or A2-1, then variant input for it must be A2.

  3. Data Processing

    Click 'Process Data' to start the analysis. The application will:
    
        Read and normalize the data
        
        Calculate IC50 values
        
        Generate plots
  
  4. Results

    After processing, you can view, per plate # and sample type:

      Raw Data
      
      Normalized Data
      
      IC50 Plots
      
      Results Table

  5. Download

    You can download the aggregated results and all generated plots.
    
    Once downloads are done you will get a confirmation message.
    
    When analyzing a large number of plates, be aware this process might take a while.
