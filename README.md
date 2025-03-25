# Virus Neutralization Data Analysis

## Access the Web Application

You can access the **Virus Neutralization Assay Data Analysis** web application through the following link:

[**Virus Neutralization Assay Data Analysis Application**](https://mohamedmire-virus-neutralization-assay-data-analysis.share.connect.posit.cloud)

## Overview

This application processes and analyzes data from high throughput luciferase-based virus neutralization assays. The luminescence plate reader generates data in Excel format, with results across three to five sheets (depending on instrument data output settings). The middle three sheets contain the assay data in various formats for further processing. Examble data file formats can be found under the Data folder; these can be used as test runs for the application.

## Excel File Structure

The Excel file contains the following sheets:

- **Firefly luc**: Contains assay run details and machine information.
- **Results**: Assay data in a 96-well plate format (less convenient for analysis).
- **Results By Well**: Common format used for data processing.
- **Results Table**: A structured format ideal for generating dictionaries.
- **Results Charts**: A blank 96-well plate template (not used by the application).

## Data Normalization

To standardize the data, we apply the following normalization formula:

**Normalization Formula**:  
`x' = ((x - MIN) / (MAX - MIN)) * 100`

Where:
- **x'** = normalized value (as a percentage)
- **x** = original value
- **MIN** = minimum value or baseline (0% neutralization)
- **MAX** = maximum value or defined maximum (100% neutralization)

## IC50 Calculations

IC50 values are calculated using the following dose-response function based on the GraphPad Prism 10 Curve Fitting Guide:

`100 / (1 + 10**((logIC50 - x) * HillSlope))`

We use the **lmfit** Python module to fit the curve to the normalized data. IC50 values are adjusted based on control sera calculations, with a minimum threshold of 10 for any IC50 values below that.

## Application Usage Guide

### 1. Data Upload
Begin by uploading the necessary files under the **'Analyze Data'** tab:
- **Sample Names (TXT file)**: Contains columns for `sample_number` and `Sample`.
- **Sample Dilution Factors (TXT file)**: Contains a column for `dilution_factors`.
- **Control Dilution Factors (TXT file)**: Contains a column for `dilution_factors`.
- **Excel Files (raw data)**: You can upload multiple Excel files (e.g., `RSV-A2-1.xlsx`).

### 2. Processing Options
- Select the method for deriving the reference value and specify the variant of interest.
- **Important**: Ensure the variant input matches the file name in both spelling and case (e.g., file name `RSV-A2-1` should match variant input `A2`). This field is pre-populated with common SARS-CoV-2 and RSV strains for convenience. 

### 3. Data Processing
Click **'Process Data'** to begin the analysis. The application will:
- Read and normalize the uploaded data.
- Calculate the IC50 values.
- Generate plots for visualization.

### 4. Results
After processing, results will be available for each plate and sample type:
- **Raw Data**
- **Normalized Data**
- **IC50 Plots**
- **Results Table**
  
**NOTE:** With the exception of the Adj_factor, all columns are rounded to four decimal places.

### 5. Download
Once processing is complete, you can download:
- **Aggregated results**
- **All generated plots**

A confirmation message will appear once your downloads are complete. 
