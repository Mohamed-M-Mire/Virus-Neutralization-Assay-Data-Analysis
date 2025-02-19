from typing import List, Dict, Tuple, Literal
from shiny import App, reactive, render, ui
from lmfit import Model, Parameters
import matplotlib.pyplot as plt
from zipfile import ZipFile
import pandas as pd
import numpy as np
import base64
import io

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Defining the functions needed to perform data manipulation and analysis  
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Compile assay metaData 
def render_metadata(
    sample_filename, 
    sample_dilution_factors_filename, 
    ctrl_dilution_factors_filename
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Processes assay metadata information from input files.

    This function reads sample names, sample dilution factors, and control sera dilution factors
    from separate files. It duplicates the sample names and reassigns sample numbers.
    It also validates that dilution factor files have the correct format.

    Args:
        sample_filename (str): Path to the text file containing sample IDs.
        sample_dilution_factors_filename (str): Path to the file with sample dilution factors.
        ctrl_dilution_factors_filename (str): Path to the file with control dilution factors.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: A tuple containing:
            - sample_names_duplicated: Duplicated sample names with reassigned sample numbers.
            - sample_dilution_factors: Sample dilution factors.
            - ctrl_dilution_factors: Control dilution factors.

    Raises:
        FileNotFoundError: If any of the input files cannot be found.
        pd.errors.EmptyDataError: If any of the input files is empty.
        ValueError: If the data in the files is not in the expected format.
    """
    try:
        sample_names= pd.read_csv(sample_filename, sep="\t", index_col="sample_number")
        sample_names_duplicated= sample_names.loc[sample_names.index.repeat(2)].reset_index()
        sample_names_duplicated['sample_number']= range(1, len(sample_names_duplicated) + 1)
        sample_names_duplicated.set_index('sample_number', inplace=True)

        sample_dilution_factors= pd.read_csv(sample_dilution_factors_filename, sep="\t")
        ctrl_dilution_factors= pd.read_csv(ctrl_dilution_factors_filename, sep="\t")

        return sample_names_duplicated, sample_dilution_factors, ctrl_dilution_factors
    except Exception as e:
        raise ValueError(f"Error processing input files: {str(e)}")

# Compile assay data
def read_raw_data(
    excel_path_name_dict: Dict[str, str],
    sample_dilutionFactors: pd.DataFrame,
    sample_names_duplicated: pd.DataFrame,
    ctrl_dilutionFactors: pd.DataFrame
    ) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    
    """
    Read and process raw data from luciferase-based virus neutralization assays.

    This function reads multiple Excel files, extracts luminescence data, and organizes it into
    structured DataFrames for further manipulations.

    Args:
        excel_path_name_dict (Dict[str, str]): A dictionary mapping file names to file paths of Excel files
            containing virus plate data. Each file should have a 'Results By Well' sheet with luminescence data.
        sample_dilutionFactors (pd.DataFrame): DataFrame containing dilution factors for samples.
        sample_names_duplicated (pd.DataFrame): DataFrame containing duplicated sample names.
        ctrl_dilutionFactors (pd.DataFrame): DataFrame containing dilution factors for control samples.

    Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: A tuple containing two dictionaries:
        - virus_plate_rawData: A dictionary where:
            - Keys (str): Plate names (derived from filenames without extension)
            - Values (pd.DataFrame): DataFrames with the following structure:
                - Index: Dilution factors (from sample_dilutionFactors["dilution_factors"])
                - Columns: Sample names and two additional columns:
                    - 'virus_cells_only': Data for virus cells only
                    - 'positive_ctrl': Data for positive control
        - postv_ctrl_rawData: A dictionary where:
            - Keys (str): Plate names (same as in virus_plate_rawData)
            - Values (pd.DataFrame): DataFrames with the following structure:
                - Index: Dilution factors (from ctrl_dilutionFactors["dilution_factors"])
                - Columns: Two columns, both named "Positive Control Sera"

    Raises:
        Exception: If there's an error processing any of the input files.

    Note:
        This function processes each Excel file to extract luminescence data, assigns appropriate
        column names based on sample information, and organizes the data into separate DataFrames
        for virus plate data and positive control data.

    Example:
    virus_plate_rawData, postv_ctrl_rawData = read_raw_data(
        excel_path_name_dict, sample_dilutionFactors, sample_names_duplicated, ctrl_dilutionFactors
        )
    """
    virus_plate_rawData= {}
    postv_ctrl_rawData= {}

    for filename, filepath in excel_path_name_dict.items():
        try:
            df= pd.read_excel(filepath, sheet_name="Results By Well", skiprows=5)
            luminescence_indices= df.index[df['Unnamed: 1'] == 'Luminescence'].tolist()
            values= [df.loc[idx, 'Unnamed: 6'] for idx in luminescence_indices]
            result_df= pd.DataFrame(np.array(values).reshape(8, 12))
            result_df.rename(columns={result_df.columns[10]: 'virus_cells_only', result_df.columns[11]: 'positive_ctrl'}, inplace=True)
            result_df.index= sample_dilutionFactors["dilution_factors"]
            virus_plate_rawData[filename]= result_df

            if filename.endswith("1"):
                new_column_names= sample_names_duplicated.loc[:10].values.flatten().tolist()
            elif filename.endswith("2"):
                new_column_names= sample_names_duplicated.loc[11:20].values.flatten().tolist()
            elif filename.endswith("3"):
                new_column_names= sample_names_duplicated.loc[21:30].values.flatten().tolist()
            elif filename.endswith("4"):
                new_column_names= sample_names_duplicated.loc[31:40].values.flatten().tolist()
            elif filename.endswith("5"):
                new_column_names= sample_names_duplicated.loc[41:50].values.flatten().tolist()
            
            virus_plate_rawData[filename].columns= new_column_names + list(virus_plate_rawData[filename].columns[10:])

            odd_vals= virus_plate_rawData[filename].positive_ctrl.values[1::2]
            even_vals= virus_plate_rawData[filename].positive_ctrl.values[::2]
            postv_ctrlData= pd.DataFrame(data=(even_vals, odd_vals)).T
            postv_ctrlData.columns= ["Positive Control Sera"]*2
            postv_ctrlData.index= ctrl_dilutionFactors["dilution_factors"]
            postv_ctrl_rawData[filename]= postv_ctrlData
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")

    return virus_plate_rawData, postv_ctrl_rawData

# Implement data normalization 
def normalize_data(
    raw_data: Dict[str, pd.DataFrame],
    postv_ctrl_rawData: Dict[str, pd.DataFrame],
    derive_reference_value_by: Literal["The last row", "The virus only wells"]
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Normalize sample and positive control data from luciferase-based virus neutralization assays.

    The normalization method based on:
        x' = ((x - MIN) / (MAX - MIN)) * 100
        Where:
            x'     = normalized value (as a percentage)
            x      = original value
            MIN    = minimum value or defined baseline (0% neutralization)
            MAX    = maximum value or defined maximum (100% neutralization)
            
    Args:
        raw_data (Dict[str, pd.DataFrame]): Dictionary containing raw assay data DataFrames.
            Keys are plate names, values are DataFrames with raw luminescence data.
        postv_ctrl_rawData (Dict[str, pd.DataFrame]): Dictionary containing raw positive control data.
            Keys are plate names, values are DataFrames with raw positive control luminescence data.
        derive_reference_value_by (Literal["row", "virus"]): Method to derive the reference low value.
            "row": Uses the average of the highest three values in the last row of the raw data.
            "virus": Uses the average of virus-only wells.

    Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: A tuple containing two dictionaries:
            - normalized_data: Normalized sample data.
                Keys are plate names, values are DataFrames with normalized data (0-100 scale).
            - normalized_postv_CTRL_data: Normalized positive control data.
                Keys are plate names, values are DataFrames with normalized data (0-100 scale).

    Raises:
        ValueError: If an invalid value is provided for derive_reference_value_by.

    Notes:
        - The normalization formula used is: ((x - reference_low) / (reference_high - reference_low)) * 100
        - reference_high is always the average of cells-only wells
        - reference_low depends on the derive_reference_value_by parameter for both sample and positive control data

    Example:
        raw_data = {...}  # Dictionary of raw assay DataFrames
        postv_ctrl_rawData = {...}  # Dictionary of raw positive control DataFrames
        normalized_data, normalized_postv_CTRL_data = normalize_data(raw_data, postv_ctrl_rawData, "row")
    """
    normalized_data= {}
    normalized_postv_CTRL_data= {}

    def normalize(x, reference_low, reference_high):
        return ((x - reference_low) / (reference_high - reference_low)) * 100

    for key, df in raw_data.items():
        cells_only_wells= df.iloc[-2:, 10]
        average_cells= np.mean(cells_only_wells)
        
        virus_only_wells= df.iloc[:6, 10]
        ave_virusWells= np.mean(virus_only_wells)

        last_row= df.iloc[-1, 0:10]
        highest_three_ave= np.mean(np.sort(last_row.values)[-3:])
        
        if derive_reference_value_by == "The last row":
            reference_low_sample= highest_three_ave
            reference_low_ctrl= np.mean(postv_ctrl_rawData[key].iloc[-1])
        elif derive_reference_value_by == "The virus only wells":
            reference_low_sample= ave_virusWells
            reference_low_ctrl= ave_virusWells
        else:
            raise ValueError("Invalid value for derive_reference_value_by. Must be 'last row' or 'virus'.")

        # Normalize sample data
        norm_data= df.iloc[:, 0:10].apply(lambda x: normalize(x, reference_low_sample, average_cells))
        normalized_data[key]= norm_data

        # Normalize positive control data
        norm_posCtrl= postv_ctrl_rawData[key].apply(lambda x: normalize(x, reference_low_ctrl, average_cells))
        normalized_postv_CTRL_data[key]= norm_posCtrl

    return normalized_data, normalized_postv_CTRL_data

# Dose response function 
def dose_response(x, logIC50, HillSlope):
    return 100 / (1 + 10**((logIC50 - x) * HillSlope))

# IC50 calculation helper function 1 0f 2
def estimate_initial_logIC50(x_data, y_data):
    idx= np.argmin(np.abs(y_data - 50))
    initial_logIC50= x_data[idx]
    print(f"Initial logIC50 estimate: {initial_logIC50}")
    return initial_logIC50

# IC50 calculation helper function 2 0f 2
def calculate_ic50(x_data, y_data):
    try:
                
        if len(x_data) != len(y_data) or len(x_data) == 0:
            raise ValueError("Invalid input data: x_data and y_data must have the same non-zero length")
        
        model= Model(dose_response)
        
        params= Parameters()
        initial_logIC50= estimate_initial_logIC50(x_data, y_data)
        params.add('logIC50', value=initial_logIC50)
        params.add('HillSlope', value=-1)
        
        print(f"Initial parameters: logIC50 = {initial_logIC50}, HillSlope = -1")
        
        result= model.fit(y_data, params, x=x_data)
        
        logIC50_fit= result.params['logIC50'].value
        HillSlope_fit= result.params['HillSlope'].value
        
        print(f"Fitted parameters: logIC50 = {logIC50_fit}, HillSlope = {HillSlope_fit}")
        print(f"Fit success: {result.success}")
        print(f"Fit message: {result.message}")
        
        IC50= 10**logIC50_fit
        print(f"Calculated IC50: {IC50}")
        
        r_squared= result.summary()["rsquared"]
        print(f"R-squared: {r_squared}")
        
        return IC50, logIC50_fit, HillSlope_fit, r_squared
    except ValueError as ve:
        print(f"\033[1;31mValue Error in calculate_ic50: {ve}\033[0m")
        return np.nan, np.nan, np.nan, np.nan
    except Exception as e:
        print(f"\033[1;31mError in calculate_ic50: {e}\033[0m")
        return np.nan, np.nan, np.nan, np.nan

# Curve plot helper function 1 of 2
def plot_dose_response_curve(ax, x_data, y_data, sample_name, ic50, logIC50, HillSlope, log_dilution_factors):
    ax.semilogx(x_data, y_data, 'o', label= f'{sample_name} (Data)')
    if not np.isnan(ic50):
        x_fit= np.logspace(min(log_dilution_factors), max(log_dilution_factors), 100)
        y_fit= dose_response(np.log10(x_fit), logIC50, HillSlope)
        ax.semilogx(x_fit, y_fit, '-', label= f'{sample_name} (Fit)')
        print(f"\033[1;32mCurve fitted for {sample_name}\033[0m")
    else:
        print(f"\033[1;31mFailed to fit curve for {sample_name}\033[0m")

# Curve plot helper function 2 of 2
def format_plot(ax, key):
    ax.set_xlabel('Log(Dilution Factor)', fontsize= 14)
    ax.set_ylabel('Response (%)', fontsize= 14)
    ax.axhline(y=50, color="gray", linestyle= '--', linewidth= 2, alpha= 0.6)
    ax.text(ax.get_xlim()[0] + (ax.get_xlim()[0]/100), 51, '50% Neutralization', color= 'black', fontsize= 10, fontweight= 'bold')
    ax.set_title(f'Dose-Response Curve - {key[:-2]}', weight= 'bold', fontsize= 16)
    ax.legend(loc= 'center left', bbox_to_anchor= (1.02, 0.5))
    ax.grid(True)
    plt.tight_layout()

# Derive IC50 values and plot results    
def process_and_plot_results(
    Norm_data: dict[str, pd.DataFrame],  
    sample_type: Literal["ctrl", "sample"]
    ) -> tuple[dict[str, pd.DataFrame], dict[str, plt.Figure]]:
    """
    Process normalized dose-response data, calculate IC50 values, and generate plots.

    This function analyzes normalized data for multiple plates, calculates IC50 values
    for each sample, creates dose-response curves, and saves the plots as image files.

    Args:
        Norm_data (dict[str, pd.DataFrame]): Normalized data for each plate.
            Keys are plate identifiers, values are DataFrames with normalized data.
        file_path (str): Directory path to save generated plots.
        sample_type (Literal["Ctrl", "sample"]): Type of sample being processed.

    Returns:
        tuple[dict[str, pd.DataFrame], dict[str, plt.Figure]]: A tuple containing:
            1. results_df: Dictionary of DataFrames with calculated values for each plate.
                Keys are plate identifiers, values are DataFrames with IC50, logIC50,
                    Hill Slope, and R-squared values for each sample.
            2. figure_dict: Dictionary of matplotlib Figure objects for each plate.
                Keys are plate identifiers, values are the corresponding Figure objects.

    Process for each plate:
    1. Calculate IC50, logIC50, Hill Slope, and R-squared values for each sample.
    2. Generate a semi-log plot of the dose-response curve for all samples.
    3. Save the plot as a PNG file in the appropriate subdirectory based on sample_type.
    4. Compile results into a DataFrame and store the Figure object.

    Note:
        - Utilizes `calculate_ic50` and `dose_response` functions for calculations.
        - Prints progress information and any errors during processing.
        - Plots include both raw data points and fitted curves (when possible).
        - Helper functions `plot_sample_data`, `format_plot`, and `save_plot` are used
            for plotting and figure management.
        - Saves plots in 'sample_results' or 'ctrl_results' subdirectory based on sample_type.
        - Creates the necessary subdirectory if it doesn't exist.

    Example usage:
        postv_ctrl_IC50_results, postv_ctrl_Figures = process_and_plot_results(
            normalized_postv_CTRL_data, 
            file_path= "./",
            sample_type= "ctrl"
            )
            
        or
        postv_ctrl_IC50_results = process_and_plot_results(
            normalized_postv_CTRL_data, 
            file_path= "./",
            sample_type= "ctrl" 
            )
        and postv_ctrl_IC50_results[0] holds dict of results df 
        and postv_ctrl_IC50_results[1] holds dict of plt figures
    """
    results_df= {}
    figure_dict= {}
    
    for key, norm_DF in Norm_data.items():
        print(f"\n\033[1;33mProcessing plate: {key}\033[0m")
        
        log_dilution_factors= np.log10(norm_DF.index)        
        results_dict= {}
        
        fig, ax= plt.subplots(figsize=(9, 6))
        
        for i in range(0, len(norm_DF.columns), 2):
            sample_name= norm_DF.columns[i]
            print(f"\nStarting\033[1;33m sample {sample_name}\033[0m IC50 calculation...")
            
            sample_data= norm_DF.iloc[:, i:i+2].mean(axis=1)
            
            ic50, logIC50, HillSlope, r_squared = calculate_ic50(log_dilution_factors, sample_data)
            
            results_dict[sample_name]= {'IC50': ic50, 'logIC50': logIC50, 'HillSlope': HillSlope, "r_squared": r_squared}
            
            plot_dose_response_curve(ax, norm_DF.index, sample_data, sample_name, ic50, logIC50, HillSlope, log_dilution_factors)

        results_df[key]= pd.DataFrame(results_dict)
        format_plot(ax, key)
        
        # Save the plot as a base64 encoded string...
        # This only way I found to avoid problems when user goes back and forth in selecting the same images
        buffer= io.BytesIO()
        fig.savefig(buffer, format='png', dpi= 300)
        buffer.seek(0)
        image_base64= base64.b64encode(buffer.getvalue()).decode('utf-8')
        figure_dict[key]= image_base64
        
        plt.close(fig)
    
    return results_df, figure_dict

# Adjust sample IC50 values in relation to control 
def adjust_ic50_values(
    postv_ctrl_IC50_results: Dict[str, pd.DataFrame],
    sample_IC50_results: Dict[str, pd.DataFrame],
    variants: List[str]
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """
    Adjust IC50 values based on positive control sera data and calculate adjusted sample IC50 values.

    Args:
        postv_ctrl_IC50_results (Dict[str, pd.DataFrame]): Dictionary of positive control sera IC50 results.
        sample_IC50_results (Dict[str, pd.DataFrame]): Dictionary of sample IC50 results.
        variants (List[str]): List of variant names to consider. These should match the characters
                            in the xls_files and have the same case.

    Returns:
        Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]: A tuple containing:
            1. Dictionary of positive control IC50 results with the calculated adjustment factors.
            2. Dictionary of adjusted sample IC50 results.

    Raises:
        KeyError: If a key is missing in postv_ctrl_IC50_results_w_adj_factor.
        Exception: For any other error during processing.
    
    Example usage:
        variants = ['UK', 'A2', 'B1', 'WA', 'KP', 'BA']
        postv_ctrl_IC50_results_w_adj_factor, Adjusted_sample_IC50_results= adjust_ic50_values(postv_ctrl_IC50_results, sample_IC50_results, variants)
    """
    # Dictionary to store average values for each variant
    ave_values= {variant: [] for variant in variants}

    for key, df_data in postv_ctrl_IC50_results.items():
        for variant in ave_values:
            if variant in key:
                ave_values[variant].append(df_data.iloc[0, :].values)
                break

    postv_ctrl_IC50_results_w_adj_factor= {}

    for key, df_data in postv_ctrl_IC50_results.items():
        for variant in ave_values:
            if variant in key:
                adj_factor= pd.DataFrame((np.mean(ave_values[variant])/df_data.iloc[0, :]).rename("Adj_factor")).T
                postv_ctrl_IC50_results_w_adj_factor[key]= pd.concat([df_data, adj_factor])
                break

    Adjusted_sample_IC50_results= {}

    for key, sample in sample_IC50_results.items():
        try:
            postv_ctrl= postv_ctrl_IC50_results_w_adj_factor[key]
            
            adj_ic50= pd.DataFrame(
                sample.iloc[0, :].values * postv_ctrl.iloc[-1, :].values, 
                columns=["Adjusted_IC50"]
            ).clip(lower=10).T  # Adjusted IC50 values are set to 10 if lower than 10
            
            adj_ic50.columns= sample.columns
            Adjusted_sample_IC50_results[key]= pd.concat([sample, adj_ic50])
            
        except KeyError:
            print(f"\033[1;31mError: Missing key {key} in postv_ctrl_IC50_results_w_adj_factor\033[0m")
        except Exception as e:
            print(f"\033[1;31mError processing key {key}: {str(e)}\033[0m")

    return postv_ctrl_IC50_results_w_adj_factor, Adjusted_sample_IC50_results

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Implementation of the user interface logic 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Function to setup the 'About' page 
def create_info_tab():
    return ui.card(
        ui.card_header(
            "Neutralization Assay Analysis",
            style="background-color: #00274C; color: #FFCB05; font-size: 25px; font-weight: bold; padding: 10px;"
        ),
        ui.h4(ui.strong("Description")),
        ui.div(
            ui.p("This application is designed to process and analyze data from luciferase-based virus neutralization assays.", style= "margin-top: -25px;"),
            ui.p("The lab makes use of the Promega Luciferase Assay System. The plate reader provides data in excel form in a specified format.", style= "margin-top: -20px;"),
            ui.p("Within the excel file, there are five sheets; the middle three sheets have the assay results in three different formats.", style= "margin-top: -20px;"),       
        ),
        ui.h4(ui.strong("Excel File Structure"), style= "margin-top: -15px;"),
        ui.div(
            ui.tags.ul(
                ui.tags.li(ui.strong("Firefly luc"), " - Contains assay run and machine information", style= "margin-top: -20px;"),
                ui.tags.li(ui.strong("Results"), " - Assay data in 96 well plate format (difficult to work with)"),
                ui.tags.li(ui.strong("Results By Well"), " - Common format used for data processing"),
                ui.tags.li(ui.strong("Results Table"), " - Data presented in a manner conducive to creating dictionaries"),
                ui.tags.li(ui.strong("Results Charts"), " - Empty 96 well plate (not used in this application)")
            )
        ),
        ui.h4(ui.strong("Data Normalization"), style= "margin-top: -25px;"),
        ui.div(
            ui.p("Normalization Formula: x' = ((x - MIN) / (MAX - MIN)) * 100", style= "margin-top: -25px;"),
            ui.p("Where:", style= "margin-top: -10px;"),
            ui.tags.ul(
                ui.tags.li("x' = normalized value (as a percentage)", style= "margin-top: -10px;"),
                ui.tags.li("x = original value"),
                ui.tags.li("MIN = minimum value or defined baseline (0% neutralization)"),
                ui.tags.li("MAX = maximum value or defined maximum (100% neutralization)")
            )
        ),
        ui.h4(ui.strong("IC50 Calculations"), style= "margin-top: -25px;"),
        ui.div(
            ui.p("The dose-response function is based on the GraphPad Prism 10 Curve Fitting Guide:", style= "margin-top: -20px;"),
            ui.p("100 / (1 + 10**((logIC50 - x) * HillSlope))"),
            ui.p("We use the lmfit python module to derive the best fit line for our normalized data."),
            ui.p("IC50 values are adjusted based on the adjustment factors derived from the respective positive control sera calculations.", style= "margin-top: -20px;"),
            ui.p("Samples with IC50 values < 10 are adjusted to 10.", style="margin-top: -20px;"),
        ),
        ui.h4(ui.strong("Application Usage Guide"), style= "margin-top: -15px;"),
        ui.div(
            ui.h5(ui.strong("1. Data Upload"), style= "color: #696969; margin-top: -20px;"),
            ui.p("Start by uploading your data files in the 'Analyze Data' tab:", style= "margin-bottom: -5px;"),
            ui.tags.ul(
                ui.tags.li("Sample Names (TXT file) - with columns: sample_number and Sample"),
                ui.tags.li("Sample Dilution Factors (TXT file) - with column: dilution_factors"),
                ui.tags.li("Control Dilution Factors (TXT file) - with column: dilution_factors"),
                ui.tags.li("Excel Files containing raw data (multiple files allowed) - Example: RSV-A2-1 or A2-1")
            ),
            ui.h5(ui.strong("2. Processing Options"), style= "color: #696969;"),
            ui.p("Select how to derive the reference value and specify the variants."),
            ui.p("NOTE: The variant input needs to match in both spelling and font case; file name is RSV-A2-1 or A2-1, then variant input for it must be A2."),
            ui.div(
                ui.h5(ui.strong("3. Data Processing"), style= "color: #696969;"),
                ui.p("Click 'Process Data' to start the analysis. The application will:", style= "margin-bottom: -1px;"),
                ui.tags.ul(
                    ui.tags.li("Read and normalize the data"),
                    ui.tags.li("Calculate IC50 values"),
                    ui.tags.li("Generate plots")
                )
            ),
            ui.div(
                ui.h5(ui.strong("4. Results"), style= "color: #696969;"),
                ui.p("After processing, you can view, per plate # and sample type:", style= "margin-bottom: -1px;"),
                ui.tags.ul(
                    ui.tags.li("Raw Data"),
                    ui.tags.li("Normalized Data"),
                    ui.tags.li("IC50 Plots"),
                    ui.tags.li("Results Table")
                ),
            ),
            ui.div(
                ui.h5(ui.strong("5. Download"), style= "color: #696969;"),
                ui.p("You can download the aggregated results and all generated plots.", style="margin-top: 0px;"),
                ui.p("Once downloads are done you will get a confirmation message.", style="margin-top: -20px;"),
                ui.p("When analyzing a large number of plates, be aware this process might take a while.", style="margin-top: -20px;"),
            ),
        ),
        ui.h4(ui.strong("Created By: Mohamed M. Mire, MS, PhD"), style="background-color: #00274C; color: #FFCB05;")
    )

app_ui= ui.page_navbar(
    ui.nav_panel("Application Guide", create_info_tab()),
    ui.nav_panel("Analyze Data",
        ui.layout_sidebar(
            ui.sidebar(
                ui.div(
                    ui.h3("Upload Data Files", style="margin-bottom: 10px;"),
                    ui.div(
                        ui.tags.style("""
                            .shiny-input-container {
                                margin-bottom: 5px !important;
                            }
                            .selectize-control {
                                margin-bottom: 5px !important;
                            }
                        """),
                        ui.input_file("sample_names", "Sample Names (TXT)", accept=[".txt"], button_label="Browse...", multiple=False),
                        ui.input_file("sample_dilution", "Sample Dilution Factors (TXT)", accept=[".txt"], button_label="Browse...", multiple=False),
                        ui.input_file("ctrl_dilution", "Control Dilution Factors (TXT)", accept=[".txt"], button_label="Browse...", multiple=False),
                        ui.input_file("excel_files", "Excel Files", accept=[".xlsx"], button_label="Browse...", multiple=True),
                        ui.input_select("derive_reference", "Derive 0% Reference from", choices=["The last row", "The virus only wells"]),
                        ui.input_text("variants", "Variants (comma-separated)", value="UK,A2,B1,WA,KP,BA"),
                        class_="shiny-input-container",
                    ),
                    ui.div(
                        ui.input_action_button("process", "Process Data", class_="btn-primary"),
                        ui.output_ui("process_status_message"),
                        ui.download_button("download_results", "Download IC50 Results", class_="btn-success"),
                        ui.output_ui("download_status_message"),
                        ui.download_button("download_all_plots", "Download All Plots", class_="btn-success"),
                        ui.output_ui("plot_download_status_message"),
                        style="margin-top: 5px; display: flex; flex-direction: column; gap: 2px;"
                    ),
                    style="display: flex; flex-direction: column; gap: 0px;"
                ),
                style="padding: 5px;",
                width=300
            ),
            ui.navset_tab(
                ui.nav_panel("Raw Data",
                    ui.card(
                        ui.card_header("Raw Data"),
                        ui.input_select("raw_data_key", "Select Plate", choices=[]),
                        ui.output_table("raw_data_table")
                    )
                ),
                ui.nav_panel("Normalized Data",
                    ui.card(
                        ui.card_header("Normalized Data"),
                        ui.input_select("normalized_data_key", "Select Plate", choices=[]),
                        ui.output_table("normalized_data_table")
                    )
                ),
                ui.nav_panel("IC50 Plots",
                    ui.card(
                        ui.card_header("IC50 Plots"),
                        ui.layout_sidebar(
                            ui.sidebar(
                                ui.input_select("plot_type", "Select Plot Type", choices=["Control", "Sample"]),
                                ui.input_select("plot_key", "Select Plate", choices=[]),
                                width="3"
                            ),
                            ui.output_ui("ic50_plot")
                        )
                    )
                ),
                ui.nav_panel("Results",
                    ui.card(
                        ui.card_header("Results"),
                        ui.input_select("results_type", "Select Results Type", choices=["Control", "Sample"]),
                        ui.output_data_frame("results_table")
                    )
                ),
                selected="IC50 Plots"
            )
        )
    ),
    title=ui.h1("Neutralization Assay Analysis", style="font-size: 35px; font-weight: bold;"),
    id="analysis_tabs",
    selected="Application Guide"
)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Implementation of the server logic 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

def server(input, output, Session):
    rv_processed_data= reactive.Value(None)
    rv_processing= reactive.Value(False)
    rv_downloading_results= reactive.Value(False)
    rv_downloading_plots= reactive.Value(False)

    @reactive.Effect
    @reactive.event(input.process)
    def process_data():
        rv_processing.set(True)
        try:
            
            # Reading input files
            sample_names_file= input.sample_names()
            sample_dilution_file= input.sample_dilution()
            ctrl_dilution_file= input.ctrl_dilution()
            sample_names_duplicated, sample_dilutionFactors, ctrl_dilutionFactors = render_metadata(
                sample_names_file[0]["datapath"],
                sample_dilution_file[0]["datapath"],
                ctrl_dilution_file[0]["datapath"]
            )

            # Reading Excel files
            excel_file_paths= [file['datapath'] for file in input.excel_files()]
            excel_file_names= [file['name'][:-5] for file in input.excel_files()]
            excel_path_name_dict= dict(zip(excel_file_names, excel_file_paths))
            virus_plate_rawData, postv_ctrl_rawData= read_raw_data(
                excel_path_name_dict,
                sample_dilutionFactors,
                sample_names_duplicated,
                ctrl_dilutionFactors
            )

            # Normalizing data
            normalized_data, normalized_postv_CTRL_data= normalize_data(
                virus_plate_rawData,
                postv_ctrl_rawData,
                input.derive_reference()
            )

            # Processing and plotting results
            postv_ctrl_IC50_results, postv_ctrl_Figures= process_and_plot_results(
                normalized_postv_CTRL_data,
                sample_type="ctrl"
            )
            sample_IC50_results, sample_IC50_Figures= process_and_plot_results(
                normalized_data,
                sample_type="sample"
            )

            # Adjusting IC50 values
            variants= [v.strip() for v in input.variants().split(',')]
            postv_ctrl_IC50_results_w_adj_factor, Adjusted_sample_IC50_results = adjust_ic50_values(
                postv_ctrl_IC50_results,
                sample_IC50_results,
                variants
            )
            aggregate_control_data= pd.concat(postv_ctrl_IC50_results_w_adj_factor, axis=1).T
            aggregate_sample_data= pd.concat(Adjusted_sample_IC50_results, axis=1).T

            rv_processed_data.set({
                "virus_plate_rawData": virus_plate_rawData,
                "normalized_data": normalized_data,
                "postv_ctrl_Figures": postv_ctrl_Figures,
                "sample_IC50_Figures": sample_IC50_Figures,
                "aggregate_control_data": aggregate_control_data,
                "aggregate_sample_data": aggregate_sample_data
            })

            ui.notification_show("Data processing complete!", type="message")
        except Exception as e:
            ui.notification_show(f"Error during data processing: {str(e)}", type="error")
        finally:
            rv_processing.set(False)

    @reactive.Effect
    @reactive.event(rv_processed_data)
    def update_dropdown_choices():
        if rv_processed_data() is not None:
            plate_keys= list(rv_processed_data()["virus_plate_rawData"].keys())
            ui.update_select("raw_data_key", choices=plate_keys)
            ui.update_select("normalized_data_key", choices=plate_keys)
            ui.update_select("plot_key", choices=plate_keys)

    @output
    @render.ui
    def process_status_message():
        if rv_processing():
            return ui.p(ui.strong("Processing data... Please wait.", style="color: blue;"))
        elif rv_processed_data() is not None:
            return ui.p(ui.strong("Data processed successfully.", style="color: green;"))
        else:
            return ui.p(ui.strong("Click 'Process Data' to begin.", style="color: red;"))

    @output
    @render.ui
    def download_status_message():
        if rv_downloading_results():
            return ui.p(ui.strong("Preparing download... Please wait.", style="color: blue;"))
        elif rv_processed_data() is not None:
            return ui.p(ui.strong("IC50 Results ready for download.", style="color: green;"))
        else:
            return ui.p(ui.strong("Process data to enable.", style="color: red;"))

    @output
    @render.ui
    def plot_download_status_message():
        if rv_downloading_plots():
            return ui.p(ui.strong("Preparing plots for download... Please wait.", style="color: blue;"))
        elif rv_processed_data() is not None:
            return ui.p(ui.strong("Plots are ready for download.", style="color: green;"))
        else:
            return ui.p(ui.strong("Process data to enable.", style="color: red;"))

    @output
    @render.table
    def raw_data_table():
        if rv_processed_data() is None or input.raw_data_key() == "":
            return None
        return rv_processed_data()["virus_plate_rawData"][input.raw_data_key()]

    @output
    @render.table
    def normalized_data_table():
        if rv_processed_data() is None or input.normalized_data_key() == "":
            return None
        return rv_processed_data()["normalized_data"][input.normalized_data_key()]

    @output
    @render.ui
    def ic50_plot():
        if rv_processed_data() is None or input.plot_key() == "":
            return None
        elif input.plot_type() == "Control":
            image_base64= rv_processed_data()["postv_ctrl_Figures"][input.plot_key()]
        else:
            image_base64= rv_processed_data()["sample_IC50_Figures"][input.plot_key()]
        
        return ui.HTML(f'<img src="data:image/png;base64,{image_base64}">') # Can further manipulate the 'canvas' dimensions and such by including a style input
        

    @output
    @render.data_frame
    def results_table():
        if rv_processed_data() is None:
            return None
        if input.results_type() == "Control":
            return rv_processed_data()["aggregate_control_data"].reset_index()
        else:
            return rv_processed_data()["aggregate_sample_data"].reset_index()

    @render.download(filename="aggregated_results.xlsx")
    def download_results():
        if rv_processed_data() is None:
            return None
        
        rv_downloading_results.set(True)
        
        output= io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            rv_processed_data()["aggregate_control_data"].to_excel(writer, sheet_name='Control Sera Results Data')
            rv_processed_data()["aggregate_sample_data"].to_excel(writer, sheet_name='Sample Results Data')
        
        output.seek(0)
        ui.notification_show("Aggregated results downloaded successfully!", type="message")
        rv_downloading_results.set(False)
        return output
    
    @render.download(filename="all_plots.zip")
    def download_all_plots():
        if rv_processed_data() is None:
            return None
        
        rv_downloading_plots.set(True)
        
        zip_buffer= io.BytesIO()
        
        with ZipFile(zip_buffer, 'w') as zipf:
            for plot_type in ["postv_ctrl_Figures", "sample_IC50_Figures"]:
                for key, image_base64 in rv_processed_data()[plot_type].items():
                    image_data= base64.b64decode(image_base64)
                    file_name= f"{plot_type}_{key}.png"
                    zipf.writestr(file_name, image_data)
        
        zip_buffer.seek(0)
        ui.notification_show("All plots downloaded successfully!", type="message")
        rv_downloading_plots.set(False)
        return zip_buffer

#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-

app= App(app_ui, server)

if __name__ == "__main__":
    app.run()

#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-#+-
