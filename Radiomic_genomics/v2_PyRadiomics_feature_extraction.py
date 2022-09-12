#!/usr/bin/env python
# coding: utf-8

# # PyRadiomics feature extraction
# This code has been created to calculate the radiomic features of a set of radiomic images of different patients. 

# In[1]:


import os

from radiomics import featureextractor 
import SimpleITK as sitk

import numpy as np
import nibabel as nib

import csv
import pandas as pd


# ***
# ## General instructions
# The main function to implement to calculate the radiomic features from the images of a set of patients is **Radiomic_features**.
# 
# **Radiomic_features()** takes 3 arguments: 1) path of the input folder, 2) path of the output folder 3) path of the parameters file.
# 
# ### In-depth analysis of the required arguments
# ***1. Path to the input folder*** The function requires two images in NIFTI format for each patient: the original and its mask. The name of the original image can be any, but the name of the mask file must contain the word "mask" in order to differentiate it. The pair of images (original and mask) for each patient must be stored in a unique folder named with the patient ID. All patient folders will be stored at the same time inside the input folder. The function needs as first argument the path of this input folder containing all patient folders.
# 
# ***2. Path of the output folder*** The function needs as second argument the path of the folder where the files with the results are to be stored. The Radiomic_features function returns one *id_patient.tsv* file for each patient with its radiomic features and another *all_patient_features.tsv* with the radiomic features of all patients joined in the same file. The first column of *all_patient_features.tsv* corresponds to the name of the radiomic features, the other columns are the values of the radiomic features for each patient, the first row being the patient ID.
# 
# ***3. Path of the parameters file*** All 4 categories of customization can be provided in a single yaml or JSON structured text file. Find attached an example. The path of the *params* files must be the last argument to the function.
# 
# ***

# In[2]:


# Create function to instantiate teh extractor

def inizialize_extractor(path_to_parameters):
    return featureextractor.RadiomicsFeatureExtractor(path_to_parameters)


# In[3]:


# Function to print the features that will be calculated
def print_enabled_features(extractor):
    radiomic_features = list(extractor.enabledFeatures.keys())
    d = {'Enabled features': radiomic_features}
    print(pd.DataFrame(d))


# In[4]:


# Fucntion to print the settings defined to the extractor
def print_extractor_settings(extractor):
    settings, values = list(extractor.settings.keys()), list(extractor.settings.values())
    d = {'Settings': settings, 'Values': values}
    print(pd.DataFrame(d))


# In[5]:


# Function to create a list of patients ids
def list_patients(path_to_patients):
    return os.listdir(path_to_patients)


# In[6]:


# Function to get the original image and the mask
def get_image_mask(patient, path_to_patients):
    path_to_patient_files = path_to_patients + '/' + patient
    files_for_patient = os.listdir(path_to_patient_files)
    if len(files_for_patient) != 2:
        print('Patient ',patient,'missing image or mask')
        return False
    else:
        for file in files_for_patient:
            if 'mask' in file:
                path_to_mask = str(path_to_patient_files + '/' + file)
                path_to_mask = path_to_mask
            else:
                path_to_original = path_to_patient_files + '/' + file
        return path_to_original, path_to_mask


# In[7]:


# Function to calculate the radiomic features of a patient
def calculate_radiomic_features(extractor, path_to_original_image, path_to_mask):
    return extractor.execute(path_to_original_image, path_to_mask)


# In[8]:


# Function to save radiomic extraction results to a .tsv
def save_results(result, path_to_output_dir, patient):
    path_to_output_file = path_to_output_dir + '/' + patient + '.tsv'
    
    # Storing the radiomic features in the .tsv file
    with open(path_to_output_file, "w") as outfile:
        # Create the .tsv file
        csvwriter = csv.writer(outfile, delimiter='\t')
        
        # Iteration over all the values and keys of the extraction results
        for key, value in result.items():
            obj_to_write = []
            obj_to_write.append(str(key))
            obj_to_write.append(str(value))
            csvwriter.writerow(obj_to_write)
    print('Successful radiomic extraction for patient:\t', patient)


# In[9]:


# Function to merge all .tsv into a unique .tsv file for all patients
def merge_tsv_files(general_file, input_file, patient_id):
    df_general = pd.read_csv(general_file,header=None, sep='\t')
    df_input = pd.read_csv(input_file, header=None, sep='\t')

    df_general[patient_id] = df_input[1]
    
    df_general.to_csv(general_file, header=None, index=False, sep='\t')


# In[10]:


# Function to generate a unique file with the radiomic features of all the patients
def generate_final_file(patients_list, path_to_output_dir):
    path_to_final_file = path_to_output_dir + '/all_patient_features.tsv'
    # Create a .tsv where to merge the values for all patients
    base_file = path_to_output_dir + '/' + patients_list[0] + '.tsv'
    df = pd.read_csv(base_file, header=None, sep='\t')
    df = df[[0]]
    df.to_csv(path_to_final_file , header=None, index=False, sep='\t')
    
    # Iteration over all the patient .tsv files to generate the global one
    for patient in patients_list:
        path_to_patient_file = path_to_output_dir + '/' + patient + '.tsv'
        merge_tsv_files(path_to_final_file, path_to_patient_file, patient)
    
    patients_list.insert(0, '')
    # Apending the header to the final .tsv file which will be the input for RadAR
    df_final = pd.read_csv(path_to_final_file, header=None, sep='\t')
    df_final.to_csv(path_to_final_file, header=patients_list, index=False, sep='\t')


# In[11]:


# Function to print the summary of the process
def print_summary(id_list):
    print('Radiomic features extraction have been finished')
    
    if len(id_list) == 0:
        print('All patients could be processed')
        
    else:
        print('All patients but ', len(id_list), ' were processed. The following ids correspond to the patients that may have some issues in the NIFTI input files', id_list)


# In[12]:


def Radiomic_features(path_to_patients, path_to_output, 
                                path_to_parameters):
    # Instatntiate the extractor
    extractor = inizialize_extractor(path_to_parameters)
    print('Instatntiate the extractor: Done\n')
    
    # Print enabled parameters
    print_enabled_features(extractor)
    print('Print enabled parameters: Done\n')
    
    # Print settings and its corresponding values
    print_extractor_settings(extractor)
    print('Print settings and its corresponding values: Done\n')
    
    # Creation of a patients id list
    patients_id_list = list_patients(path_to_patients)
    print('Creation of a patients id list: Done\n')
    
    # List of patients that could not be processed
    patients_with_error = []
    print('list for errors: Created\n')
    
    # List of patients successfuly processed
    patients_processed = []
    
    # Iteration over all patients to calculate their radiomic features
    for patient in patients_id_list:
        try:
            # Get path to image and mask file
            image, mask = get_image_mask(patient, path_to_patients)
            print('image:', image, '\nmask:', mask)
        
            # Calculate radiomic features
            result = calculate_radiomic_features(extractor, image, mask)
            print('\nRadiomic features claculated for patient:',patient)
            # Save radiomic features of the patient to a .tsv file
            save_results(result, path_to_output, patient)
            
            # Append patient id to the processed patients list
            patients_processed.append(patient)
            
        except:
            patients_with_error.append(patient)
            print('\nError occurred for patient:', patient)
    
    # Create the final .tsv with the radiomic features of all patients
    if len(patients_id_list) != len(patients_with_error):
        generate_final_file(patients_processed, path_to_output)
    
    print_summary(patients_with_error)


# ## Example execution

# In[13]:


# input arguments
input_files = '/home/carlos/ISGlobal/im_mk_same_geometry/patients2'
output_files = '/home/carlos/ISGlobal/im_mk_same_geometry/results4'
parameters = '/home/carlos/ISGlobal/im_mk_same_geometry/Params.yaml'


# In[14]:


Radiomic_features(input_files, output_files, parameters)

