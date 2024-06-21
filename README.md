# Mouth-Based IMU Feature Engineering

This repository contains scripts (MATLAB and python variants) for processing and extracting features from mouth-based IMU sensors, as inspired by Gabler 2020. The features extracted are intended to be used in future classification models, with the specific classification depending on the project's aims.

## Features Extracted

The script extracts and summarizes the following types of features:

1. **Statistical Summaries**
   - Mean, min, max, standard deviation (SD), variance, and the ratio of linear/angular resultant.
2. **Pulse Size Summaries**
   - Width and area under the curve (AUC) of resultant traces.
3. **Power Spectral Density (PSD) Summaries**
   - PSD features including values at which 95% power within the signal is achieved and frequency content in 20Hz bins up to 200Hz+.
4. **Kinematic-Based Summaries**

## Data Assumptions

- The data is formatted in individual `.csv` files and saved within the file directory.
- Each `.csv` file should have the following column structure:
  1. Time
  2. Rotational acceleration (x, y, z, resultant)
  3. Rotational velocity (x, y, z, resultant)
  4. Linear acceleration (x, y, z, resultant)

## Script Overview

1. **Data Import**
   - Imports data from `.csv` files and extracts specific axes.
   - Assumes data is raw.

2. **Statistical Summaries**
   - Calculates trial mean, min, max, SD, variance, and range for all variables.

3. **Pulse Size Summaries**
   - Defines pulse duration using a threshold of trial mean as the start and end points.
   - Calculates AUC for linear and rotational acceleration.

4. **PSD Summaries**
   - Calculates power spectral density summaries and features.

5. **Output Data Table**
   - Constructs output data table with all extracted features.
   - Normalizes the values for additional features.
   - Exports the data to a `.csv` file.
