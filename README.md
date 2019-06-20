# EEG-ISC
This is matlab code to compute intersubject correlation (ISC) in EEG using correlated component analysis, 
specified for analyses presented in the project "Learning desire is predicted by similar neural processing of naturalistic educational materials".

You will need EEGLAB and EEG data collected from NeuroScan to run Step1_preprocess_demo.m.

You will need following files to run Step2_ISC_demo.m:
1. runisc.m -- EEG-ISC specific code
2. topoplot.m -- stand-alone version of EEGlab's popular display function (does not require EEGlab). 
3. Neuroscan64.loc -- Neuroscan location file for topoplot
4. notBoxPlot.m -- stand-alone version of Rob Campbell's scatter plot (does not require suplementary functions)
5. v12.mat -- EEG data with 60 electrodes from 15 subjects while watching the course clip No.12. (time-points * channels * subjects)
