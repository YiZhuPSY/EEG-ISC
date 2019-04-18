# EEG-ISC
This is matlab code to compute intersubject correlation (ISC) in EEG using correlated component analysis, 
specified for analyses presented in Zhu, Y., Pan, Y., & Hu, Y. (2019). Learning desire is predicted by similar neural processing of naturalistic educational materials. eNeuro.

You will need EEGLAB to run Step1_preprocess_demo;

You will need these files to run Step2_ISC_demo:
1. runisc.m -- EEG specific code
2. topoplot.m -- stand-alone version of EEGlab's popular display function (does not require EEGlab). 
3. Neuroscan64.loc -- Neuroscan location file for topoplot
4. notBoxPlot.m -- stand-alone version of Rob Campbell's scatter plot (does not require suplementary functions)
5. v12.mat -- Yi Zhu's EEG data with 64 electrodes from 15 subjects while watching the course clip No.12. 
(v1.mat, v2.mat, v3.mat, v4.mat, v5.mat, v6.mat, v7.mat, v8.mat, v9.mat, v10.mat, v11.mat, v12.mat, v13.mat, v14.mat, v15.mat).
