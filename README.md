# EEG-ISC
This is matlab code to compute intersubject correlation (ISC) in EEG using correlated component analysis, 
specified for analyses presented in the project "Learning desire is predicted by similar neural processing of naturalistic educational materials" carried out by Yi Zhu, Yafeng Pan and Yi Hu.

You will need EEGLAB and curry7 data collected from NeuroScan to run Step1_preprocess_demo.m.

You will need following files to run Step2_ISC_demo.m:
1. runisc.m -- EEG-ISC specific code
2. topoplot.m -- stand-alone version of EEGlab's popular display function (does not require EEGlab). 
3. Neuroscan64.loc -- Neuroscan location file for topoplot
4. notBoxPlot.m -- stand-alone version of Rob Campbell's scatter plot (does not require suplementary functions)
5. v1.mat -- Yi Zhu's EEG data with 60 electrodes from 15 subjects while watching the course clip No.1. (time-points * channels * subjects) (v1.mat, v2.mat, v3.mat, v4.mat, v5.mat, v6.mat, v7.mat, v8.mat, v9.mat, v10.mat, v11.mat, v12.mat, v13.mat, v14.mat, v15.mat, you could ask authors for the data files)


(c) Yi Zhu, zhuyi860574@gmail.com
