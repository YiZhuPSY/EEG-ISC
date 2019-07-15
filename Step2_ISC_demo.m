% This code could be applied to calculate inter-subject correlation (ISC) in EEG 
% using correlated component analysis, specified for analyses presented in the project 
% "Learning desire is predicted by similar neural processing of naturalistic educational materials"
 
% To run this code appropriately, you need to set following files to the path
% runisc.m -- EEG-ISC specific code
% topoplot.m -- stand-alone version of EEGlab's popular display function (does not require EEGlab).
% Neuroscan64.loc -- Neuroscan location file for topoplot
% notBoxPlot.m -- stand-alone version of Rob Campbell's scatter plot (does not require suplementary functions)

% Use sampledata -- v12.mat to have a try!!!
% v12.mat -- EEG data while watching the first ten seconds of course clip No.12.(time-points * channels * subjects)

% Due to the data size limitation for upload, full-length data is available from the authors on reasonable request.

%% load parameters
clear;clc

allsub_No = [1:15]; 
badsub_No = []; 
sub_No = setdiff(allsub_No,badsub_No); 
n = length(sub_No);

%% ISC analysis for each video
% initialize following variables
Nchan = 60; % number of channels
Ncomp = 3; % number of component
Nvid = 1; % number of videos

isc = zeros(Nvid,1);
isc_percomp = zeros(Ncomp,Nvid);
isc_persubject = zeros(Nvid,n); 
isc_persubject_percomp = zeros(Ncomp,Nvid,n);  

w = zeros(Nvid,Nchan,Ncomp); 
a = zeros(Nvid,Nchan,Ncomp);

datafile = 'v12.mat';

[ISC,ISC_persubject,W,A] = runisc(datafile); 

isc(1,1) = sum(ISC(1:3,1));
isc_persubject(1,:) = sum(ISC_persubject(1:3,:));
    
    for compi = 1:3
        isc_percomp(compi,1) = ISC(compi,1);
        isc_persubject_percomp(compi,1,:) = ISC_persubject(compi,:);
        w(1,:,compi) = W(:,compi);
        a(1,:,compi) = A(:,compi);
    end
    
save('ISC_v12.mat','isc','isc_persubject','isc_percomp','isc_persubject_percomp','w','a')

%% topoplot
Ncomp = 3;

if ~exist('topoplot') | ~exist('notBoxPlot')
    warning('Get display functions topoplot, notBoxPlot where you found this file or on the web');
else
    for i=1:Ncomp
        subplot(1,Ncomp,i);
        topoplot(squeeze(a(1,:,i)),'Neuroscan64.loc'); % video 12
        set(gca,'clim',[-0.6 0.6])
        %colorbar
    end
end