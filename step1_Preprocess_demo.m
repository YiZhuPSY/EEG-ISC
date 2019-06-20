% This code was based on EEGLAB funtions, which could be applied to preprocess 
% the 64-channel EEG data recorded from NeuroScan system via curry7 software,
% including high-pass filter(>1hz), notch-filter(50hz), downsample(to 250hz),
% epoch extraction, eye-movement artifact rejection,badchannel rejection,
% and outlier rejection.

% This code might not support older versions of EGGLAB (version 14.1.0) and MATLAB (R2016b). 

% To run this code appropriately, you need to set EEGLAB to the path, which can be
% downloaded from https://sccn.ucsd.edu/eeglab/download.php.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNZIP THE DOWNLOADED FILE "EEG-ISC-master.zip" INTO
% (e.g. C:\Users\hp\Desktop\EEG-ISC-master)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use SampleData.set&SampleData.fdt to have a try!!! 

clear;clc
addpath ('C:\Users\hp\Desktop\EEG-ISC-master')
%% loading parameters
%allsub_No = [1:15]; 
%badsub_No = []; 
%sub_No = setdiff(allsub_No,badsub_No); 
%n = length(sub_No); 

%video =1:15; 
%duration = [100;81;149;129;133;67;123;57;121;133;132;170;215;170;129]; 
%mark = [3:17]';  % mark=video+2

eogchannels = [63,64]; 
m_eogchannels = [33,43,63,64]; 
fs = 250; % sampling rate

%% load curry7 data recorded from NeuroScan 

% EEG = loadcurry('Acquisition 01.dap', 'CurryLocations', 'False'); 
% The size of the raw data (i.e. "Acquisition 01.dap") was too large to upload, 
% which is available from the authors on reasonable request

%% filtering (high pass 1hz)¡¢notching (50hz)¡¢downsampling (250hz) 
% filtering
% EEG = pop_eegfiltnew(EEG, [],1,3300,1,[],0);

% notching
% EEG = pop_eegfiltnew(EEG, 49,51,3300,1,[],0);

% downsampling to 250hz
% EEG = pop_resample( EEG, 250); 
% EEG = pop_saveset( EEG, 'filename','s1.set');

%% Epoch for ISC
%for videoi = 1:15 
%    EEG = pop_loadset('filename','s1.set'); 
%%% epoching
%    EEG = pop_epoch( EEG, { num2str(mark(videoi,1)) }, [0 duration(videoi,1)], 'newname',...
%        ['s1v',num2str(videoi),'.set'], 'epochinfo', 'yes');
%    EEG = pop_saveset( EEG, 'filename',['s1v',num2str(videoi),'.set'],...
%        'filepath',Epoch4ISC_dir);     
%    clear EEG
%end

%% Regressing out eye-movenments
EEG = pop_loadset('filename','SampleData.set');

data = double(EEG.data)';    
data = data - data(:, eogchannels) * (data(:, eogchannels) \ data); 
EEG.data = single(data)';
EEG = pop_saveset( EEG, 'filename','regressEOG.set'); 

        
%% Clearing badchannels 
% EEGLAB -- Automatic Channel Rejection -- 'spec'£¬5
[~,indelec] = pop_rejchan(EEG, 'elec',[1:32 34:42 44:62] ,'threshold',5,'norm','on','measure','spec');
 
badchannels = indelec';  
misplaced_channels = badchannels;
index1 = find(misplaced_channels > 32 & misplaced_channels < 42); 
index2 = find(misplaced_channels > 41);
misplaced_channels(index1) = misplaced_channels(index1)+1;
misplaced_channels(index2) = misplaced_channels(index2)+2;
badchannels = misplaced_channels;
        
EEG.data(badchannels,:) = 0;
EEG = pop_saveset( EEG, 'filename','rejbadcha.set');

%% Rejecting outliers (over 4 std)

[D,T] = size(EEG.data); %D:channels£»T£ºtimepoint
data = double(permute(EEG.data,[2,1])); 
        
stdThresh = 4;
% data(abs(data)>stdThresh*repmat(std(data),[T 1])) = NaN;  
data(abs(data-repmat(mean(data,1),[T 1]))>stdThresh*repmat(std(data),[T 1])) = NaN; 

% remove 40ms before and after
% h = [1; zeros(round(0.04*fs)-1,1)];
% remove 50ms before and after
h = [1; zeros(round(0.05*fs)-1,1)]; 
data = filter(h,1,flipud(filter(h,1,flipud(data))));  

data(isnan(data))=0;  % Mark outliers as 0, to avoid NaN 

EEG.data = single(permute(data,[2,1]));
EEG = pop_saveset( EEG, 'filename','rejoutliers');     