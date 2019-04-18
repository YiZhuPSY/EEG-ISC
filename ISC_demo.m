%% load parameters
clear;clc

allsub_No = [1:15]; 
badsub_No = [];
sub_No = setdiff(allsub_No,badsub_No); 
n = length(sub_No); 

fs = 250; 


duration = [100;81;149;129;133;67;123;57;121;133;132;170;215;170;129];
allvideo = 1:15;


%% ISC for single video 

disp('3. It is time to run ISC!')

load_dir = 'G:/Research/data/6. dataforISC_5sp_5SD_40ms/full duration';
save_dir = 'G:/Research/data/8. ResultsforISC_5sp_5SD_40ms/full duration';

% initialize following variables
Nchan = 60; % number of channels
Ncomp = 3; % number of component
Nvid = numel(allvideo); % number of videos
isc = zeros(Nvid,1);
isc_percomp = zeros(Ncomp,Nvid);
isc_persubject = zeros(Nvid,n); 
isc_persubject_percomp = zeros(Ncomp,Nvid,n);  

w = zeros(Nvid,Nchan,Ncomp); 
a = zeros(Nvid,Nchan,Ncomp);

for videoi = 1:15
     videoi
     
     cd(load_dir)
     datafile = ['v',num2str(videoi),'.mat'];
     
     [ISC,ISC_persubject,~,W,A] = runisc(datafile); %setpath runisc.m
     
     isc(videoi,1) = sum(ISC(1:3,1));
     isc_persubject(videoi,:) = sum(ISC_persubject(1:3,:));
    
    for compi = 1:3
        isc_percomp(compi,videoi) = ISC(compi,1);
        isc_persubject_percomp(compi,videoi,:) = ISC_persubject(compi,:);
        w(videoi,:,compi) = W(:,compi);
        a(videoi,:,compi) = A(:,compi);
    end
end

cd(save_dir)
save('ISC_allvideo.mat','isc','isc_persubject','isc_percomp','isc_persubject_percomp','w','a')

