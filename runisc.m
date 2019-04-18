function [ISC,ISC_persubject,ISC_persecond,W,A] = runisc(datafile)
%
% This is matlab code to compute intersubject correlation (ISC) in EEG
% using correlated component analysis, as described in the references below
% (description closest to the present code is in Cohen et al, 2016).
%
% You will need these files to run this code:
%    topoplot.m    -- EEGlab's popular display function. 
%    BioSemi64.loc -- BioSemi location file for topoplot
%    notBoxPlot.m  -- Rob Campbell's scatter plot
%    EEGVolume.mat -- Samantha Cohen's EEG data with 64 electrodes from 18
%                     subjects while watching this video:
%                     https://www.youtube.com/watch?v=DEo-rw2TeAU.
%
% A stand-alone version of these files can be found on the same location as
% this file. Once you have all this the code can be run as:
%
% [ISC,ISC_persubject,ISC_persecond,W,A] = isceeg('EEGVolume.mat');
%
% To understand the format of the data file take a look at lines maked with
% ***. You will have to edit these lines if your data is in a different
% format.
%
% The input to this function is EEG data (volumes of samples by channels by
% subjects). The output are component projection vectors W, scalp
% projections A, and ISC_persubject and ISC_persecond, i.e. the
% intersubject correlation computed resolved for each subject, or resolved
% in time. Components are sorted by how much correlation between subjects
% they capture. We typically only look at the first few, and often just take
% the sum of the first 3 ISC. This is not done in this code.
%
% Parameters of the algorithm are currently set to reasonable defaults.
%
% Curently this is coded assuming a single continuous dynamic stimulus (the
% example data was recored for a video of 197 seconds duration). In
% practice I strongly recommend running this by combining data from
% multiple stimuli to get more robust estimates of the required within- and
% between-subjects covariances. Comments on how to change the code to do
% this are marked with +++
%
% Code for preprocessing EEG and generating surrogate data may be useful in
% other applications and has been separated as sub-functions. 
%
% Surrogate data can be used to test for statistical significance of ISC
% values. However, the code is not currently used for this purpose and is
% only provided as example to generate a single surrogate. See ###
%
% References:
%
% Samantha Cohen, Lucas C Parra, Memorable audiovisual narratives
% synchronize supramodal neural responses, in review at eNeuro
%
% Agustin Petroni, Samantha Cohen, Nicolas Langer, Simon Henin, Tamara
% Vanderwal, Michael P. Milham, Lucas C. Parra, Children, particularly
% males, exhibit higher intersubject correlation of neural responses to
% natural stimuli than adults", in review
%
% Jason Ki, Simon Kelly, Lucas C. Parra, "Attention strongly modulates
% reliability of neural responses to naturalistic narrative stimuli.鈥?% Journal of Neuroscience, 36 (10), 3092-3101.
%
% Jacek P. Dmochowski, Matthew A. Bezdek, Brian P. Abelson, John S.
% Johnson, Eric H. Schumacher, Lucas C. Parra, 鈥淎udience preferences are
% predicted by temporal reliability of neural processing鈥? Nature
% Communication, 5567, July 2014.
%
% Jacek P. Dmochowski, Paul Sajda, Joao Dias, Lucas C. Parra, 鈥淐omponents
% of ongoing EEG with high correlation point to emotionally-laden attention
% -- a possible marker of engagement?鈥? Frontiers in Human Neuroscience,
% 6:112, April 2012.

% (c) Lucas C Parra, parra@ccny.cuny.edu
%
% Revision history (in case you have an earlier version):
% 07/15/16. version 0.00: incomplete; got to start somewhere.
% 07/19/16. version 0.01: complete; some bugs in formulas removed.
% 07/21/16. version 0.02: cleaned up HP filter
% 08/02/16. version 0.03: removed filter bug in ourlier detection.
%                         added sample code to generate surrogate data
% 08/03/16. version 0.04: added sample data, organized into function calls
% 08/05/16. version 0.05: fixed bug in pooling over subjects, added some displays

% some ISC processing parameters
gamma = 0; % shrinkage parameter; smaller gamma for less regularization  %%%% 可更改
Nsec  = 5;  % time-window (in seconds) over which to compute time-reposeved ISC
Ncomp = 3;  % number of components to dispaly (all D are computed)
    
% If you have no data, use ours for demo. 
% if nargin<1, datafile = 'EEGVolume.mat'; disp('Using demo data from Cohen et al.'); end 
if nargin<1, datafile = 'v4.mat'; disp('analyzing video 4'); end 

if exist(datafile) 
    %load(datafile,'X','fs','eogchannels','badchannels');  
    %load(datafile,'X','fs','eogchannels','badchannels','m_eogchannels'); 
    load(datafile,'X','fs');
else
    warning('Can not find data file. Using random data.') 
    fs=256;                 % *** Samling rate needed for filtering
    X=randn(30*fs,64,15);   % *** EEG volume, recommend T>100000, D>16, N>10 
    badchannels=cell(15,1); % *** List of bad channels (one vector of indice per subject).
    eogchannels=[];         % *** Index of EOG channels, if any; will regress them out and exclude
end; 

% T samples, D channels, N subjects
[T,D,N] = size(X);  

% standard eeg preprocessing (see function below) 
%X = preprocess(X,eogchannels,badchannels,fs);

% discard the eog channels; we already used them in preprocess
% X = X(:,setdiff(1:D,eogchannels),:); D = size(X,2);
% X = X(:,setdiff(1:D,m_eogchannels),:); D = size(X,2);

% now start the ISC code proper

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); % check this line!!!

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw         
Rw_reg = (1-gamma)*Rw + gamma*mean(eig(Rw))*eye(size(Rw));

% +++ If multiple stimuli are available, then Rw and Rb should be averaged over
% stimuli here prior to computing W and A +++

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);

% compute forward model ("scalp projections") A
%A=Rw_reg*W*inv(W'*Rw_reg*W);
A=Rw*W*inv(W'*Rw*W);

% +++ If multiple stimuli are available, then Rij as computed for each stimulus
% should be used in the following to compute ISC_persubject, and
% ISC_persecond +++

% Compute ISC resolved by subject, see Cohen et al.
for i=1:N
    Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij(:,:,i,i)+Rij(:,:,j,j)); end; end
    Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij(:,:,i,j)+Rij(:,:,j,i)); end; end
    ISC_persubject(:,i) = diag(W'*Rb*W)./diag(W'*Rw*W);
end

% Compute ISC resolved in time
for t = 1:floor((T-Nsec*fs)/fs)
    Xt = X((1:Nsec*fs)+(t-1)*fs,:,:);
    Rij = permute(reshape(cov(Xt(:,:)),[D N  D N]),[1 3 2 4]);
    Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
    Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects
    ISC_persecond(:,t) = diag(W'*Rb*W)./diag(W'*Rw*W);
end

% show some results
%if ~exist('topoplot') | ~exist('notBoxPlot')
%    warning('Get display functions topoplot, notBoxPlot where you found this file or on the web');
%else
%    for i=1:Ncomp
%       subplot(Ncomp,1,i);
%        topoplot(A(:,i),'Neuroscan64.loc','electrodes','on'); title(['a_' num2str(i)])
    
%    end
    %colorbar
%    subplot(2,2,3); notBoxPlot(ISC_persubject(1:Ncomp,:)'); xlabel('Component'); ylabel('ISC'); title('Per subjects');
%   subplot(2,2,4); plot(ISC_persecond(1:Ncomp,:)'); xlabel('Time (s)'); ylabel('ISC'); title('Per second');
end

% ### Run all the code above with phase-randomized Xr to get chance values
% of ISC measures under the null hypothesis of no ISC.
%Xr = phaserandomized(X);
%%%%零假设的ISC值计算







