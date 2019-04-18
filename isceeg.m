function [ISC,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile)
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
% synchronize supramodal neural responses, eNeuro, 3(6), November 2016.
%
% Agustin Petroni, Samantha Cohen, Nicolas Langer, Simon Henin, Tamara
% Vanderwal, Michael P. Milham, Lucas C. Parra, Children, particularly
% males, exhibit higher intersubject correlation of neural responses to
% natural stimuli than adults", in review
%
% Jason Ki, Simon Kelly, Lucas C. Parra, "Attention strongly modulates
% reliability of neural responses to naturalistic narrative stimuli.鈥�
% Journal of Neuroscience, 36 (10), 3092-3101.
%
% Jacek P. Dmochowski, Matthew A. Bezdek, Brian P. Abelson, John S.
% Johnson, Eric H. Schumacher, Lucas C. Parra, 鈥淎udience preferences are
% predicted by temporal reliability of neural processing鈥�, Nature
% Communication, 5567, July 2014.
%
% Jacek P. Dmochowski, Paul Sajda, Joao Dias, Lucas C. Parra, 鈥淐omponents
% of ongoing EEG with high correlation point to emotionally-laden attention
% -- a possible marker of engagement?鈥�, Frontiers in Human Neuroscience,
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
% 10/27/17. version 0.06: outlier rejection now based on interquartile difference
% 12/06/17. version 0.07: Jens Madsen: added automatic bad-channel rejection based on power

% some ISC processing parameters
gamma = 0.1; % shrinkage parameter; smaller gamma for less regularization
Nsec  = 5;  % time-window (in seconds) over which to compute time-reposeved ISC
Ncomp = 3;  % number of components to dispaly (all D are computed)
    
% If you have no data, use ours for demo. 
if nargin<1, datafile = 'EEGVolume.mat'; disp('Using demo data from Cohen et al. 2016'); end 

if exist(datafile) 
    load(datafile,'X','fs','eogchannels','badchannels');
    for i=1:size(X,3), badchannels{i} = -1; end
else
    warning('Can not find data file. Using random data.') 
    fs=256;                 % *** Samling rate needed for filtering
    X=randn(30*fs,64,15);   % *** EEG volume, recommend T>100000, D>16, N>10 
    % *** List of bad channels (one vector of indice per subject). 
    % Here set to each to -1 for automatic bad-channel detection
    for i=1:size(X,3), badchannels{i} = -1; end
    eogchannels=[];         % *** Index of EOG channels, if any; will regress them out and exclude
end; 

% standard eeg preprocessing (see function below). Will discard EOG channels 
X = preprocess(X,eogchannels,badchannels,fs);

% T samples, D channels, N subjects
[T,D,N] = size(X);  

% now start the ISC code proper

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); 

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
A=Rw*W/(W'*Rw*W);

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
if ~exist('topoplot') | ~exist('notBoxPlot')
    warning('Get display functions topoplot, notBoxPlot where you found this file or on the web');
else
    for i=1:Ncomp
        subplot(2,Ncomp,i);
        topoplot(A(:,i),'BioSemi64.loc','electrodes','off'); title(['a_' num2str(i)])
    end
    subplot(2,2,3); notBoxPlot(ISC_persubject(1:Ncomp,:)'); xlabel('Component'); ylabel('ISC'); title('Per subjects');
    subplot(2,2,4); plot(ISC_persecond(1:Ncomp,:)'); xlabel('Time (s)'); ylabel('ISC'); title('Per second');
end

% ### Run all the code above with phase-randomized Xr to get chance values
% of ISC measures under the null hypothesis of no ISC.
Xr = phaserandomized(X);




% -------------------------------------------------------------------------
function X = preprocess(X,eogchannels,badchannels,fs)
% All the usual EEG preprocessing, except epoching and epoch rejection as
% there are not events to epoch for natural stimuli. duh! Instead, bad data
% is set = 0 in the continuous stream, which makes sense when computing
% covariance matrices but maybe not for other purposes. Bad channels are
% removed for all the indice given in badchannels cell array (each subject
% has its own vector of indice). None are removed if this is set to []. If
% it is set to -1, channels are removed based on outlies in power.

debug = 0;     % turn this on to show data before/after preprocessing. 
kIQD=4;        % multiple of interquartile differences to mark as outliers samples
kIQDp=3;       % multiple of interquartile differences to mark as outliers channels
HPcutoff =0.5; % HP filter cut-off frequequency in Hz

% pick your preferred high-pass filter
[z,p,k]=butter(5,HPcutoff/fs*2,'high'); sos = zp2sos(z,p,k);

[T,D,N]=size(X); 

% if it is not EOG, then it must be EEG channel
eegchannels = setdiff(1:D,eogchannels);

% Preprocess data for all N subjects
for i=1:N
    
    data = X(:,:,i);

    % remove starting offset to avoid filter trancient
    data = data-repmat(data(1,:),T,1);
    
    % show the original data
    if debug, subplot(2,1,1); imagesc((1:T)/fs,1:D,data'); title(['Subject ' num2str(i)]); end

    % high-pass filter
    data = sosfilt(sos,data);          
    
    % regress out eye-movements;
    data = data - data(:,eogchannels) * (data(:,eogchannels)\data);     

    % detect outliers above stdThresh per channel; 
    data(abs(data)>kIQD*repmat(diff(prctile(data,[25 75])),[T 1])) = NaN;
    
    % remove 40ms before and after;
    h=[1; zeros(round(0.04*fs)-1,1)];    
    data = filter(h,1,flipud(filter(h,1,flipud(data))));
    
    % Mark outliers as 0, to avoid NaN coding and to discount noisy channels
    data(isnan(data))=0;

    % Find bad channels based on power ourliers, if not specified "by hand"
    if badchannels{i} == -1, 
        logpower = log(std(data)); Q=prctile(log(std(data(:,eegchannels))),[25 50 75]);
        badchannels{i} = find(logpower-Q(2)>kIQDp*(Q(3)-Q(1)));  
    end
    
    % zero out bad channels
    data(:,badchannels{i})=0; 
    
    % show the result of all this
    if debug, subplot(2,1,2); imagesc((1:T)/fs,1:D,data'); caxis([-100 100]); xlabel('Time (s)'); drawnow; end

    X(:,:,i) = data;
    
end

% remove the eog channels as we have used them already
X = X(:,eegchannels,:);


% -------------------------------------------------------------------------
function Xr = phaserandomized(X);
% Generate phase randomized surrogate data Xr that preserves spatial and
% temporal correlation in X, following Prichard D, Theiler J. Generating 
% surrogate data for time series with several simultaneously measured 
% variables. Physical review letters. 1994 Aug 15;73(7):951.

[T,D,N] = size(X);

Tr = round(T/2)*2; % this code only works if T is even; make it so
for i = 1:N
    Xfft = fft(X(:,:,i),Tr); % will add a zero at the end if uneven length
    Amp = abs  (Xfft(1:Tr/2+1,:)); % original amplitude
    Phi = angle(Xfft(1:Tr/2+1,:)); % orignal phase
    Phir = 4*acos(0)*rand(Tr/2-1,1)-2*acos(0); % random phase to add
    tmp(2:Tr/2,:) = Amp(2:Tr/2,:).*exp(sqrt(-1)*(Phi(2:Tr/2,:)+repmat(Phir,1,D))); % Theiler's magic
    tmp = ifft([Xfft(1,:); tmp(2:Tr/2,:); Xfft(Tr/2+1,:); conj(tmp(Tr/2:-1:2,:))]); % resynthsized keeping it real
    Xr(:,:,i) = tmp(1:T,:,:); % grab only the original length
end

