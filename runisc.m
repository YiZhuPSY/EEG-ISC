function [ISC,ISC_persubject,W,A] = runisc(datafile)

% some ISC processing parameters
gamma = 0; % shrinkage parameter; smaller gamma for less regularization 
Nsec  = 5;  % time-window (in seconds) over which to compute time-reposeved ISC
Ncomp = 3;  % number of components to dispaly (all D are computed)
    
% If you have no data, use ours for demo. 
if nargin<1, datafile = 'v12.mat'; disp('analyzing video 12'); end 

if exist(datafile) 
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

% now start the ISC code proper

% compute cross-covariance between all subjects i and j
Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]); % check this line!!!

% compute within- and between-subject covariances
Rw =       1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

% shrinkage regularization of Rw         
Rw_reg = (1-gamma)*Rw + gamma*mean(eig(Rw))*eye(size(Rw));

% compute correlated components W using regularized Rw, sort components by ISC
[W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);

% compute forward model ("scalp projections") A
%A=Rw_reg*W*inv(W'*Rw_reg*W);
A=Rw*W*inv(W'*Rw*W);


% Compute ISC resolved by subject, see Cohen et al.
for i=1:N
    Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij(:,:,i,i)+Rij(:,:,j,j)); end; end
    Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij(:,:,i,j)+Rij(:,:,j,i)); end; end
    ISC_persubject(:,i) = diag(W'*Rb*W)./diag(W'*Rw*W);
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
% Xr = phaserandomized(X);
