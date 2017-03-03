close all; clear all; clc;

%% ========= (2) RECEIVED SIGNAL =========

[a,Fs]=wavread('BruitBlanc'); 

%% Paramètres modifiables

M=1;                            % Number of sources
N=2;                            % Number of microphones (on doit avoir M<=N-1)
dist = .15;                     % Distance between adjacent microphones (en m)
%c = 346.287;                   % Speed of sound in air (en m/s)
c=1500;                         % Speed of sound in water
f = 7500;                      % Signal frequency (en Hz)

lfft=1024*1;                    % number of data points for FFT in a snapshot
K= floor(length(a)/lfft);       % Number of frequency snapshots (multiple of lfft)
L = K*lfft;                     % Number of data snapshots recorded by receiver
y=a(1:L,:);                     % Troncated signal (contient nbre entier de snapshots)
                     
T=1/Fs;                    
r = [(-(N-1)/2:(N-1)/2).',zeros(N,2)]; % Assume uniform linear array (on peut avoir élévation si on fait array rectangulaire avec min 3 hydrophones)


%% Forme complexe du signal reçu
df = Fs/lfft/1;                 % frequency grid size
F = 0:df:Fs/1-df;
for ih=1:N
    for iv=1:K
        pos=(iv-1)*lfft+1;
        tmp=y(pos:pos+lfft-1,ih);
        X(:,ih,iv)=fft(tmp);
    end
end
% Find the closest frequency from the discrete frequency vector resulted from FFT
[mf,mi]=min(abs(F-f)); 
f0=F(mi);

for ih=1:N
    for iv=1:K
        X0(ih,iv)=X(mi,ih,iv); % signal complexe
    end
end

%% Common parameters for all algo
Rxx = X0*X0'/L;
% Search directions
AzSearch = (0:1:180).'; % Azimuth values to search
ElSearch = zeros(size(AzSearch)); % Simple 1D example (à changer si on cherche aussi élévation)

% Corresponding points on array manifold to search
kSearch = pi*[cosd(AzSearch).*cosd(ElSearch), ...
          sind(AzSearch).*cosd(ElSearch), sind(ElSearch)].';
ASearch = exp(-1j*r*kSearch);

%% Bartlett or Conventional Beamforming

%Bartlett spectrum
for i=1:length(AzSearch)
    Z=ASearch(:,i)'*Rxx*ASearch(:,i);
    Zbf(i)=abs(Z);
end

%% Capon

% Capon spectrum
for i=1:length(AzSearch)
    Z=ASearch(:,i)'*inv(Rxx)*ASearch(:,i);
    Zc(i)=abs(1/Z);
end


%% MUSIC

% Eigendecompose
[E,D] = eig(Rxx);
[lambda,idx] = sort(diag(D)); % Vector of sorted eigenvalues
E = E(:,idx); % Sort eigenvalues accordingly
En = E(:,1:end-M); % Noise eigenvectors (ASSUMPTION: M IS KNOWN)


% MUSIC spectrum
for i=1:length(AzSearch)
    Z=(ASearch(:,i)'*En)*En'*ASearch(:,i);
    Zm(i)=abs(1/Z);
end

%% Plot

figure();
plot(AzSearch,10*log10(Zbf/max(Zbf)),AzSearch,10*log10(Zc/max(Zc)),AzSearch,10*log10(Zm/max(Zm)));
title('Comparison between 3 different aglorithms of localization');
xlabel('Azimuth (degrees)');
ylabel('Spectrum (dB)');
legend('Bartlett','Capon','MUSIC');
grid on; axis tight;

