close all; clear all; clc;

%% ======= (1) TRANSMITTED SIGNALS ======= %

% Signal source directions
az = [35;39;127]; % Azimuths
el = zeros(size(az)); % Simple example: assume elevations zero
M = length(az); % Number of sources
L = 200; % Number of data snapshots recorded by receiver
c = 346.287;                    % Speed of sound in air
K = 100;                        % Number of frequency snapshots
Fs  = 44100;                    % Sampling frequency
T=1/Fs;                    
dist = .03;                     % Distance between adjacent microphones

m = randn(M,L); % Example: normally distributed random signals

%% ========= (2) RECEIVED SIGNAL =========

% Wavenumber vectors (in units of wavelength/2)
k = pi*[cosd(az).*cosd(el), sind(az).*cosd(el), sind(el)].';

% Array geometry [rx,ry,rz]
N = 10; % Number of antennas
r = [(-(N-1)/2:(N-1)/2).',zeros(N,2)]; % Assume uniform linear array

% Matrix of array response vectors
A = exp(-1j*r*k);

% Additive noise
sigma2 = 0.01; % Noise variance
n = sqrt(sigma2)*(randn(N,L) + 1j*randn(N,L))/sqrt(2);

% Received signal
x = A*m + n;
%%
plot(abs(x(1,:)))

%% Common parameters for all algo

% Sample covariance matrix
Rxx = x*x'/L;
% Search directions
AzSearch = (0:1:180).'; % Azimuth values to search
ElSearch = zeros(size(AzSearch)); % Simple 1D example

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
%Z = sum(abs(ASearch'*En).^2,2);
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

