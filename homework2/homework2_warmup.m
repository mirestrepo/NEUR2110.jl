cd('C:\Dropbox\workDocs\Teaching\NEUR2110_Statistical_Neuroscience_Fall_2017\Homework_2_Lect_3&4\')
% White noise: power spectrum (periodogram and multitaper; integrated power and variance )==========================
clear all, close all
Fs=2^10;%sampling rate
dt = 1/Fs;
T=2;%Time window
df = 1/T; %frequency resolution
fNyquist=Fs/2;
N=T*Fs;
t=[0:N-1]'*dt;

F = [0:df:fNyquist];%frequency grid

% Generate a realization from GWN (zero-mean, variance=sigma^2)
sigma2=10;
X = sqrt(sigma2)*randn(N,1); 

% Periodogram (rectangular window): Discrete Fourier transform (via Fast Fourier Transform, FFT)
X = X - mean(X); %subtract mean: zero DC shift
Xf = fft(X);
% One-sided power spectrum (periodogram)
Sp = dt^2 * 1/T * abs(Xf).^2; % dt^2 * 1/T * Xf.*conj(Xf);
Sp = 2 * Sp(1:N/2+1); %One-sided spectrum (i.e. only positive frequencies) for even N

% Power spectrum using the Hann taper
H=hanning(N);
Xf=fft(H.*X);    
Sh = dt^2 * 1/T * abs(Xf).^2; % dt^2 * 1/T * Xf.*conj(Xf);
Sh = 2 * Sh(1:N/2+1); %One-sided spectrum (i.e. only positive frequencies) for even N

% multitaper
R = 10;%bandwidth in Hz
NFFT=N;
nTapers=[];
removeTemporalMean=true;
removeEnsembleMean=true;
% [Pxx, Pyy, Pxy, XYphi, Cxy, F, nTapers]= multitaperSpectrum(X,Y,Fs,bandWidth,NFFT,removeTemporalMean,RemoveEnsembleMean,nTapers)
[Smt, ~, ~, ~, ~, ~, nTapers]= multitaperSpectrum(X,X,Fs,R,NFFT,removeTemporalMean,removeEnsembleMean,nTapers);

figure(1),clf,set(gcf,'color',[1 1 1])
subplot(311)
plot(t,X,'k');
xlabel('Time (s)')
ylabel('X_t')
title(['Gaussian white noise; zero mean; variance = ',num2str(sigma2)])
box off
subplot(312)
plot(F,Sp,'k',F,Sh,'g',F,Smt,'r.-')
xlim([F(1) F(end)])
xlabel('Frequency (Hz)')
ylabel('Power (^2/Hz)')
title(['Integrated power (variance): ',num2str(sum(Sp)*df)])
box off
subplot(313)
% X is zero mean, thus DC component for the periodogram is zero. For dB we do not plot it
plot(F(2:end),10*log10(Sp(2:end)/max(Sp)),'k',F(2:end),10*log10(Sh(2:end)/max(Sh)),'g',F,10*log10(Smt/max(Smt)),'r.-')
xlabel('Frequency (Hz)')
ylabel('dB')
legend('Periodogram','Hann taper','Multitaper')
xlim([F(1) F(end)])
box off


% Confidence intervals based on the chi2 CDF
alpha=0.05;
nTrials=1;
nTapers=5;
S = ...% power spectrum (a vector)
dof=2*nTapers*nTrials;
q1=chi2inv(alpha/2,dof);
q2=chi2inv(1-alpha/2,dof);
CI(1,:)=dof*S/q2;
CI(2,:)=dof*S/q1;


% Spectrogram: skeleton code ================================================
N = ... % total number of samples
n = ... % number of samples in each moving time window
% t is time for each of the N samples in the LFP variable
j=0;
for k=n+1:round(n/2):N-n
    j=j+1;
    tt(j)=t(k);
    x=LFP(k-n+1:k);x=x(:);% multitaperSpectrum.m expects column vectors ...    
    % periodogram (rectangular taper)
    % ... compute periodogram power spectrum and store    
    Sp(:,j) = ... %  store peridogram for the j window; One-sided spectrum (i.e. only positive frequencies)            
    
    % Hanning taper
    % ... compute Hann taper power spectrum and store
    Sh(:,j) = ... store Hann taper power spectrum for the j window; One-sided spectrum (i.e. only positive frequencies)            
        
    % Multitaper    
    % ... compute multitaper power and store
    
    Smt(:,j) = ... % store multitaper power spectrum for the j window; One-sided spectrum (i.e. only positive frequencies)            
end


% Theoretical spectrum for an univariate AR(p) =============================
% Example: AR(6)
clear all
order=6;
A=[3.9515 -7.8885 9.7340 -7.7435 3.8078 -0.9472];
sigma2=1;
burnin = 1000;
Fs = 2^10; dt = 1/Fs;
T=1; N=T*Fs;
X = zeros(N+burnin,1);
for t=order+1:N+burnin
    X(t) = A * X(t-1:-1:t-order) + sqrt(sigma2) * randn;
end

% X=filter(1,[1 -A],sigma* randn(N+burnin,1)); % Alternative implementation using the function filter.m
X=X(burnin+1:end);
t=[0:N-1]*dt;

powerspect =@(f,sigma2,a) sigma2 ./ abs(1 + a * exp(-1i*pi*f*[2:2:length(a)*2]'))^2;
df=0.001;
f=-0.5:df:.5;
delta=1; %delta in seconds; To get the units right (1/Hz)...
for k =1:length(f)        
    Sf(k) = delta * 1/N * powerspect(f(k),sigma2,-A); % without multiplication by 2 to match var(X)
    %Sf(k) = 2* delta * 1/N * powerspect(f(k),sigma2,-A); % 2 * delta^2 * 1/T = 2 * delta^2 * 1/(N * delta) =  2 * delta * 1/N
    %Sf(k) = delta * powerspect(f(k),sigma2,-A); %This is what Babadi/Brown plot in Figure 1c, no normalization by N
end
figure(1),clf,set(gcf,'color',[1 1 1])
% f * Fs = [0 ... Nyquist freq]
plot(f*Fs,10*log10(Sf),'k','linewidth',2)
axis([0 Fs/2 -80 30])
grid on
xlabel('Frequency [Hz]')
ylabel('dB')

[var(X) xcov(X,0,'unbiased') sum(Sf)]% The sample variance varies around the theoretical value




