% % Random number generation, mean and variance ===============
rng(0);
N=100;
% Uniform [0,1]
X=rand(N,1);
N=100;
sigma2=2;
% Gaussian rv (zero-mean, variance sigma2
X=sqrt(sigma2)*randn(N,1);
mu=mean(X);
s2= var(X);

% Gaussian random vectors, covariance and correlation matrices ======
N=10000;
mu=zeros(3,1);
C = [4 0 1; 0 5 0; 1 0 2];
X = mvnrnd(mu,C,N);
CHat = cov(X);
RHat = corrcoef(X);

% Autocovariance function (biased, unbiased, non-normalized, normalized)==
N = 300;
nlags = N-1;
variance = 3;
sigma = sqrt(variance);
X = sigma*randn(N,1);
figure(1)
subplot(311)
[acf1,lags]=xcov(X,nlags,'biased');
plot(lags,acf1,'k')
ylabel('Autocovariance (biased)')
axis tight
subplot(312)
[acf2,lags]=xcov(X,nlags,'unbiased');
plot(lags,acf2,'k')
ylabel('Autocovariance (unbiased)')
axis tight
subplot(313)
[acf3,lags]=xcov(X,nlags,'coeff');
plot(lags,acf3,'k')
xlabel('Lag')
ylabel('Autocorrelation')
axis tight

% Solution of linear systems of equations ==============
N  = 1000;
p  = 3;
mu = 0.5;
B  = [.8 -2 .3]';
X  = randn(N,p);
Y  = mu + X*B + randn(N,1);
X  = [ones(N,1) X];
BHat = inv(X'*X) * (X'*Y) %Via matrix inversion
BHat = X \ Y              %Via Matlab backslash operator
[Q,R] = qr(X,0);          %Via QR decomposition (R is upper triangular)
BHat = R\(Q'*Y)

% Univariate AR(p) process simulation and estimation ==============
clear all
N = 1000;
X = zeros(N,1);
a = .78;
variance = .25;
sigma = sqrt(variance);
for t=2:N
    X(t) = a*X(t-1) + sigma*randn;
end
figure(1),clf,set(gcf,'color',[1 1 1])
plot(1:N,X,'k')
xlabel('Time t','fontsize',12),ylabel('X_t','fontsize',12)
box off

% Estimation via conditional MLE (Ordinary Least-Squares - OLS) ======
Y=X(2:end);
Z=X(1:end-1);
aHat = Y\Z;

% Estimation via the Burg maximum entropy method ============
[H,sigma2Hat,K] = arburg(X,3); 
aHat = -H(2:end)

% Convolution ==========================================
% Simple convolution example: flip, shift, multiply and sum (dot product)
% Note: for large sequences, convolutions are preferably implemented as multiplication in the
% Fourier domain.
clear all, close all
N=1000;
h=ones(10,1);
h=flip(h);
m=length(h);
E=randn(N,1);
E = [zeros(m,1);E];
X=zeros(N+m,1);
j=0;
for k = m:N+m
    j=j+1;
    X(k) = h'*E(j:j+m-1);
end
X=X(m+1:end);% In this example and in Homework 1, we do not normalize X 
E=E(m+1:end);
figure(1),clf,set(gcf,'color',[1 1 1])
plot(1:N,E,'k',1:N,X,'r')

