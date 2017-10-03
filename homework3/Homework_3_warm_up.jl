# Simulate a 5-area (~20 Hz) beta network instantiated by an AR(3) model =============================
# cd("Homework_3_Lect_5&6")

xDim=5 # state dimensionå
order=3# p=3 order of the model
T=3# Trial duration in seconds
nTrials=200
Fs=200# sampling rate in per second
dt=1/Fs
burnin=1000 # for transient removal
df=0.001#frequency resolution to be used in the parametric AR estimation
nF=1/df+1#number of frequencies evaluated in the parametric AR estimation

N=T*Fs

# The AR(3) matrices 
A = zeros(xDim,xDim,order)
a=sqrt(2)
A[1,1,1] =  0.95*a
A[1,1,2] = -0.9025
A[2,1,2] =  0.5
A[3,1,3] = -0.4
A[4,1,2] = -0.5
A[4,4,1] =  0.25*a
A[4,5,1] =  0.25*a
A[5,4,1] = -0.25*a
A[5,5,1] =  0.25*a

# EXERCISE (1a): Stability/Stationarity Build the companion matrix (augmented state AR(1)) and check stability
m,n,p = size(A)
pn = (p-1)*m
Ac = [reshape(A,m,p*n); eye(pn) zeros(pn,m)]# companion matrix



# Sampling from the VAR(3) model
# The noise covariance (Using DIAGONAL COVARIANCE MATRIX FOR SIMPLICITY ...)
SIGMA =diag([0.60.5 0.3 0.3 0.6])
mu = zeros(xDim,1)
for r=1:nTrials
    r
    E=mvnrnd(mu,SIGMA,N+burnin)" #Multivariate GWN sequence
    X(:,1:order,r)=E(:,1:order)    
    for k=order+1:N+burnin
        x=0
        for j=1:order
            x=x+squeeze(A(:,:,j))*squeeze(X(:,k-j,r))
        end
        X(:,k,r)=x+E(:,k)
    end
end
X=X(:,burnin+1:end,:)

# Computing the spectral matrix for the sampled VAR(p) data via multitaper ===========
NFFT=N
bandWidth = ...
removeTemporalMean=true
RemoveEnsembleMean=true
nTapers=[]
clear S
xDim = 5
for j=1:xDim
    x=squeeze(X(j,:,:)) 
    for k=j+1:xDim
        y=squeeze(X(k,:,:))        
        [Sxx, Syy, Sxy, ~, ~, F, nTapers]= multitaperSpectrum(x,y,Fs,bandWidth,NFFT,removeTemporalMean,RemoveEnsembleMean,nTapers)                 
        S(j,j,:)=Sxx        
        S(j,k,:)=Sxy #Cross-spectrum
        S(k,j,:)=conj(Sxy)         
    end    
end
S(j,j,:)=Syy

# To compute the partial coherence =======================
# S is the spectral matrix S(channel j, channel k, frequency), j,k = 1, 2, ... xDim
for j=1:xDim
    for k=j+1:xDim
        i1=setdiff([1:xDim],j)
        i2=setdiff([1:xDim],k)
        for f=1:length(F)
            Mjk = det(S(i1,i2,f))
            Mjj = det(S(i1,i1,f))
            Mkk = det(S(i2,i2,f))            
            Cp(j,k,f) = abs(Mjk)/real(sqrt(Mjj*Mkk))#to avoid numerical issues, force it to be real 
        end                
    end
end



# Estimating VAR(p) from data ===================
[Ah,SIGMAh,Eh] = var_maxent(X,order)
DSIG = det(SIGMA)# residuals covariance matrix determinant
if DSIG <= 0
    fprintf(2,"  WARNING: residuals covariance not positive definite\n")
end
M=N*nTrials
L= -(M/2)*log(DSIG) #max loglikelihood
aic = -2*L + 2*order*xDim^2*(M/(M-order-1)) # Note AIC without correction = -2*L + 2*order*xDim^2
bic = -2*L + order*log(M)

# Transfer function H(f) and spectral matrix S(f) from fitted model
# F=[0:nF-1]*df*Fs/2
delta=1
F=linspace(0,0.5,nF)# use normalized frequency (cycles/sample) first 
j=0
for f=F
    j=j+1
    H = eye(xDim) # identity matrix
    for m=1:order
        H=H-squeeze(Ah(:,:,m))*exp(-1i*m*2*pi*f)        
    end
    H = inv(H)
    #S(:,:,j) = H*SIGMAh*ctranspose(H) 
    S(:,:,j) = H*SIGMAh*H"
    #Note: the transpose " " " will also recognize the complex data and perform the conjugate transpose
    #Note also, however, that ." implements only the transpose, not the conjugate transpose
end
S = 2 * delta^2 * 1/(N*delta) * S #Normalize the spectral matrix
for k = 1:xDim, S(k,k,:) = real(S(k,k,:))end #make sure the diagonal has real numbers (power spectrum)
F=F*Fs# Change frequency to Hz



