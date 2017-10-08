"""
    function var_maxent(X,p)
Fits Vector AR model (VAR) to multi-trial, multivariate time series data
NOTE: May be inneficient!!! Direct translation from MATLAB code

### Arguments

* X - multi-trial time series data 
      [n,m,N] = size(X); for n = number of variables (e.g.
       channels), m = number of sample per trial, N = number of trials.
* p  - model order (number of lags)

### output

* A          VAR coefficients matrix
* SIGMA      residuals covariance matrix
* E          residuals time series

Returns VAR coefficients |A| and (optionally) residuals covariance matrix
|SIG| and serially uncorrelated residuals |E| for the |p|-lag autoregression

(where  [[ii_Sigma.png]] = |SIG|) of a stationary multivariate process
|X|. |X| may contain single- or multi-trial multivariate time series
data. LWR algorithm [1].

[1] M. Morf, A. Viera, D. T. L. Lee and T. Kailath, "Recursive Multichannel
Maximum Entropy Spectral Estimation", _IEEE Trans. Geosci. Elec._, 16(2), 1978.
This function is based on the MVGC toolbox (Barnett and Seth, 2014)
"""
function var_maxent(X,p)

    n,m,N = size(X);
    assert(p < m);
    p1 = p+1;

    A   = NaN; # assure a "bad" return value if anything goes wrong (see routine "isbad")
    SIGMA = NaN; # assure a "bad" return value if anything goes wrong (see routine "isbad")
    E   = NaN; # assure a "bad" return value if anything goes wrong (see routine "isbad")

    X = demean_time(X); # no constant term

    q1n = p1*n;

    I = eye(n);

    # store lags
    XX = zeros(n,p1,m+p,N);
    for k = 0:p
        XX[:,k+1,k+1:k+m,:] = X; # k-lagged observations
    end

    # initialise recursion
    AF = zeros(n,q1n); # forward  AR coefficients
    AB = zeros(n,q1n); # backward AR coefficients (reversed compared with Morf"s treatment)

    k  = 1;            # model order is k-1
    kn = k*n;
    M  = N*(m-k);
    kf = 1:kn;         # forward  indices
    kb = q1n-kn+1:q1n; # backward indices

    XF = reshape(XX[:,1:k,k+1:m,:],kn,M);
    XB = reshape(XX[:,1:k,k:m-1,:],kn,M);

    CXF = chol(XF*XF');
    # if cholp, return; end # show-stopper! MATLAB check for A positive-definite. Julia? 

    CXB= chol(XB*XB');
    # if cholp, return; end # show-stopper! MATLAB check for A positive-definite. Julia? 


    AF[:,kf] = CXF'\I;
    AB[:,kb] = CXB'\I;

    AFPREV = zeros(size(AF))
    EF = zeros(kn, M)

    # and loop
    while k <= p

        EF = AF[:,kf]*reshape(XX[:,1:k,k+1:m,:],kn,M); # forward  prediction errors
        EB = AB[:,kb]*reshape(XX[:,1:k,k:m-1,:],kn,M); # backward prediction errors

        CEF = chol(EF*EF');
    #     if cholp, return; end  # show-stopper!

        CEB = chol(EB*EB');
    #     if cholp, return; end  # show-stopper!

        R = CEF'\(EF*EB')/CEB; # normalised reflection coefficients

        RF = chol(I-R*R');
    #     if cholp, return; end  # show-stopper!

        RB = chol(I-R'*R);
    #     if cholp, return; end  # show-stopper!

        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
        kf = 1:kn;
        kb = q1n-kn+1:q1n;

        AFPREV = AF[:,kf];
        ABPREV = AB[:,kb];

        AF[:,kf] = RF'\(AFPREV-R*ABPREV);
        AB[:,kb] = RB'\(ABPREV-R'*AFPREV);

    end

    E   = AFPREV[:,1:n]\EF;   # residuals
    SIGMA = (E*E')/(M-1);       # residuals covariance matrix
    E   = reshape(E,n,m-p,N); # put residuals back into per-trial form
    A = reshape(-AF[:,1:n]\AF[:,n+1:end],n,n,p); # so A(:,:,k) is the k-lag coefficients matrix
   
    return A,SIGMA,E

end

"""
    function demean(X; normalise=false)
    
Temporal demean of time series data

### Arguments

* X - multi-trial time series data
* normalise::Bool  - normalise (temporal) variance of each variable to 1 (default: false)

### output

* Y - demeaned time series data
"""
function demean_time(X::Array{Float64,3}; normalise=false)
    n,m,N = size(X);
    U = ones(1,N*m);
    Y = X.-mean(X,2);
    if normalise
        Y = Y./std(Y,2);
    end
    return Y
end


