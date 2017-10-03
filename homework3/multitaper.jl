# FrequencyDomainAnalysis.jl
# Tools for spectral density estimation and analysis of phase relationships
# between sets of signals.

# Copyright (C) 2013   Simon Kornblith

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

export multitaper, psd, xspec, coherence

# Perform tapered FFT of continuous signals
# A is samples x channels x trials
function multitaper(A::Union{AbstractVector{T}, AbstractMatrix{T}, AbstractArray{T,3}},
                             stats:: Vector{S}, fs::Real=1;
                             tapers::Union{Vector{T}, Matrix{T}}=dpss(size(A, 1), 4), nfft::Int=nextpow2(size(A, 1)),
                             freqrange::Range{Int}=0:-1) where {T<:Real, S<:TransformStatistic}

    # We can't guarantee order of kwargs, so this has to be here
    if freqrange == 0:-1
        freqrange = 1:(nfft >> 1 + 1)
    end

    n = size(A, 1)
    nout = nfft >> 1 + 1
    ntapers = size(tapers, 2)
    nchannels = size(A, 2)

    multiplier = convert(T, sqrt(2/fs))
    firstsamplemultiplier = freqrange[1] == 1 ? convert(T, sqrt(1/fs)) : multiplier
    lastsamplemultiplier = freqrange[end] == nout && iseven(nfft) ?
                           convert(T, sqrt(1/fs)) : multiplier

    for stat in stats
        init(stat, length(freqrange), nchannels, ntapers)
    end

    # dtype = outputtype(T)
    dtype = Float64
    fftin = zeros(Complex{dtype}, nfft, size(A, 2))
    fftout = Array{Complex{dtype}}(nout, size(A, 2))
    #fftview = sub(fftout, freqrange, 1:size(A,2))
    fftview = zeros(Complex{dtype}, length(freqrange), size(A, 2))

    # p = FFTW.plan_fft(fftin, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    for j = 1:size(A, 3), i = 1:ntapers
        @inbounds for k = 1:nchannels, l = 1:n
            fftin[l, k] = A[l, k, j]*tapers[l, i]
        end

        # FFTW.execute(dtype, p.plan)
        fftout=fft(fftin, 1)
        # Copy FFT from one buffer to another, including only
        # frequencies in freqrange
        @inbounds for k = 1:size(fftout, 2)
            # Correct edges of FFT so that all power is in the same
            # units. This is necessary because the FFT divides power
            # associated with real signals between positive and
            # negative frequencies, except for power at 0 and FMAX. We
            # don't compute the negative frequencies, but we do need
            # to scale the FFT so that the power is correct.
            fftview[1, k] = fftout[freqrange[1], k]*firstsamplemultiplier
            for l = 2:size(fftview, 1)-1
                fftview[l, k] = fftout[freqrange[l], k]*multiplier
            end
            fftview[end, k] = fftout[freqrange[end], k]*lastsamplemultiplier
        end

        for stat in stats
            accum(stat, fftview, i)
        end
    end

    [finish(stat) for stat in stats]
end

# Calling with a single type or instance
multitaper(A::Union{AbstractVector{T}, AbstractMatrix{T}, AbstractArray{T,3}},
                    stat::S, args...; kw...) where {T<:Real, S <: TransformStatistic} =
    multitaper(A, [stat] , args...; kw...)[1]

# Convenience functions
psd{T<:Real}(A::AbstractArray{T}, args...; kw...) =
    multitaper(A, (PowerSpectrum{T}(),), args...; kw...)[1]
xspec{T<:Real}(A::Vector{T}, B::Vector{T}, args...; kw...) =
    multitaper(hcat(A, B), (CrossSpectrum{T}(),), args...; kw...)[1]
coherence{T<:Real}(A::Vector{T}, B::Vector{T}, args...; kw...) =
    multitaper(hcat(A, B), (Coherence{T}(),), args...; kw...)[1]


# using DSP
# using FFTW
# include("transform_stats.jl")
# include("multitaper.jl")
# x = 0:63
# signal = zeros(length(x), 1, 35)
# for i = 1:size(signal, 3)
# 	signal[:, i] = cos.(0.2pi*x+rand()*2pi)
# end
# c = multitaper([signal signal], Coherence())
