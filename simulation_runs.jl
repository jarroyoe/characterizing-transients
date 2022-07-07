using Plots

p0 = 10
n0 = 1
deltap = 0.9
sigma = 2.67
deltan = 0.8
gamman = 1
R = 2
nstar = R/(1-deltan)-1
star = (1-deltap)*(1+sigma*nstar^2)/nstar
gammap = 1

#Simulation function
function simulations(T,p0,n0;gammap = gammap,deltap = deltap, sigma = sigma, deltan = deltan, gamman = gamman, R = R)
    N = zeros(Float64, T, 1)
    logN = N
    P = zeros(Float64, T, 1)
    logP = P
    N[1] = n0
    P[1] = p0
    # run the simulation
    for i in 1:(T - 1)
        logN[i + 1] = log(N[i] * (deltan * exp(-gamman * P[i] / (1 + sigma * N[i]^2)) + R * exp(-P[i]) / (1 + N[i])))
        logP[i + 1] = log(P[i] * (deltap + gammap * N[i] / (1 + sigma * N[i]^2)))
        N[i + 1] = exp(logN[i + 1])
        P[i + 1] = exp(logP[i + 1])
    end

    return P,N
end

#Observed long transient near resource extinction
observedMs = zeros(100)
theorMs = zeros(100)
p0 = range(1,50,length=100)

for j in 1:100
    P,N = simulations(T,p0[j],1)

    observedMs[j] = try
                        findfirst([N[k]>N[k-1] for k in 2:T])
                    catch
                        T
                    end
    theorMs[j] = log(1/p0[j])/log(deltap)
end
bigO = (observedMs./theorMs)[end]
scatter(p0,observedMs)
plot!(p0,theorMs*bigO)

observedMs = zeros(100)
theorMs = zeros(100)
epsilons = range(1e-3,1e-1,length=100)
gammap = 8
lambda1 = deltap + gammap*nstar/(1+sigma*nstar^2)
T = 10000
P = zeros(T)
N = zeros(T)

for j in 1:100
    P,N = simulations(T,epsilons[j],nstar-epsilons[j];gammap=gammap)

    observedMs[j] = try
                        findfirst([P[k]>1 for k in 1:T])
                    catch
                        T
                    end
    theorMs[j] = log(1/epsilons[j])/log(lambda1)
end
bigO = (observedMs./theorMs)[end]
scatter(epsilons,observedMs)
plot!(epsilons,theorMs*bigO)
