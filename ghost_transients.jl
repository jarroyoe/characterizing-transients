using GLM,Plots,DataFrames
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

#Ghost transient times
function plus_coexist(
    gammap=gammap,
    deltan=deltan,
    sigma=sigma,
    gamman=gamman,
    R=R,
    deltap=deltap
)
    a = gammap * (1 + sqrt(1 - 4 * sigma * (1 - deltap)^2 / gammap^2)) / (2 * (1 - deltap) * sigma)
    function f(x)
        x - log(R / ((1 + a) * (1 - deltan * exp(-gamman * x / (1 + sigma * a^2)))))
    end

    u = fzero(f, 1)
    return [a u]
end

function minus_coexist(
    gammap=gammap,
    deltan=deltan,
    sigma=sigma,
    gamman=gamman,
    R=R,
    deltap=deltap
)
    a = gammap * (1 - (sqrt(1 - 4 * sigma * (1 - deltap)^2 / gammap^2))) / (2 * (1 - deltap) * sigma)
    function f(x)
        x - log(R / ((1 + a) * (1 - deltan * exp(-gamman * x / (1 + sigma * a^2)))))
    end

    u = fzero(f, 1)
    return [a u]
end

function transient_time(gammap=0.9 * star; deltan=deltan, sigma=sigma, gamman=gamman, R=R, deltap=deltap, T=T)
    # define initial condition
    lowerBound = minus_coexist(gammap,
    deltan,
    sigma,
    gamman,
    R,
    deltap)
    upperBound = plus_coexist(gammap,
    deltan,
    sigma,
    gamman,
    R,
    deltap)
    
    N0 = 1.01 * lowerBound[1]
    P0 = 0.99 * lowerBound[2]

    N = zeros(Float64, T, 1)
    logN = N
    P = zeros(Float64, T, 1)
    logP = P
    N[1] = N0
    P[1] = P0
    # run the simulation
    for i in 1:(T - 1)
        logN[i + 1] = log(N[i] * (deltan * exp(-gamman * P[i] / (1 + sigma * N[i]^2)) + R * exp(-P[i]) / (1 + N[i])))
        logP[i + 1] = log(P[i] * (deltap + gammap * N[i] / (1 + sigma * N[i]^2)))
        N[i + 1] = exp(logN[i + 1])
        P[i + 1] = exp(logP[i + 1])
    end

    return findfirst(>(upperBound[1]), N)[1]
end

epsilons = [0.0088,0.009,0.01,0.04,0.06,0.08,0.1]
times = transient_time.((1 .- epsilons) .* star)
linearRegression = lm(@formula(logtimes ~ logepsilons), 
    DataFrame(logtimes=log.(times), logepsilons=log.(epsilons)))
coefficients=DataFrame(A=exp(coef(linearRegression)[1]),
    B=coef(linearRegression)[2])
scatter(epsilons,times,legend=false,label="",   
    yguidefontsize=12,
    xguidefontsize=12,
    ytickfontsize=12,
    xtickfontsize=12)

x_epsilons = 1e-4:1e-4:0.1
plot!(x_epsilons,coefficients.A[1].*x_epsilons.^coefficients.B[1],
    legend=true,
    label=string("τ=",round(coefficients.A[1],digits=2),
                "ε^",round(coefficients.B[1],digits=2)),
    legendfontsize=12)
xlabel!("ε")
ylabel!("τ")