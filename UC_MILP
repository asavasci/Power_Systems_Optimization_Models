using JuMP, Gurobi

#=
Reference Paper:
R. A. Jabr, “Tight polyhedral approximation for mixed-integer linear programming
unit commitment formulations,” Transmission Distribution IET Generation,
vol. 6, no. 11, pp. 1104–1111, Nov. 2012.
=#


## D_t : Net load in period t ∈ T, ( MW )
#
D = [700,750,850,950,1000,1100,1150,1200,1300,1400,1450,1500,
     1400,1300,1200,1050,1000,1100,1200,1400,1300,1100,900,800];

## Cost function parameters for 10-unit system
#
a = [1000,970,700,680,450,370,480,660,665,670]; # Fixed cost of running generator g ∈ G, ( $/h )
b = [16.19,17.26,16.60,16.50,19.70,22.26,27.74,25.92,27.27,27.79]; # ( $/MWh)
c = [0.00048,0.00031,0.002,0.00211,0.00398,0.00712,0.00079,0.00413,0.00222,0.00173]; # ( $/MW^2h)

## Maximum and Minimum production limits of generators
#
Pmax = [682.5,682.5,195,195,243,120,127.5,82.5,82.5,82.5]; # ( MW )
Pmin = [225,225,30,30,37.5,30,37.5,15,15,15]; # ( MW )

## Start up and Ramp up rate of generators
#
StartUpRate = [337.5,337.5,45,45,56.25,45,56.25,22.5,22.5,22.5]; # ( MW )
RampUpRate  = [405,405,54,54,67.5,54,67.5,27,27,27]; # ( MW )

## Shut down and Ramp down rate of generators
#
ShutDwnRate = [337.5,337.5,45,45,56.25,45,56.25,22.5,22.5,22.5]; # ( MW )
RampDwnRate = [405,405,54,54,67.5,54,67.5,27,27,27];  # ( MW )

## Minimum up and Minimum down time of generators
#
UT = [8,8,5,5,6,3,3,1,1,1]; # (h)
DT = [8,8,5,5,6,3,3,1,1,1]; # (h)

## Start-up cost of generator
#
StartUpCost = [4500,5000,550,560,900,170,260,30,30,30]; # ( $/h )
ShutDwnCost = zeros(10); # ( $/h )

# Reserve
rt = 0.1

T = 24; # Number of time periods
G = 10; # Number of generators


T =  length(D); # Number of time slots:
NG = length(Pmax) # Number of generators
NL = 4  # Number of segments in piecewise linear approximation

uc = Model(Gurobi.Optimizer) # UC model


# defining variables
@variable(uc, P[1:NG,1:T] >= 0)
@variable(uc, delta[1:NL, 1:NG, 1:T] >= 0)
@variable(uc, U[1:NG,1:T], Bin) # on/off variable
@variable(uc, V[1:NG,1:T], Bin) # Start up decision
@variable(uc, W[1:NG,1:T], Bin) # Shut dwn decision

@constraint(uc, Power_Balance[t=1:T], sum(P[g,t] for g in 1:NG) == D[t] + 0.1D[t] )


@constraint(uc, GenUpperBound[t=1:T,g=1:NG], Pmin[g]*U[g,t] <= P[g,t] )
@constraint(uc, GenLowerBound[t=1:T,g=1:NG], P[g,t] <= Pmax[g]*U[g,t] )

@constraint(uc, Rampup[t=2:T,g=1:NG], P[g,t] - P[g,t-1] <= StartUpRate[g]*V[g,t] + RampUpRate[g]*U[g,t-1]  )
@constraint(uc, Rampdw[t=2:T,g=1:NG], P[g,t-1] - P[g,t] <= ShutDwnRate[g]*W[g,t] + RampDwnRate[g]*U[g,t]  )

@constraint(uc, [g=1:NG, t=1:T], P[g,t] == sum( delta[l,g,t] for l=1:NL) + Pmin[g]*U[g,t] )

@constraint(uc, UpTime[g=1:NG,t=UT[g]:T], sum(V[g,k] for k=t-UT[g]+1:t ) <= U[g,t] )
@constraint(uc, DwTime[g=1:NG,t=DT[g]:T], sum(W[g,k] for k=t-DT[g]+1:t ) <= 1 - U[g,t] )

@constraint(uc, Logic[t=2:T, g=1:NG], U[g,t-1] - U[g,t] + V[g,t] - W[g,t] == 0)

# --- Piece-wise linearization

x = zeros(NL+1,NG)
slope = zeros(NL,NG)

for g in 1:NG

    x[:,g] = range(Pmin[g],stop=Pmax[g],length=NL+1);  # liner segments start - end point
    slope[:,g] = c[g] * (x[2:end,g] + x[1:end-1,g]) + b[g]*ones(NL,1)  # slopes for each linear segment

end

@constraint(uc, d_1[g=1:NG, t=1:T], delta[1,g,t] <= x[2,g] - Pmin[g] )
@constraint(uc, d_1N[l=2:NL-1, g=1:NG, t=1:T], delta[l,g,t] <= x[l,g] - x[l-1,g] )
@constraint(uc, d_N[g=1:NG, t=1:T], delta[NL,g,t] <= Pmax[g] - x[NL,g] )

# ---

# Approximated Cost Function
A(g) = a[g]+b[g]*Pmin[g]+c[g]*Pmin[g]^2;

@objective(uc, Min, sum([A(g)*U[g,t] + sum( delta[l,g,t]*slope[l,g] for l=1:NL ) + StartUpCost[g]*V[g,t] + ShutDwnCost[g]*W[g,t] for t=1:T for g=1:NG]) )


optimize!(uc)

# --- Problem Formulation Ends

println("\n\n")

# Piecewise Cost Value
println("Piecewise Cost = ", JuMP.objective_value(uc))

# Quadratic Cost Value
qc = sum( value(U[g,t])*a[g] + b[g]*value(P[g,t]) + c[g]*value(P[g,t])^2+ StartUpCost[g]*value(V[g,t]) + ShutDwnCost[g]*value(W[g,t])  for g = 1:NG for t=1:T )
println("Quadratic Cost = ", qc)
