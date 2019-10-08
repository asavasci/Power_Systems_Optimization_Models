clearconsole()

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
MinUpTime = [8,8,5,5,6,3,3,1,1,1]; # (h)
MinDwnTime = [8,8,5,5,6,3,3,1,1,1]; # (h)

## Start-up cost of generator
#
StartUpCost = [4500,5000,550,560,900,170,260,30,30,30]; # ( $/h )
ShutDwnCost = zeros(10); # ( $/h )



T = 24; # Number of time periods
G = 10; # Number of generators

using JuMP
using CPLEX


## Model and solver decleration
#
optimizer = CPLEX.Optimizer
uc = Model(with_optimizer(optimizer))


## Problem variables
#
@variable(uc, P[1:G,1:T]) # Production amount of generator g ∈ G in period t ∈ T
@variable(uc, u[1:G,1:T], Bin) # Status of generator g ∈ G in period t ∈ T -- 1 if generator g is ON in period t; 0 otherwise
@variable(uc, y[1:G,1:T], Bin) # Start up decision of generator g ∈ G in period t ∈ T  -- 1 if u[g,(t−1)] = 0 and u[g,t] = 1; 0 otherwise
@variable(uc, z[1:G,1:T], Bin) # Shut down decision of generator g ∈ G in period t ∈ T -- 1 if u[g,(t−1)] = 1 and uit = 0; 0 otherwise

## Problem constraints
#
@constraint(uc, Power_Balance_Constraint[t=1:T],
                sum(P[g,t] for g=1:G ) == D[t] )


@constraint(uc, GenerationLimitsMax_Constraint[t=1:T,g=1:G],
                P[g,t] <= Pmax[g]*u[g,t] )


@constraint(uc, GenerationLimitsMin_Constraint[t=1:T,g=1:G],
                Pmin[g]*u[g,t] <= P[g,t] )


@constraint(uc, MinUpTime_Constraint[g=1:G, t=2:T, tau=t+1:min(t+MinUpTime[g],T)],
                u[g,t] - u[g,t-1] <= u[g,tau])

@constraint(uc, MinDwnTime_Constraint[g=1:G, t=2:T, tau=t+1:min(t+MinDwnTime[g],T)],
                u[g,t-1] - u[g,t] <= 1 - u[g,tau])

@constraint(uc, StartUp_Constraint[g=1:G, t=2:T],
                u[g,t] - u[g,t-1] <= y[g,t])

@constraint(uc, ShutDwn_Constraint[g=1:G, t=2:T],
                u[g,t-1] - u[g,t] <= z[g,t])

@constraint(uc, RampStartUp_Constraint[g=1:G, t=2:T],
                P[g,t] - P[g,t-1] <= StartUpRate[g]*y[g,t] + RampUpRate[g]*u[g,t-1])

@constraint(uc, RampShutDwn_Constraint[g=1:G, t=2:T],
                P[g,t-1] - P[g,t] <= ShutDwnRate[g]*z[g,t] + RampDwnRate[g]*u[g,t])



## Objective Function
#
@objective(uc, Min, sum(a[g]*u[g,t] + b[g]*P[g,t] + c[g]*P[g,t]*P[g,t]
                    + StartUpCost[g]*y[g,t] + ShutDwnCost[g]*z[g,t] for g=1:G, t=1:T ))


optimize!(uc)
