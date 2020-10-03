# A Household's Utility Max. Model in Ch. 3 (HHMAX,SEQ=274)

# Load Packages
using Ipopt, JuMP, DataFrames


# Creat JuMP Model
model = Model(Ipopt.Optimizer)

## Define indices

goods = [1,2] # Dimension of goods
factors = [1,2]    # Dimension of factors


####################################### STEP 1 PARAMETERS ###########################

alpha = [0.2, 0.8] # share parameter in utility function
px = [1, 2] # price of the i-th good
pf = [2, 1] # price of the h-th factor
FF = [10, 20] # Factor endowment

######################## STEP 2 VARIABLES ######################

@variables model begin
  X[i = goods]
  UU
end


############### STEP 3 EQUATIONS/CONSTRAINTS #

@NLconstraints model begin
  EQ_X[i = goods], X[i] == alpha[i] * sum(pf[h] * FF[h] for h in factors) / px[i] # household demand function
  EQ_UU, UU == prod(X[i] * alpha[i] for i in goods) # utility function
end

# Model Solver

@NLobjective(model, Max, 0)
@time optimize!(model)

termination_status(model)
primal_status.(model)
result_count(model)
value(UU)
value.(X)
