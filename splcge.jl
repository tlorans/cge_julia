# A Simple CGE Model in Ch. 5 (SPLCGE,SEQ=275)

# Load Packages
using Ipopt, JuMP, DataFrames, XLSX, GLM, Statistics

#data = DataFrame(XLSX.readtable("/Users/thomaslorans/Documents/CGE/Textbook/SAM_splcge.xlsx", "SAM")...)

# Creat JuMP Model
model = Model(Ipopt.Optimizer)


## Define indices

name_goods = ["BRD","MLK"] # Dimension of goods
name_factors = ["CAP","LAB"]    # Dimension of factors
goods = [1,2]
factors = [1,2]
# SAM
sam_content = [missing     missing     missing  missing     15
               missing     missing     missing  missing     35
               5.0         20.0        missing  missing     missing
               10.0        15.0        missing  missing     missing
               missing     missing     25       25          missing]

names = ["BRD", "MLK", "CAP", "LAB", "HOH"]

SAM = DataFrame(sam_content)
rename!(SAM, names)


####################################### STEP 1 PARAMETERS ###########################

X0 = Array{Float64}(SAM[1:2,:HOH])     # household consumption of the i-th good
F0 = Array{Float64}(SAM[3:4, name_goods]) # the h-th factor input by the j-th firm
Z0 =  sum(F0, dims = 1) # output of the j-th good
FF = Array{Float64}(SAM[5, name_factors]) # factor endowment of the h-th factor

alpha =  X0 / sum(X0) # share parameter in utility function
beta  =  F0 ./ sum(F0, dims = 1)  # share parameter in production function
b = [Z0[i] / prod(F0[:,i].^beta[:,i]) for i in goods] # scale parameter in production function

px0 = [1, 1] # demand price of the i-th good
pz0 = [1, 1] # supply price of the i-th good
pf0 = [1, 1] # the h-th factor price


######################## STEP 2 VARIABLES ######################
@variables model begin
  X[i = goods], (start = X0[i]) # household consumption of the i-th good
  F[j = factors, i=goods], (start = F0[j,i]) # the h-th factor input by the j-th firm
  Z[i = goods], (start = Z0[i]) # output of the j-th good
  px[i = goods], (start = px0[i]) # demand price of the i-th good
  pz[i = goods], (start = pz0[i]) # supply price of the i-th good
  pf[j = factors], (start = pf0[j]) # the h-th factor price
  UU # objective function
end


############### STEP 3 EQUATIONS/CONSTRAINTS ########################

@NLconstraints model begin
  EQ_X[i = goods], X[i] == alpha[i] * sum(pf[h] * FF[h] for h in factors) / px[i] # household demand function
  EQ_pz[i = goods], Z[i] == b[i] * prod(F[j,i]^beta[j,i] for j in factors) # production function
  EQ_F[j = factors, i = goods], F[j,i] == beta[j,i] * pz[i] * Z[i] / pf[j] # factor demand function
  EQ_px[i = goods], X[i] == Z[i] # good market clearing condition
  EQ_pf[j = factors], sum(F[j,i] for i in goods) == FF[j] # factor market clearing condition
  EQ_Z[i = goods], px[i] == pz[i] # price equation
  UU == prod(X[i]^alpha[i] for i in goods)
end

# Price of labor as numeraire
JuMP.fix(pf[2], 1; force=true);

# setting lower bounds to avoid division per zero
for i in goods
  JuMP.set_lower_bound(X[i], 0.001)
  JuMP.set_lower_bound(Z[i], 0.001)
  JuMP.set_lower_bound(pf[1], 0.001)
  JuMP.set_lower_bound(px[i], 0.001)
  JuMP.set_lower_bound(pz[i], 0.001)
end

# Model Solver

@NLobjective(model, Max, UU)

@time optimize!(model)

termination_status(model)
primal_status.(model)
result_count(model)
value.(X)
value.(alpha)
value.(pf)
value.(F)
value.(Z)
value.(px)
value.(pf)
objective_value(model)
