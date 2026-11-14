#***********************************************************************************************
# Numerical simulation of a fully coupled magneto-mechanical model for Dowsil-based hMREs, 
# incorporating explicit representations of the magnet and the surrounding air.
#
# Solved using an Staggered strategy:
# -> Problem 1 magnetic problem solid+vacuum
# -> Problem 2 mechanical problem (solid movement)
# -> Problem 3 mechanical problem (vacuum movement)
#
# Topology and material optimization in ultra-soft magnetoactive structures: making advantage of residual anisotropies
# Carlos Perez-Garcia, Rogelio Ortigosa, Jesus Martinez-Frutos Daniel Garcia-Gonzalez
#***********************************************************************************************

# Load required Packages
using HyperFEM
using Gridap, GridapGmsh, GridapSolvers, DrWatson
using Gridap.FESpaces
using Gridap.CellData

# Initialize problem
pname = "Vacuum_with_magnet"
meshfile = "Short_magnet.msh"

simdir = datadir("sims", pname)
setupfolder(simdir)

# Load mesh file
geomodel = GmshDiscreteModel(datadir("models", meshfile))

#********************
# Constitutive models
#********************

# Magnet
μ0=4π*1e-7
αr_max = 0.165 / μ0
model_magnet = Magnetic(μ0=μ0, αr=αr_max, χe=0.05)

# Magnetic model vacuum
model_vacuum_mag   = IdealMagnetic2D(μ0=μ0, χe=0.0)

# Magneto-Mechanical model for Solid hMREs movement
params      =  [10.0, 12480.64495232286, 1.0, 2.0, 5195.545287237134, 0.2602717127043121,  2.0]
modelmech   =  NonlinearMooneyRivlin2D(λ=(params[1] + params[2]) * 1e2, μ1=params[1], μ2=params[2], α1=params[3], α2=params[4]) + 
               TransverseIsotropy2D(μ=params[5], α1=params[6], α2=params[7])
modelhard   =  HardMagnetic2D(μ0=μ0, αr=10e-3/μ0, χe=0.0, χr=8.0)
model_solid =  MagnetoMechModel(magneto=modelhard, mechano=modelmech)


# Mechanical model for Vacuum movement
model_vacuum_mech_ = NonlinearMooneyRivlin2D_CV(λ=1.0, μ1=1.0, μ2=0.0, α1=6.0, α2=1.0, γ=6.0)
model_vacuum_mech  = HessianRegularization(mechano=model_vacuum_mech_, δ=1e-6)

#**********************
# Finite Elements
#**********************
order_solid  = 1       # Linear FEs for solid
order_vacuum = 1       # Linear FEs for vacuum
reffeφ        = ReferenceFE(lagrangian, Float64, order_solid)
reffeu        = ReferenceFE(lagrangian, VectorValue{2,Float64}, order_solid)
reffeu_vacuum = ReferenceFE(lagrangian, VectorValue{2,Float64}, order_vacuum)

#**********************
# Domains
#**********************
Ωmagnet     = Interior(geomodel, tags=["Magnet"])
Ωdomain     = Triangulation(geomodel, tags=["Vacuum", "Beam", "Magnet"])
Ωvacuum     = Triangulation(geomodel, tags=["Vacuum"])
Ωsolid      = Triangulation(geomodel, tags=["Beam"])

dΩdomain      = Measure(Ωdomain, 2 * order_solid)
dΩsolid       = Measure(Ωsolid,  2 * order_solid)
dΩvacuum_mag  = Measure(Ωvacuum, 2 * order_solid)
dΩvacuum_mec  = Measure(Ωvacuum, 2 * order_vacuum)
dΩmagnet      = Measure(Ωmagnet, 2 * order_solid)

# Solid <-> Vacuum interface
Γsolid_int = BoundaryTriangulation(Ωsolid, tags="Interface")
dΓs = Measure(Γsolid_int, 4 * order_solid)

#************************************
# Dirichlet boundary conditions
#************************************

# Problem 1: Dirichlet conditions magnetic field
evolφ(Λ) = Λ
Dφ = BoundaryConditions(BoundaryCondition("phi_Dirichlet", 0.0, evolφ))

# Problem 2: Dirichlet conditions solid
evolu(Λ) = 1.0
Du_solid = BoundaryConditions(
  BoundaryCondition("Fixed_Solid", [0.0, 0.0], evolu),
  BoundaryCondition("Fixed_Solid_x", [0.0, 0.0], evolu; mask=(true, false)))

# Problem 3: Dirichlet conditions vacuum movement: below, with the solid–vacuum coupling

#******************************************
# Finite Element Spaces and state variables
# (⁺) current stage
# (⁻) previous stage
#******************************************

# Problem 1: FE Spaces and state variables
Vφ = TestFESpace(Ωdomain, reffeφ, Dφ, conformity=:H1)
Uφ⁺ = TrialFESpace(Vφ, Dφ)
Uφ⁻ = TrialFESpace(Vφ, Dφ)
φh⁺ = FEFunction(Uφ⁺, zero_free_values(Uφ⁺))
φh⁻ = FEFunction(Uφ⁻, zero_free_values(Uφ⁻))

# Problem 2: FE Spaces and state variables
Vu_solid = TestFESpace(Ωsolid, reffeu, Du_solid, conformity=:H1)
Uu_solid⁺ = TrialFESpace(Vu_solid, Du_solid)
Uu_solid⁻ = TrialFESpace(Vu_solid, Du_solid)
uh_solid⁺ = FEFunction(Uu_solid⁺, zero_free_values(Uu_solid⁺))
uh_solid⁻ = FEFunction(Uu_solid⁻, zero_free_values(Uu_solid⁻))

# ******************************************************
# Solid–vacuum coupling
# State of the solid between two staggered stages, prescribed on the
# interface of the vacuum. HyperFEM locates the nodes of the interface
# in the mesh of the solid once (they coincide: a copy of dof values)
#******************************************************

# Interface of the solid (space of the Maxwell traction)
Usolid_int = FESpace(Γsolid_int, reffeu)

uhsolid = StageBlend(uh_solid⁻, uh_solid⁺)      # uhsolid(Λ) = uh_solid⁻ + (uh_solid⁺ - uh_solid⁻) * Λ

# Problem 3: Dirichlet conditions vacuum movement (the interface follows the solid: uhsolid(Λ))
Du_vacuum = BoundaryConditions(
  BoundaryCondition("Fixed_Air", [0.0, 0.0], evolu),
  BoundaryCondition("Fixed_Air_x", [0.0, 0.0], evolu; mask=(true, false)),
  BoundaryCondition("Interface", Coupling(uhsolid)))

# Problem 3: FE Spaces and state variables 
Vu_vacuum = TestFESpace(Ωvacuum, reffeu_vacuum, Du_vacuum, conformity=:H1)
Uu_vacuum⁺ = TrialFESpace(Vu_vacuum, Du_vacuum)
Uu_vacuum⁻ = TrialFESpace(Vu_vacuum, Du_vacuum)
uh_vacuum⁺ = FEFunction(Uu_vacuum⁺, zero_free_values(Uu_vacuum⁺))
uh_vacuum⁻ = FEFunction(Uu_vacuum⁻, zero_free_values(Uu_vacuum⁻))

# ******************************************************
#    Maxwell stress in vacuum & traction force in interface 
#******************************************************

V_maxwell = TestFESpace(Ωvacuum, reffeu_vacuum)
Pmaxw = FEFunction(V_maxwell, zero_free_values(V_maxwell))
Pmaxw_ = get_free_dof_values(Pmaxw)

V_traction = TestFESpace(Γsolid_int, reffeu)
U_traction = TrialFESpace(V_traction)
tmaxw = interpolate_everywhere(Pmaxw, Usolid_int)
tmaxw_ = get_free_dof_values(tmaxw)

T_maxwell = Transfer(V_maxwell, Usolid_int)      # vacuum → interface of the solid, built once
Massoperator = StaticLinearModel((u, v) -> ∫(u ⋅ v)dΓs, U_traction, V_traction, NothingBC())
 
# ******************************************************
#               Weak Forms
#******************************************************

# Derivatives of the energy functions
Ψs, ∂Ψs∂F, ∂Ψs∂H0, ∂Ψs∂FF, ∂Ψs∂H0F, ∂Ψs∂H0H0 = model_solid()          
Ψv, ∂Ψv∂F, ∂Ψv∂H0, ∂Ψv∂FF, ∂Ψv∂H0F, ∂Ψv∂H0H0 = model_vacuum_mag()   
Ψvm, ∂Ψvm∂F, ∂Ψvm∂FF = model_vacuum_mech()  

# Kinematic functions
F, H, J = get_Kinematics(Kinematics(Mechano, Solid))
ℋ₀     = get_Kinematics(Kinematics(Magneto, Solid))

# Staggered evolucion of state variable
uhvacuum(Λ) = uh_vacuum⁻ + (uh_vacuum⁺ - uh_vacuum⁻) * Λ
φh(Λ) = φh⁻ + (φh⁺ - φh⁻) * Λ

# Magnetization field in beam
V_Nbeam = TestFESpace(Ωsolid, reffeu)
Nhbeam  = interpolate_everywhere((x)->VectorValue(1.0, 0.0), V_Nbeam)
# Magnetization field in magnet
V_Nmagnet = TestFESpace(Ωmagnet, reffeu)
Nhmagnet  = interpolate_everywhere((x)->VectorValue(0.0, 1.0), V_Nmagnet)

# Remanent magnetisation of the magnet, ramped by the staggered loop: called before every step with the
# stage τ⁻ ∈ [0, 1] of the last accepted step and the increment ∆τ; equal steps, step n = τ⁻/∆τ
function αr_update(τ⁻, ∆τ)
  model_magnet.αr = αr_max * ∆τ * (round(τ⁻ / ∆τ) + 1)   # HyperFEM dereferences αr itself (PR #158)
end
 
# Problem 1: residual and jacobian
res_mag(Λ) = (φ, vφ) -> -1.0 * ∫((∇(vφ) ⋅ (∂Ψs∂H0 ∘ (F ∘ (∇(uhsolid(Λ))'), ℋ₀ ∘ (∇(φ)), Nhbeam))))dΩsolid -
                        ∫((∇(vφ) ⋅ (∂Ψv∂H0 ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φ))))))dΩvacuum_mag -
                        ∫((∇(vφ) ⋅ (model_magnet(Λ)[2] ∘ (ℋ₀ ∘ (∇(φ)), Nhmagnet))))dΩmagnet

jac_mag(Λ) = (φ, dφ, vφ) -> ∫(∇(vφ)' ⋅ ((∂Ψs∂H0H0 ∘ (F ∘ (∇(uhsolid(Λ))'), ℋ₀ ∘ (∇(φ)), Nhbeam)) ⋅ ∇(dφ)))dΩsolid +
                            ∫(∇(vφ)' ⋅ ((∂Ψv∂H0H0 ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φ)))) ⋅ ∇(dφ)))dΩvacuum_mag +
                            ∫(∇(vφ)' ⋅ ((model_magnet(Λ)[3] ∘ (ℋ₀ ∘ (∇(φ)), Nhmagnet)) ⋅ ∇(dφ)))dΩmagnet

# Problem 2: residual and jacobian
res_mech(Λ) = begin
    # Assemble Maxwell stress in vacuum
    assemble_vector!((v) -> ∫(∇(v)' ⊙ (∂Ψv∂F ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φh(Λ))))))dΩvacuum_mec, Pmaxw_, V_maxwell)
    # nodal values on the interface of the solid
    transfer!(tmaxw_, T_maxwell, Pmaxw)
    Massoperator(tmaxw_)
     (u, v) -> ∫((∇(v)' ⊙ (∂Ψs∂F ∘ (F ∘ (∇(u)'), ℋ₀ ∘ (∇(φh(Λ))), Nhbeam))))dΩsolid + ∫((v ⋅ tmaxw))dΓs
end

jac_mech(Λ) = (u, du, v) -> ∫(∇(v)' ⊙ ((∂Ψs∂FF ∘ (F ∘ (∇(u)'), ℋ₀ ∘ (∇(φh(Λ))), Nhbeam)) ⊙ (∇(du)')))dΩsolid

# Problem 3: residual and jacobian
res_vacmech(Λ) = (u, v) -> ∫((∇(v)' ⊙ (∂Ψvm∂F ∘ (F ∘ (∇(u)')))))dΩvacuum_mec
jac_vacmech(Λ) = (u, du, v) -> ∫(∇(v)' ⊙ ((∂Ψvm∂FF ∘ (F ∘ (∇(u)'))) ⊙ (∇(du)')))dΩvacuum_mec


# ******************************************************
#             Computational problems
#******************************************************

# Problem 1
nls_mag = NewtonRaphson(LUSolver(); maxiter=10, atol=1.e-9, rtol=1.e-8, verbose=true)
comp_model_mag = StaticNonlinearModel(res_mag, jac_mag, Uφ⁺, Vφ, Dφ; nls=nls_mag, xh=φh⁺)

# Problem 2
nls_mech = NewtonRaphson(LUSolver(); maxiter=20, atol=1.e-6, rtol=1.e-2, verbose=true)
comp_model_mech = StaticNonlinearModel(res_mech, jac_mech, Uu_solid⁺, Vu_solid, Du_solid; nls=nls_mech, xh=uh_solid⁺)

# Problem 3 
α = CellState(1.0, dΩvacuum_mec)
linesearch = Injectivity_Preserving_LS(α, Uu_vacuum⁺, Vu_vacuum; maxiter=50, αmin=1e-16, ρ=0.5, c=0.95)
nls_vacmech = Newton_RaphsonSolver(LUSolver(); maxiter=10, rtol=2, verbose=true, linesearch=linesearch)
comp_model_vacmech = StaticNonlinearModel(res_vacmech, jac_vacmech, Uu_vacuum⁺, Vu_vacuum, Du_vacuum; nls=nls_vacmech, xh=uh_vacuum⁺)

# Staggered resolution:  Problem1-> Problem2-> Problem3
comp_model = StaggeredModel((comp_model_mag, comp_model_mech, comp_model_vacmech), (φh⁺, uh_solid⁺, uh_vacuum⁺), (φh⁻, uh_solid⁻, uh_vacuum⁻))

# ******************************************************
#               Run solver
#******************************************************

args_mag = Dict(:stepping => (nsteps=1, maxbisec=5), :ProjectDirichlet=>true)
args_mech = Dict(:stepping => (nsteps=1, maxbisec=5))
args_vacmech = Dict(:stepping => (nsteps=1, maxbisec=5), :ProjectDirichlet=>true)
args = (args_mag, args_mech, args_vacmech)

solve!(comp_model; stepping=StaggeredStepping(nsteps=20, max_substeps=2, maxbisec=1), presolver=αr_update, kargsolve=args)

writevtk(Ωsolid,simdir*"/hMRE2",cellfields=["uh"=>uh_solid⁺,"-∇(φh)" => -∇(φh⁺)]) # MRE deformation
writevtk(Ωvacuum,simdir*"/vacuum",cellfields=["uh"=>uh_vacuum⁺, "-∇(φh)" => -∇(φh⁺)]) # vacuum
writevtk(Ωmagnet,simdir*"/magnet",cellfields=["-∇(φh)" => -∇(φh⁺)]) # Magnet