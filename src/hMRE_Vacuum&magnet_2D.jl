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
using GridapSolvers.NonlinearSolvers
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
params      =  [0.010000000052968838, 12480.64495232286, 1.0, 2.0, 5195.545287237134, 0.2602717127043121,  2.0]
modelmech   =  NonlinearMooneyRivlin2D(λ=(params[1] + params[2]) * 1e2, μ1=params[1]*1e3, μ2=params[2], α1=params[3], α2=params[4]) + 
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
Γair_int   = BoundaryTriangulation(Ωvacuum, tags="Interface")
Γsolid_int = BoundaryTriangulation(Ωsolid, tags="Interface")
Γsf        = InterfaceTriangulation(Ωsolid, Ωvacuum)
nΓsf       = get_normal_vector(Γsf)
dΓsf       = Measure(Γsf, 4 * order_solid)


#************************************
# Dirichlet boundary conditions
#************************************

# Problem 1: Dirichlet conditions magnetic field
evolφ(Λ) = Λ
dir_φ_tags = ["phi_Dirichlet"]
dir_φ_values = [0.0]
dir_φ_timesteps = [evolφ]
Dφ = DirichletBC(dir_φ_tags, dir_φ_values, dir_φ_timesteps)

# Problem 2: Dirichlet conditions solid
evolu(Λ) = 1.0
dir_u_tags_solid = ["Fixed_Solid", "Fixed_Solid_x"]
dir_u_values_solid = [[0.0, 0.0], [0.0, 0.0]]
dir_u_timesteps_solid = [evolu, evolu]
Du_solid = DirichletBC(dir_u_tags_solid, dir_u_values_solid, dir_u_timesteps_solid)

# Problem 3: Dirichlet conditions vacuum movement
dir_u_tags_air      = ["Fixed_Air", "Fixed_Air_x", "Interface"]
dir_u_values_air    = [[0.0, 0.0], [0.0, 0.0], DirichletCoupling]
dir_u_timesteps_air = [evolu, evolu, evolu]
Du_vacuum = DirichletBC(dir_u_tags_air, dir_u_values_air, dir_u_timesteps_air)

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
Vu_solid = TestFESpace(Ωsolid, reffeu, Du_solid, conformity=:H1, dirichlet_masks=[(true, true), (true, false)])
Uu_solid⁺ = TrialFESpace(Vu_solid, Du_solid)
Uu_solid⁻ = TrialFESpace(Vu_solid, Du_solid)
uh_solid⁺ = FEFunction(Uu_solid⁺, zero_free_values(Uu_solid⁺))
uh_solid⁻ = FEFunction(Uu_solid⁻, zero_free_values(Uu_solid⁻))

# Problem 3: FE Spaces and state variables 
Vu_vacuum = TestFESpace(Ωvacuum, reffeu_vacuum, Du_vacuum, conformity=:H1, dirichlet_masks=[(true, true), (true, false), (true, true)])
Uu_vacuum⁺ = TrialFESpace(Vu_vacuum, Du_vacuum)
Uu_vacuum⁻ = TrialFESpace(Vu_vacuum, Du_vacuum)
uh_vacuum⁺ = FEFunction(Uu_vacuum⁺, zero_free_values(Uu_vacuum⁺))
uh_vacuum⁻ = FEFunction(Uu_vacuum⁻, zero_free_values(Uu_vacuum⁻))

# ******************************************************
# Solid–vacuum coupling
# Function to interpolate staggered solid displacement
# increments to the vacuum region
#******************************************************

# Interface solid-vacuum coupling
Uair_int = FESpace(Γair_int, reffeu_vacuum)
Usolid_int = FESpace(Γsolid_int, reffeu)

uhsolid(Λ) = uh_solid⁻ + (uh_solid⁺ - uh_solid⁻) * Λ

uhsolid_int = interpolate_everywhere(uh_solid⁺, Usolid_int)
uhair_int = interpolate_everywhere(Interpolable(uhsolid_int), Uair_int)
uhsolid_int_(Λ) = Interpolable(interpolate_everywhere!(uhsolid(Λ), get_free_dof_values(uhsolid_int), uhsolid_int.dirichlet_values, Usolid_int))
uhair_int_(Λ) = interpolate_everywhere!(uhsolid_int_(Λ), get_free_dof_values(uhair_int), uhair_int.dirichlet_values, Uair_int)

InterpolableBC!(Uu_vacuum⁺, Du_vacuum, "Interface", uhair_int_)

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

function αr_update(τ,∆τ)
model_magnet.αr[] = αr_max * ∆τ * (τ + 1)
end
 
# Problem 1: residual and jacobian
res_mag(Λ) = (φ, vφ) -> -1.0 * ∫((∇(vφ) ⋅ (∂Ψs∂H0 ∘ (F ∘ (∇(uhsolid(Λ))'), ℋ₀ ∘ (∇(φ)), Nhbeam))))dΩsolid -
                        ∫((∇(vφ) ⋅ (∂Ψv∂H0 ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φ))))))dΩvacuum_mag -
                        ∫((∇(vφ) ⋅ (model_magnet(Λ)[2] ∘ (ℋ₀ ∘ (∇(φ)), Nhmagnet))))dΩmagnet

jac_mag(Λ) = (φ, dφ, vφ) -> ∫(∇(vφ)' ⋅ ((∂Ψs∂H0H0 ∘ (F ∘ (∇(uhsolid(Λ))'), ℋ₀ ∘ (∇(φ)), Nhbeam)) ⋅ ∇(dφ)))dΩsolid +
                            ∫(∇(vφ)' ⋅ ((∂Ψv∂H0H0 ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φ)))) ⋅ ∇(dφ)))dΩvacuum_mag +
                            ∫(∇(vφ)' ⋅ ((model_magnet(Λ)[3] ∘ (ℋ₀ ∘ (∇(φ)), Nhmagnet)) ⋅ ∇(dφ)))dΩmagnet

# Problem 2: residual and jacobian
res_mech(Λ) = (u, v) -> ∫((∇(v)' ⊙ (∂Ψs∂F ∘ (F ∘ (∇(u)'), ℋ₀ ∘ (∇(φh(Λ))), Nhbeam))))dΩsolid -
                        ∫((v.⁺ ⋅ ((∂Ψv∂F ∘ (F ∘ (∇(uhvacuum(Λ))'), ℋ₀ ∘ (∇(φh(Λ))))).⁻ ⋅ nΓsf.⁺)))dΓsf

jac_mech(Λ) = (u, du, v) -> ∫(∇(v)' ⊙ ((∂Ψs∂FF ∘ (F ∘ (∇(u)'), ℋ₀ ∘ (∇(φh(Λ))), Nhbeam)) ⊙ (∇(du)')))dΩsolid

# Problem 3: residual and jacobian
res_vacmech(Λ) = (u, v) -> ∫((∇(v)' ⊙ (∂Ψvm∂F ∘ (F ∘ (∇(u)')))))dΩvacuum_mec
jac_vacmech(Λ) = (u, du, v) -> ∫(∇(v)' ⊙ ((∂Ψvm∂FF ∘ (F ∘ (∇(u)'))) ⊙ (∇(du)')))dΩvacuum_mec


# ******************************************************
#             Computational problems
#******************************************************

# Problem 1
nls_mag = NewtonSolver(LUSolver(); maxiter=10, atol=1.e-9, rtol=1.e-8, verbose=true)
comp_model_mag = StaticNonlinearModel(res_mag, jac_mag, Uφ⁺, Vφ, Dφ; nls=nls_mag, xh=φh⁺)

# Problem 2
nls_mech = NewtonSolver(LUSolver(); maxiter=20, atol=1.e-6, rtol=1.e-2, verbose=true)
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

solve!(comp_model; stepping=(nsteps=20, nsubsteps=2, maxbisec=1), presolver=αr_update, kargsolve=args)

writevtk(Ωsolid,simdir*"/hMRE2",cellfields=["uh"=>uh_solid⁺]) # MRE deformation
writevtk(Ωvacuum,simdir*"/vacuum",cellfields=["uh"=>uh_vacuum⁺, "-∇(φh)" => -∇(φh⁺)]) # vacuum
writevtk(Ωmagnet,simdir*"/magnet",cellfields=["-∇(φh)" => -∇(φh⁺)]) # Magnet