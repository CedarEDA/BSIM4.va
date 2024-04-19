using DAECompiler, CedarSim, Test, OrdinaryDiffEq, DataFrames, BSIM4, CedarSim.Accessors
using OrdinaryDiffEq.LineSearches
using CedarSim.SpectreEnvironment
using FiniteDiff
using NgSpice
using Random
include(joinpath(Base.pkgdir(BSIM4), "test", "common.jl"))

# Need to copy over the modelcards as long as `.include` is not inlined in our spice parser,
# so we pass these `modelcard` variables in to `ngspice_dc_solve()`.  Eventually,
# we'll inline all `.include` and `.lib` statements when doing a spice alter, which
# will obviate the need for this.
nmos_modelcard = joinpath(repo_root, "test", "bsim_group", "nmos.modelcard")
pmos_modelcard = joinpath(repo_root, "test", "bsim_group", "pmos.modelcard")

nmos_ast, nmos_circuit = load_test_model(
    joinpath(repo_root, "test", "bsim_group", "nmos.spice"),
)
pmos_ast, pmos_circuit = load_test_model(
    joinpath(repo_root, "test", "bsim_group", "pmos.spice"),
)

# Do a large sweep over all of the products
nmos_sweep = ProductSweep(v_drain=0.0: 0.02: 1.2, v_gate=0.0: 0.2: 1.2, v_bulk=0.0:-0.3:-1.2, temp=[CtoK(-55.), CtoK(100.)], NF=[1.,5.], gmin=[0.])
pmos_sweep = ProductSweep(v_drain=0.0:-0.02:-1.2, v_gate=0.0:-0.2:-1.2, v_bulk=0.0: 0.3: 1.2, temp=[CtoK(-55.), CtoK(100.)], NF=[1.,5.], gmin=[0.])

function make_dc_prob(::Type{DAEProblem}, cs, circuit)
    return DAEProblem(cs.sys, nothing, nothing, (0.0, 1.0), circuit; initializealg=CedarDCOp(), jac=true)
end
function make_dc_prob(::Type{ODEProblem}, cs, circuit)
    return ODEProblem(cs.sys, nothing, (0.0, 1.0), circuit; initializealg=CedarDCOp(), jac=true)
end
solver(::DAEProblem) = DFBDF(;nlsolve=NLNewton(relax=BackTracking()))
solver(::ODEProblem) = Rosenbrock23()

function jac_diffractor(prob::DAEProblem, sol)
    J = zeros(length(prob.u0), length(prob.u0))
    u = sol.u
    du = sol.du
    γ = 1.0
    t = 0.0
    prob.f.jac(J, du, u, prob.p, γ, t)
    return J
end
function jac_finitediff(prob::DAEProblem, sol)
    J = zeros(length(prob.u0), length(prob.u0))
    J2 = copy(J)
    u = sol.u
    du = sol.du
    γ = 1.0
    t = 0.0
    #=dG/du =# FiniteDiff.finite_difference_jacobian!(J, (_out, _u)->prob.f.f(_out, du, _u, prob.p, t), u)
    #=dG/ddu=# FiniteDiff.finite_difference_jacobian!(J2, (_out, _du)->prob.f.f(_out, _du, u, prob.p, t), du)
    J .+= γ*J2
    return J
end
function jac_diffractor(prob::ODEProblem, sol)
    J = zeros(length(prob.u0), length(prob.u0))
    u = sol.u
    t = 0.0
    prob.f.jac(J, u, prob.p, t)
    return J
end
function jac_finitediff(prob::ODEProblem, sol)
    J = zeros(length(prob.u0), length(prob.u0))
    u = sol.u
    γ = 1.0
    t = 0.0
    #=dG/du =# FiniteDiff.finite_difference_jacobian!(J, (_out, _u)->prob.f.f(_out, _u, prob.p, t), u)
    return J
end

function compare_jacobians(probtype, circuit, sweep; abstol=1e-9)
    # Run the DC solves
    cs = CircuitSweep(circuit, sweep)
    prob = make_dc_prob(probtype, cs, circuit)

    # Get an initial `u0`/`du0`, then use that to compare jacobians
    for sim in cs
        prob = remake(prob, p=sim)
        sol = OrdinaryDiffEq.init(prob, solver(prob))

        # Collect jacobians from diffractor and finitediff, ensure they are equal:
        J_diffractor = jac_diffractor(prob, sol)
        J_finitediff = jac_finitediff(prob, sol)

        @testset "$(sim.params)" begin
            max_divergence = maximum(abs.(J_diffractor .- J_finitediff))
            # This threshold is so high because FiniteDiff can run into a lot of trouble.
            @test max_divergence <= abstol*10000
        end
    end
end

@testset "Jacobian comparison" begin
    @testset "nmos" begin
        @testset "DAEProblem" begin
            compare_jacobians(DAEProblem, nmos_circuit, nmos_sweep)
        end
        @testset "ODEProblem" begin
            compare_jacobians(ODEProblem, nmos_circuit, nmos_sweep)
        end
    end

    @testset "pmos" begin 
        @testset "DAEProblem" begin
            compare_jacobians(DAEProblem, pmos_circuit, pmos_sweep)
        end
        @testset "ODEProblem" begin
            compare_jacobians(ODEProblem, pmos_circuit, pmos_sweep)
        end
    end
end
