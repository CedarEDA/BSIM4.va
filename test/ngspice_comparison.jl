using DAECompiler, CedarSim, Test, OrdinaryDiffEq, DataFrames, BSIM4, CedarSim.Accessors
using CedarSim.SpectreEnvironment
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
nmos_sweep = ProductSweep(v_drain=0.0: 0.02: 1.2, v_gate=0.0: 0.2: 1.2, v_bulk=0.0:-0.3:-1.2, temp=[-55., 100.], NF=[1.,5.], gmin=[0.])
pmos_sweep = ProductSweep(v_drain=0.0:-0.02:-1.2, v_gate=0.0:-0.2:-1.2, v_bulk=0.0: 0.3: 1.2, temp=[-55., 100.], NF=[1.,5.], gmin=[0.])

function compare_sweep(ast, circuit, modelcard, sweep; abstol=1e-9)
    # Build a mapping from dataframe column name to the NgSpice probe name:
    ngspice_probes = Dict(
        "gate voltage" => "v(gate)",
        "drain voltage" => "v(drain)",
        "drain current" => "i(vd)",
    )

    # Cedarsim "probes" are symbols that index into `sys`
    cedarsim_probes = Dict(
        "gate voltage" => "vg.V",
        "drain voltage" => "vd.V",
        "drain current" => "vd.I",
    )

    # Split v_drain and v_gate to be an inner sweep for NgSpice.
    outer, inner = split_axes(sweep, [:v_drain, :v_gate])

    # The NgSpice `.dc` operation operates not on `.param` objects, but on actual voltage nodes,
    # so re-name our inner sweep's axes to refer to `vd` and `vg`:
    ngspice_dc_aliases = Dict("v_drain" => "vd", "v_gate" => "vg")

    # Run the DC solves
    df_ngspice = ngspice_dc_solve(ast, outer, inner, ngspice_probes, ngspice_dc_aliases, [modelcard])
    df_cedarsim = cedarsim_dc_solve(circuit, sweep, cedarsim_probes; abstol)

    # Ensure that the ordering of all parameters is identical:
    @testset "parameter order" begin
        @test df_ngspice[:, collect(sweepvars(sweep))] == df_cedarsim[:, collect(sweepvars(sweep))]
    end

    # Ensure all probes are within tolerance
    for probe in Symbol.(keys(cedarsim_probes))
        @testset "$(probe)" begin
            @test all(isapprox.(df_cedarsim[:, probe], df_ngspice[:, probe]; atol=abstol*10))
        end
    end
    return df_cedarsim, df_ngspice
end

@testset "ngspice comparison" begin
    @testset "nmos" begin 
        compare_sweep(nmos_ast, nmos_circuit, nmos_modelcard, nmos_sweep)
    end

    @testset "pmos" begin 
        compare_sweep(pmos_ast, pmos_circuit, pmos_modelcard, pmos_sweep)
    end
end
