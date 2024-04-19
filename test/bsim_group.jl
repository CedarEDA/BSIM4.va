using DAECompiler, CedarSim, Test, OrdinaryDiffEq, DataFrames, BSIM4, CedarSim.Accessors
using CedarSim.SpectreEnvironment
using Random
include(joinpath(Base.pkgdir(BSIM4), "test", "common.jl"))

# Load our models, one for the nmos, one for the pmos variant
nmos_ast, nmos_circuit = load_test_model(
    joinpath(repo_root, "test", "bsim_group", "nmos.spice"),
)
pmos_ast, pmos_circuit = load_test_model(
    joinpath(repo_root, "test", "bsim_group", "pmos.spice"),
)

# The BSIM4 test suite is broken into two parts; nmos and pmos
# Each part contains a sweep over multiple parameters.
# The upstream tests separate these into 14 disparate simulation outputs.
# We mimic that here, to make it easier to compare against other simulators.
test_suites = Dict(
    :nmos => [
        "test1"  => (nmos_circuit, ProductSweep(v_drain=0.0:0.02:1.18, v_gate=0.2:0.2:1.2, NF=[5.], gmin=[0.])),
        "test2"  => (nmos_circuit, ProductSweep(v_drain=0.0:0.02:1.18, v_gate=0.0:0.2:1.2, temp=[-55.], gmin=[0.])),
        "test3"  => (nmos_circuit, ProductSweep(v_drain=0.0:0.02:1.18, v_gate=0.0:0.2:1.2, temp=[100.], gmin=[0.])),
        "test4"  => (nmos_circuit, ProductSweep(v_gate= 0.0:0.02:1.18, v_drain=0.05:0.5:1.2, gmin=[0.])),
        "test5"  => (nmos_circuit, ProductSweep(v_gate=-0.6:0.02:1.18, v_bulk=0.0:-0.3:-1.2, v_drain=[0.1], gmin=[0.])),
        "test6"  => (nmos_circuit, ProductSweep(v_gate= 0.6:0.02:1.18, v_bulk=0.0:-0.3:-1.2, v_drain=[0.1], temp=[-55.], gmin=[0.])),
        "test7"  => (nmos_circuit, ProductSweep(v_gate=-0.6:0.02:1.18, v_bulk=0.0:-0.3:-1.2, v_drain=[0.1], temp=[100.], gmin=[0.])),
    ],
    :pmos => [
        "test8"  => (pmos_circuit, ProductSweep(v_drain=0.0:-0.02:-1.18, v_gate=-0.2:-0.2:-1.2, NF=[5.], gmin=[0.])),
        "test9"  => (pmos_circuit, ProductSweep(v_drain=0.0:-0.02:-1.18, v_gate=0.0:-0.3:-1.2, temp=[-55.], gmin=[0.])),
        "test10" => (pmos_circuit, ProductSweep(v_drain=0.0:-0.02:-1.18, v_gate=0.0:-0.2:-1.2, temp=[100.], gmin=[0.])),
        "test11" => (pmos_circuit, ProductSweep(v_gate= 0.6:-0.02:-1.18, v_drain=-0.1:-0.3:-1.2, gmin=[0.])),
        "test12" => (pmos_circuit, ProductSweep(v_gate= 0.6:-0.02:-1.18, v_bulk=0.0:0.3:1.2, v_drain=[-0.1], gmin=[0.])),
        "test13" => (pmos_circuit, ProductSweep(v_gate= 0.6:-0.02:-1.18, v_bulk=0.0:0.3:1.2, v_drain=[-0.1], temp=[-55.], gmin=[0.])),
        "test14" => (pmos_circuit, ProductSweep(v_gate= 0.6:-0.02:-1.18, v_bulk=0.0:0.3:1.2, v_drain=[-0.1], temp=[100.], gmin=[0.])),
    ]
)

# Cedarsim "probes" are symbols that index into `sys`
cedarsim_probes = Dict(
    "gate voltage" => "vg.V",
    "drain voltage" => "vd.V",
    "drain current" => "vd.I",
)
const abstol = 1e-9


function parse_out_file(out_file::String)
    vals = Float64[]
    for line in readlines(out_file)
        m = match(r"^(\d+)\s+[\d\.e+\-]+\s+([\d\.e+\-]+)\s*$", line)
        if m !== nothing
            if parse(Int64, m.captures[1]) != length(vals)
                error("missed a line!  We've parsed $(length(vals)) but we should have parsed $(m.captures[1])")
            end
            push!(vals, parse(Float64, m.captures[2]))
        end
    end
    return vals
end

@testset "bsim_group" begin
    for (circuit_type, test_sweeps) in test_suites, (test_name, (circuit, sweep)) in test_sweeps
        @testset "$(circuit_type) - $(test_name)" begin
            df_cedarsim = cedarsim_dc_solve(circuit, sweep, cedarsim_probes; abstol)
            drain_current = parse_out_file(joinpath(@__DIR__, "bsim_group", "$(test_name).out"))

            max_divergence = maximum(abs.(df_cedarsim[:, Symbol("drain current")] .- drain_current))
            if max_divergence > abstol*10
                @error("$(circuit_type) - $(test_name)", max_divergence, atol=abstol*10)
            end
            @test max_divergence < abstol*10 broken=(test_name ∈ ("test1", "test8"))
        end
    end
end
