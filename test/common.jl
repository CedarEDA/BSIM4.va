using DataFrames
using CedarSim
using CedarSim: SweepFlattener, ParamLens

const repo_root = Base.pkgdir(BSIM4)

"""
    load_test_model(modelcard_path, model_path)

Load a BSIM4 test model, provided via the modelcard path (e.g. `"bsim_group/nmos.modelcard"`)
and a model path (e.g. `"bsim_group/nmos.spice"`).  Returns the AST and the circuit function,
so that the user can alter the spice AST, and compile it for use with CedarSim.
"""
function load_test_model(model_path)
    m = Module()
    Core.eval(m, quote
        using CedarSim, BSIM4
        using CedarSim.SpectreEnvironment
    end)

    # Load model card into that module
    #Base.include(m, SpcFile(modelcard_path, true))

    # Load spice circuit into that module
    circuit_ast = Core.eval(m, :(CedarSim.SpectreNetlistParser.SPICENetlistParser.SPICENetlistCSTParser.parsefile($(model_path))))
    circuit_code = CedarSim.make_spectre_circuit(circuit_ast, [dirname(model_path)])
    circuit = Core.eval(m, circuit_code)
    
    # Call the circuit once, to lint it:
    invokelatest(circuit)

    # Return this constructed circuit function, and the parsed source
    return circuit_ast, circuit
end

"""
    ngspice_sweep(ps)

Convert a `ProductSweep` into an NgSpice sweep representation, e.g.:

    ngspice_sweep(ProductSweep(v_drain=0.0:0.02:1.2, v_gate=0.2:0.2:1.2))

gives:

    "v_drain 0.0 1.2 0.02 v_gate 0.2 1.2 0.2"

This is used by `ngspice_dc_command()`, which is used to emit the correct DC inner
sweep command for our testset.
"""
function ngspice_sweep(ps, aliases)
    # First, make sure that they're all ranges
    for sweep in ps.iterator.iterators
        if !isa(sweep.values, AbstractRange)
            throw(ArgumentError("Cannot emit an NgSpice inner sweep with non-range axis '$(sweep.selector)'"))
        end
    end

    # Construct a string representation of these sweeps:
    ngspice_repr(sweep) = "$(get(aliases, string(sweep.selector), string(sweep.selector))) $(first(sweep.values)) $(last(sweep.values)) $(step(sweep.values))"
    return join(ngspice_repr.(ps.iterator.iterators), " ")
end

function ngspice_dc_command(io::IO, inner_sweep, aliases, probes)
    println(io, """
    .dc $(ngspice_sweep(inner_sweep, aliases))
    .print dc $(join(probes, " "))
    .end
    """)
end

function ngspice_dc_solve(ast, outer_sweep::SweepFlattener, inner_sweep::SweepFlattener, probes::Dict{String,String}, aliases::Dict{String,String}, files_to_copy::Vector{String})
    dfs = []
    for params in outer_sweep
        mktempdir() do dir
            spice_path = joinpath(dir, "test.spice")
            open(spice_path, write=true) do io
                # For some reason, `alter` seems to be dropping the initial comment
                println(io, "* Generated spice")
    
                # Write out the altered spice code
                nt_params = NamedTuple(params)
                alter(io, ast, ParamLens(nt_params))

                # Special-case temperature
                if hasfield(typeof(nt_params), :temp)
                    println(io, ".options temp=$(nt_params.temp)")
                end
                
                # Add `.dc` commands:
                ngspice_dc_command(io, inner_sweep, aliases, collect(values(probes)))
            end

            # Copy in extra files, like model cards, since our generated spice
            # still contains `.include` and `.lib` commands.  :(
            for file in files_to_copy
                cp(file, joinpath(dir, basename(file)))
            end

            # Run the simulation in a clean environment
            NgSpice.cmd("destroy all")
            NgSpice.cmd("source $(spice_path)")
            NgSpice.cmd("run")
    
            # Pull out every probe we asked for
            probe_map = Dict(Symbol(k) => NgSpice.getrealvec(v) for (k, v) in probes)
            
            # Also add in all our params:
            for (name, val) in params
                probe_map[name] = repeat([val], length(first(values(probe_map))))
            end
            
            nested_inner_params = CedarSim.nest_param_list.(collect(inner_sweep)[:])
            for name in sweepvars(inner_sweep)
                probe_map[name] = Float64[getfield(p, name) for p in nested_inner_params]
            end
            push!(dfs, DataFrame(probe_map...))
        end
    end
    return vcat(dfs...)
end

"""
    cedarsim_dc_solve(circuit, sweep, probes::Dict{String,String})

Take `circuit`, solve it for all 
"""
function cedarsim_dc_solve(circuit, sweep::SweepFlattener, probes::Dict{String,String}; abstol=1e-9)
    lens(probe::String) = Accessors.opticcompose(PropertyLens.(Symbol.(split(probe, ".")))...)
    cs = CircuitSweep(circuit, sweep; debug_config=(;verify_ir_levels = true))
    sols = dc!(cs; abstol)[:]

    rows = Dict[]
    for (sol, params) in zip(sols, sweep)
        # Extract the probes from the solution object
        probe_vals = Dict{Symbol,Float64}(Symbol(name) => sol[lens(probe)(cs.sys)] for (name, probe) in probes)
        # Add in our simulation values
        for (name, val) in params
            probe_vals[name] = val
        end
        push!(rows, probe_vals)
    end
    return DataFrame(rows)
end


"""
    CtoK(C)

Convert a temperature in Celsius to Kelvin
"""
CtoK(C) = C + 273.15

"""
    CtoK(C)

Convert a temperature in Kelvin to Celsius
"""
KtoC(K) = K - 273.15
