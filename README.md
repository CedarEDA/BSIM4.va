# Verilog-A translation of BSIM4

This is a Verilog-A implementation of the BSIM4 model,
hand translated from the original C code. This translation in particular features:

- Integrated support for multiple versions of BSIM4 (via the `VERSION` model parameter)
- Familiar structure similar to more recent CMC models
- Permissive license and known origin

The model in this repository is the default model used
by CedarSim when the netlist calls for BSIM4.

> [!WARNING]
> Translation and validation is not fully complete. Some lesser-used parameter options have not yet been translated or are known to have minor bugs. While this model is used by default in Cedar and is known to match (to high numerical precision) results from other simulators on tested PDKs, not all parameter options have been validated.

## License / Contributing

The Cedar EDA platform is dual-licensed under a commercial license and CERN-OHL-S v2. In addition, some packages (including this one) are also available under the MIT license. Please see the LICENSE file for more
information and the LICENSE.FAQ.md file for more information on how to
use Cedar under the CERN-OHL-S v2 license.

We are accepting PRs on all Cedar repositories, although you must sign the Cedar Contributor License Agreement (CLA) for us to be able to incorporate your changes into the upstream repository. Additionally, if you would like to make major or architectural changes, please discuss this with us *before* doing the work. Cedar is a complicated piece of software, with many moving pieces, not all of which are publicly available. Because of this, we may not be able to take your changes, even if they are correct and useful (so again, please talk to us first).
