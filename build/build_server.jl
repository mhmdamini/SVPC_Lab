# =============================================================================
# Compile the web-interface server (app/, package CutsApp) into a standalone
# program with PackageCompiler. The result bundles Julia, every package and
# the HiGHS solver, so users need nothing installed.
#
#   julia --project=build -e 'using Pkg; Pkg.instantiate()'
#   julia --project=build build/build_server.jl
#
# Output: build/server/bin/cuts-server(.exe). The Electron shell in desktop/
# ships this folder and starts the program with PORT=<free port>.
# Takes 20–40 minutes and a few GB of temporary disk space.
# =============================================================================

using PackageCompiler

const ROOT = dirname(@__DIR__)
const APP = joinpath(ROOT, "app")
const OUT = joinpath(ROOT, "build", "server")

# Make sure the app environment (and its path dependency on CuttingPlane) is installed.
import Pkg
Pkg.activate(APP)
Pkg.instantiate()
Pkg.activate(@__DIR__)

create_app(APP, OUT;
           executables = ["cuts-server" => "julia_main"],
           precompile_execution_file = joinpath(APP, "precompile_app.jl"),
           include_lazy_artifacts = true,
           force = true)

println("\nServer compiled to ", OUT)
