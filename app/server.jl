# =============================================================================
# Web interface for CuttingPlane.jl (development entry point)
#
#   julia --project=app app/server.jl          # then open http://127.0.0.1:8000
#   PORT=9000 julia --project=app app/server.jl
#
# The server itself lives in the CutsApp package (app/src/); the desktop app
# runs the same code, compiled with PackageCompiler (see build/).
# =============================================================================

using CutsApp

CutsApp.serve(; port = parse(Int, get(ENV, "PORT", "8000")))
