"""
    CutsApp

Web interface for `CuttingPlane`: an HTTP server that serves the single-page UI
(`app/index.html`) and a JSON API (`POST /api/solve`).

* Development: `julia --project=app app/server.jl`, then open http://127.0.0.1:8000.
* Desktop app: `julia_main` is the entry point compiled by PackageCompiler
  (`build/build_server.jl`); the Electron shell in `desktop/` starts it.
"""
module CutsApp

using CuttingPlane
using HTTP
using JSON3
using JuMP
using LinearAlgebra

include("geometry.jl")   # vertex enumeration for the 2-D / 3-D plots
include("api.jl")        # parse → generate cuts → JSON
include("iterate.jl")    # cutting-plane loop, branch and bound / cut → JSON

# ---- Static files -------------------------------------------------------------

# The page is embedded at compile time so a compiled app needs no source tree;
# while developing, the file on disk is served instead so edits show up on reload.
const INDEX_HTML_PATH = normpath(joinpath(@__DIR__, "..", "index.html"))
include_dependency(INDEX_HTML_PATH)
const INDEX_HTML_EMBEDDED = read(INDEX_HTML_PATH)

index_html() = isfile(INDEX_HTML_PATH) ? read(INDEX_HTML_PATH) : INDEX_HTML_EMBEDDED

# Plotly is served from a local file when the desktop app provides one
# (CUTS_PLOTLY_JS), so the app works offline; otherwise from Plotly's CDN.
const PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"

function plotly_js(_)
    path = get(ENV, "CUTS_PLOTLY_JS", "")
    isfile(path) || return HTTP.Response(302, ["Location" => PLOTLY_CDN])
    return HTTP.Response(200, ["Content-Type" => "application/javascript",
                               "Cache-Control" => "max-age=86400"], read(path))
end

# ---- HTTP plumbing ------------------------------------------------------------

json_response(status, body) =
    HTTP.Response(status, ["Content-Type" => "application/json"], JSON3.write(body))

# Request body as bytes (HTTP.jl 1.x stores a Vector{UInt8}, 2.x wraps it in a body object).
body_bytes(body::AbstractVector{UInt8}) = body
body_bytes(body) = hasproperty(body, :data) ? body.data : read(body)

function handle_solve(http_req::HTTP.Request)
    try
        req = JSON3.read(body_bytes(http_req.body))
        return json_response(200, solve_request(req))
    catch err
        @warn "Request failed" exception = (err, catch_backtrace())
        return json_response(400, Dict("error" => sprint(showerror, err)))
    end
end

function handle_iterate(http_req::HTTP.Request)
    try
        req = JSON3.read(body_bytes(http_req.body))
        return json_response(200, iterate_request(req))
    catch err
        @warn "Request failed" exception = (err, catch_backtrace())
        return json_response(400, Dict("error" => sprint(showerror, err)))
    end
end

function router()
    r = HTTP.Router()
    HTTP.register!(r, "GET", "/", _ -> HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], index_html()))
    HTTP.register!(r, "GET", "/vendor/plotly.min.js", plotly_js)
    HTTP.register!(r, "GET", "/health", _ -> HTTP.Response(200, "ok"))
    HTTP.register!(r, "POST", "/api/solve", handle_solve)
    HTTP.register!(r, "POST", "/api/iterate", handle_iterate)
    return r
end

"""
    start(; host = "127.0.0.1", port = 8000) -> server

Start the server in the background and return it (`close(server)` stops it).
"""
start(; host::AbstractString = "127.0.0.1", port::Integer = 8000) = HTTP.serve!(router(), host, port)

"""
    serve(; host = "127.0.0.1", port = 8000)

Start the server and block until it is stopped (Ctrl+C).
"""
function serve(; host::AbstractString = "127.0.0.1", port::Integer = 8000)
    server = start(; host, port)
    println("Disjunctive cuts interface running at http://$host:$port  (Ctrl+C to stop)")
    flush(stdout)
    wait(server)
end

"""
    exit_when_stdin_closes()

Stop the process once standard input reaches end-of-file. The desktop app
keeps a pipe to our stdin open; the OS closes it when the app exits for any
reason (including a crash), so the server never outlives its window.
"""
exit_when_stdin_closes() = @async begin
    try
        read(stdin)
    catch
    end
    exit(0)
end

"""
    julia_main() -> Cint

Entry point of the compiled server: listens on `127.0.0.1:\$PORT` (default 8000).
With `CUTS_WATCH_STDIN=1` (set by the desktop app) it exits when stdin closes.
"""
function julia_main()::Cint
    try
        get(ENV, "CUTS_WATCH_STDIN", "") == "1" && exit_when_stdin_closes()
        serve(; port = parse(Int, get(ENV, "PORT", "8000")))
    catch err
        err isa InterruptException && return 0
        Base.invokelatest(Base.display_error, err, catch_backtrace())
        return 1
    end
    return 0
end

end # module
