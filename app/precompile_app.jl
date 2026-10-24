# Warm-up workload run by PackageCompiler while building the desktop app.
# Everything executed here is compiled ahead of time, so the first request in
# the app is fast. It exercises every cut family, 2-D and 3-D geometry, LP and
# MPS parsing, error handling, and the HTTP routes.

using CutsApp
using CutsApp: HTTP, JSON3

const PROBLEMS = [
    """
    Minimize
     obj: x1 + 2 x2
    Subject To
     c1: -2 x1 - 2 x2 <= -0.5
     c2: 2 x1 - x2 <= 1.5
    Bounds
     x1 <= 1
     x2 <= 1
    General
     x1 x2
    End
    """,
    """
    Maximize
     obj: 5 x + 4 y + 3 z
    Subject To
     c1: 2 x + 3 y + z <= 5.5
     c2: 4 x + y + 2 z <= 11.3
     c3: 3 x + 4 y + 2 z <= 10.7
    General
     x y z
    End
    """,
]

port = 18_765
server = CutsApp.start(; port)
try
    HTTP.get("http://127.0.0.1:$port/")
    HTTP.get("http://127.0.0.1:$port/health")
    for text in PROBLEMS, norm in ("auto", "l1")
        body = JSON3.write(Dict("text" => text, "format" => "lp", "normalization" => norm,
                                "cuts" => ["gfc", "sic", "gmi", "vpc", "vpc_mono"], "var_index" => nothing))
        HTTP.post("http://127.0.0.1:$port/api/solve", ["Content-Type" => "application/json"], body)
    end
    # Solve tab: every algorithm on both problems.
    for text in PROBLEMS, algo in ("cutting_planes", "branch_and_bound", "branch_and_cut")
        body = JSON3.write(Dict("text" => text, "format" => "lp", "algorithm" => algo, "methods" => ["gmi", "vpc"],
                                "max_rounds" => 5, "split_rule" => "all", "root_cut_rounds" => 2, "node_cut_rounds" => 1))
        HTTP.post("http://127.0.0.1:$port/api/iterate", ["Content-Type" => "application/json"], body)
    end
    # An invalid problem exercises the error path (the server answers 400).
    HTTP.post("http://127.0.0.1:$port/api/solve", ["Content-Type" => "application/json"],
              JSON3.write(Dict("text" => "not a model", "format" => "lp")); status_exception = false)
finally
    close(server)
end
