using HTTP
using JSON3
using LinearAlgebra
using DifferentialEquations
using NLsolve
using Sockets

const TRI_HEIGHT = 0.86602540378
const DEFAULTS_IMPRODUCTIVO = Dict(
    "vN" => 7.0,
    "vM" => 2.0,
    "vr" => 3.0,
    "vsigma" => 1.0,
    "vd" => 0.5,
    "vu" => 1.0,
    "vpc" => -1.0
)
const DEFAULTS_PRODUCTIVO = Dict(
    "vN" => 7.0,
    "vM" => 3.0,
    "vr" => 3.0,
    "vsigma" => 1.0,
    "vd" => 0.6,
    "vu" => 1.0,
    "vpc" => -0.5
)

# Variable order: yy[1]=Cooperators (x), yy[2]=Defectors (y), yy[3]=Loners (z)
function safe_afu2improductivo!(dy, yy, p, t)
    x, y, z = yy
    z_safe = clamp(z, 0.0, 1.0)
    N, c, r, sigma, d, u, pc = p
    n = N - c
    Q = pc * u * d * r

    if abs(1 - z_safe) < 1e-12
        a_inner = 1.0
        term_x_over_1minz_times_1mina = x * (n - 1) / 2
    else
        a_inner = (1 - z_safe^n) / (n * (1 - z_safe))
        term_x_over_1minz_times_1mina = (x / (1 - z_safe)) * (1 - a_inner)
    end

    py = sigma * z_safe^(max(0, n - 1)) + (1 - d * u) * r * term_x_over_1minz_times_1mina + c * a_inner
    px = sigma * z_safe^(max(0, n - 1)) + (r - Q + c) * a_inner + (r - Q) * term_x_over_1minz_times_1mina +
         (1 - r - c) * z_safe^(max(0, n - 1)) + Q - 1
    pz = sigma
    pprom = x * px + y * py + z * pz

    dy[1] = x * (px - pprom)
    dy[2] = y * (py - pprom)
    dy[3] = z * (pz - pprom)
    return nothing
end

function safe_afu2productivo!(dy, yy, p, t)
    x, y, z = yy
    z_safe = clamp(z, 0.0, 1.0)
    N, c, r, sigma, d, u, pc = p
    n = N - c
    Q = pc * u * d * r

    if abs(1 - z_safe) < 1e-12
        a_inner = 1.0
        term_x_over_1minz_times_1mina = x * (n - 1) / 2
    else
        a_inner = (1 - z_safe^n) / (n * (1 - z_safe))
        term_x_over_1minz_times_1mina = (x / (1 - z_safe)) * (1 - a_inner)
    end

    py = sigma * z_safe^(max(0, n - 1)) + (1 - d * u) * r * term_x_over_1minz_times_1mina +
         c * (1 - d * u) * a_inner
    px = sigma * z_safe^(max(0, n - 1)) + (r - Q + c) * a_inner + (r - Q) * term_x_over_1minz_times_1mina +
         (1 - r - c) * z_safe^(max(0, n - 1)) + Q - 1 + c * x * (n - 1) * z_safe^(max(0, n - 2))
    pz = sigma
    pprom = x * px + y * py + z * pz

    dy[1] = x * (px - pprom)
    dy[2] = y * (py - pprom)
    dy[3] = z * (pz - pprom)
    return nothing
end

function find_equilibria(p, dyn!)
    equilibria = Vector{Vector{Float64}}()
    starts = [
        [0.33, 0.33, 0.33],
        [0.99, 0.005, 0.005], [0.005, 0.99, 0.005], [0.005, 0.005, 0.99],
        [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
        [0.1, 0.1, 0.8], [0.8, 0.1, 0.1], [0.1, 0.8, 0.1], [0.5, 0.2, 0.3]
    ]

    for u0 in starts
        function f!(F, x_vec)
            dyn!(F, x_vec, p, 0.0)
            F[3] = sum(x_vec) - 1.0
        end
        res = nlsolve(f!, u0, autodiff=:forward)
        if converged(res)
            sol_x = res.zero
            if all(sol_x .>= -1e-5) && all(sol_x .<= 1.00001)
                sol_x = max.(0.0, min.(1.0, sol_x))
                sol_x ./= sum(sol_x)
                if isempty(equilibria) || !any(e -> norm(e - sol_x) < 1e-4, equilibria)
                    push!(equilibria, sol_x)
                end
            end
        end
    end
    return equilibria
end

function trajectories_data(p, dyn!; tspan=(0.0, 20.0), step=0.05, saveat=0.2)
    v_ini = Vector{Vector{Float64}}()
    for _x in 0:step:1, _y in 0:step:(1 - _x)
        _z = 1 - _x - _y
        if _z < 0.99
            push!(v_ini, [_x, _y, _z])
        end
    end

    trajectories = Vector{Dict{String, Vector{Float64}}}()
    for u0 in v_ini
        prob = ODEProblem(dyn!, u0, tspan, p)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=saveat)
        x_traj = [u[2] + 0.5 * u[1] for u in sol.u]
        y_traj = [TRI_HEIGHT * u[1] for u in sol.u]
        push!(trajectories, Dict("x" => x_traj, "y" => y_traj))
    end
    return trajectories
end

function vector_field_data(p, dyn!; step=0.04)
    vf_x = Float64[]
    vf_y = Float64[]
    for _x in 0.0:step:1.0, _y in 0.0:step:(1.0 - _x)
        _z = 1.0 - _x - _y
        if _z >= -1e-10
            dy = zeros(3)
            dyn!(dy, [_x, _y, _z], p, 0.0)
            len = sqrt(dy[1]^2 + dy[2]^2 + dy[3]^2)
            if len > 1e-4
                vis = min(len * 5.0, 0.04)
                ndx, ndy = dy[1] / len * vis, dy[2] / len * vis

                sx = _y + 0.5 * _x
                sy = TRI_HEIGHT * _x
                ex = (_y + ndy) + 0.5 * (_x + ndx)
                ey = TRI_HEIGHT * (_x + ndx)

                vx = ex - sx
                vy = ey - sy
                vmag = sqrt(vx^2 + vy^2)
                if vmag > 1e-6
                    angle = atan(vy, vx)
                    head_len = vmag * 0.35
                    recoil = 2.61799
                    rx = ex + head_len * cos(angle + recoil)
                    ry = ey + head_len * sin(angle + recoil)
                    lx = ex + head_len * cos(angle - recoil)
                    ly = ey + head_len * sin(angle - recoil)
                    append!(vf_x, [sx, ex, rx, ex, lx, NaN])
                    append!(vf_y, [sy, ey, ry, ey, ly, NaN])
                end
            end
        end
    end
    return Dict("x" => to_nullable(vf_x), "y" => to_nullable(vf_y))
end

function stability_grid(p, dyn!; grid_size=120)
    x_grid = collect(range(0, 1, length=grid_size))
    y_grid = collect(range(0, TRI_HEIGHT, length=grid_size))
    z_grid = fill(NaN, grid_size, grid_size)
    h = 1e-6

    function get_f(va, vc)
        vx = va
        vy = vc
        vz = 1.0 - vx - vy
        if vz < 0
            vz = 0.0
        end
        dy_out = zeros(3)
        dyn!(dy_out, [vx, vy, vz], p, 0.0)
        return dy_out[1], dy_out[2]
    end

    for (i, Y) in enumerate(y_grid), (j, X) in enumerate(x_grid)
        a = Y / TRI_HEIGHT
        c = X - 0.5 * a
        b = 1.0 - a - c
        if a >= -1e-8 && b >= -1e-8 && c >= -1e-8
            s = a + b + c
            if abs(s - 1.0) > 1e-8
                a /= s
                b /= s
                c /= s
            end
            if a < 0 || b < 0 || c < 0
                z_grid[i, j] = NaN
                continue
            end

            fa, fc = get_f(a, c)
            fa_da, fc_da = get_f(a + h, c)
            fa_dc, fc_dc = get_f(a, c + h)
            J = [(fa_da - fa) / h (fa_dc - fa) / h; (fc_da - fc) / h (fc_dc - fc) / h]
            z_grid[i, j] = maximum(real(eigen(J).values))
        else
            z_grid[i, j] = NaN
        end
    end

    return x_grid, y_grid, z_grid
end

function simplex_boundary()
    return Dict("x" => [0.0, 1.0, 0.5, 0.0], "y" => [0.0, 0.0, TRI_HEIGHT, 0.0])
end

function equilibria_xy(equilibria)
    eq_x = [e[2] + 0.5 * e[1] for e in equilibria]
    eq_y = [TRI_HEIGHT * e[1] for e in equilibria]
    return Dict("x" => eq_x, "y" => eq_y)
end

function to_nested(z_grid)
    rows = Vector{Vector{Union{Nothing, Float64}}}(undef, size(z_grid, 1))
    for i in 1:size(z_grid, 1)
        row = Vector{Union{Nothing, Float64}}(undef, size(z_grid, 2))
        for j in 1:size(z_grid, 2)
            val = z_grid[i, j]
            row[j] = isfinite(val) ? val : nothing
        end
        rows[i] = row
    end
    return rows
end

function to_nullable(v::Vector{Float64})
    out = Vector{Union{Nothing, Float64}}(undef, length(v))
    for i in eachindex(v)
        out[i] = isfinite(v[i]) ? v[i] : nothing
    end
    return out
end

function normalize_mode(value)
    if value === nothing
        return "improductivo"
    end
    mode = lowercase(String(value))
    return mode == "productivo" ? "productivo" : "improductivo"
end

function get_body_value(body, key)
    if body isa AbstractDict
        if haskey(body, key)
            return body[key]
        end
        sym = Symbol(key)
        if haskey(body, sym)
            return body[sym]
        end
    end
    return nothing
end

function params_from_body(body, defaults)
    merged = copy(defaults)
    for (k, v) in pairs(body)
        key = String(k)
        if haskey(merged, key)
            merged[key] = float(v)
        end
    end
    return [
        merged["vN"],
        merged["vM"],
        merged["vr"],
        merged["vsigma"],
        merged["vd"],
        merged["vu"],
        merged["vpc"]
    ]
end

function compute_payload(p, dyn!)
    equilibria = find_equilibria(p, dyn!)
    boundary = simplex_boundary()
    eq_xy = equilibria_xy(equilibria)
    trajectories = trajectories_data(p, dyn!)
    vector_field = vector_field_data(p, dyn!)
    x_grid, y_grid, z_grid = stability_grid(p, dyn!)

    payload = Dict{String, Any}()
    payload["simplex"] = Dict(
        "boundary" => boundary,
        "trajectories" => trajectories,
        "equilibria" => eq_xy,
        "vector_field" => vector_field
    )
    payload["heatmap"] = Dict(
        "x_grid" => x_grid,
        "y_grid" => y_grid,
        "z_grid" => to_nested(z_grid),
        "boundary" => boundary,
        "equilibria" => eq_xy
    )
    return payload
end

function json_response(status, payload)
    headers = [
        "Content-Type" => "application/json",
        "Access-Control-Allow-Origin" => "*",
        "Access-Control-Allow-Headers" => "Content-Type",
        "Access-Control-Allow-Methods" => "POST, OPTIONS, GET"
    ]
    return HTTP.Response(status, headers, JSON3.write(payload))
end

function handler(req)
    path = HTTP.URI(req.target).path
    if req.method == "OPTIONS"
        return json_response(200, Dict("ok" => true))
    elseif req.method == "GET" && path == "/api/health"
        return json_response(200, Dict("status" => "ok"))
    elseif req.method == "POST" && path == "/api/data"
        body = isempty(req.body) ? Dict{String, Any}() : JSON3.read(String(req.body))
        mode = normalize_mode(get_body_value(body, "mode"))
        defaults = mode == "productivo" ? DEFAULTS_PRODUCTIVO : DEFAULTS_IMPRODUCTIVO
        dyn! = mode == "productivo" ? safe_afu2productivo! : safe_afu2improductivo!
        p = params_from_body(body, defaults)
        payload = compute_payload(p, dyn!)
        payload["mode"] = mode
        return json_response(200, payload)
    else
        return json_response(404, Dict("error" => "Not found"))
    end
end

function main()
    port = parse(Int, get(ENV, "PORT", "8000"))
    println("Julia backend running on http://0.0.0.0:$port")
    HTTP.serve(handler, ip"0.0.0.0", port)
end

main()
