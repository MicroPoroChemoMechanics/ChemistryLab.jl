# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The routes of ChemistryLab Studio: the model the page edits, the script it
# stands for, running it, the files it is saved to and read from, the catalog
# of ChemistryLab and the library of published formulations.

"""
    StudioState()

The state of a session: the model, the file the script is bound to (`nothing`
until saved, and after opening a script edited by hand), and the catalog read
from ChemistryLab once the compute process is up.
"""
mutable struct StudioState
    model::StudioModel
    path::Union{Nothing, String}
    catalog::Union{Nothing, Dict{String, Any}}
    const lock::ReentrantLock
end

StudioState() = StudioState(default_model(), nothing, nothing, ReentrantLock())

"What the page needs after any change: the model, its script, its problems, the bound file."
function state_payload(st::StudioState)
    m, path = lock(() -> (st.model, st.path), st.lock)
    return Dict{String, Any}(
        "model" => to_dict(m), "script" => generate(m).source,
        "problems" => to_dict.(problems(m)), "path" => path, "stem" => script_stem(m.title),
    )
end

bad_request(text::AbstractString) = json_response(400, Dict("error" => text))

function set_model!(st::StudioState, m::StudioModel; path::Union{Nothing, String} = nothing, bind::Bool = false)
    lock(st.lock) do
        st.model = m
        bind && (st.path = path)
    end
    return nothing
end

function post_model(st::StudioState, req::HTTP.Request)
    body = read_json(req)
    raw = get(body, "model", nothing)
    raw isa AbstractDict || return bad_request("the request carries no model")
    m = try
        model_from_dict(raw)
    catch err
        err isa ArgumentError || rethrow()
        return bad_request(err.msg)
    end
    set_model!(st, m)
    return json_response(200, state_payload(st))
end

function post_run(session::Session, st::StudioState)
    m = lock(() -> st.model, st.lock)
    ps = problems(m)
    isempty(ps) || return json_response(422, Dict("ok" => false, "problems" => to_dict.(ps)))
    script = generate(m)
    outcome = compute!(session, :(Main.StudioWorker.run_compute($(script.compute))))
    outcome.ok || return json_response(200, Dict("ok" => false, "error" => outcome.error))
    return json_response(200, outcome.value)
end

"""
    script_path(raw) -> String

The file a script is saved to: the path the page gives, `~` expanded, with the
`.jl` extension added when it is missing.
"""
function script_path(raw::AbstractString)
    path = abspath(expanduser(strip(raw)))
    endswith(path, ".jl") || (path *= ".jl")
    return path
end

function post_save(st::StudioState, req::HTTP.Request)
    body = read_json(req)
    bound = lock(() -> st.path, st.lock)
    raw = get(body, "path", bound)
    raw isa AbstractString || return bad_request("choose a file name to save the script to")
    path = script_path(raw)
    isdir(dirname(path)) || return bad_request("the folder $(dirname(path)) does not exist")
    overwrite = get(body, "overwrite", false) === true
    (isfile(path) && path != bound && !overwrite) && return json_response(409, Dict("exists" => true, "path" => path))
    source = lock(() -> generate(st.model).source, st.lock)
    # written next to the target and moved, so that a failure leaves the old file whole
    tmp = path * ".part"
    write(tmp, source)
    mv(tmp, path; force = true)
    lock(() -> (st.path = path), st.lock)
    return json_response(200, state_payload(st))
end

const MAX_SCRIPT_BYTES = 5 * 2^20

function post_open(st::StudioState, req::HTTP.Request)
    raw = get(read_json(req), "path", nothing)
    raw isa AbstractString || return bad_request("choose a file to open")
    path = abspath(expanduser(raw))
    isfile(path) || return bad_request("no file $path")
    filesize(path) <= MAX_SCRIPT_BYTES || return bad_request("$path is too large to be a studio script")
    rb = try
        read_script(read(path, String))
    catch err
        # a damaged model block: say so, and leave the current model and file alone
        err isa ArgumentError || err isa ErrorException || rethrow()
        return json_response(422, Dict("error" => "The model block of $path cannot be read: $(sprint(showerror, err))"))
    end
    rb.model === nothing && return json_response(422, Dict("error" => rb.message))
    set_model!(st, rb.model; path = rb.edited ? nothing : path, bind = true)
    payload = state_payload(st)
    payload["message"] = rb.message
    return json_response(200, payload)
end

"""
    browse(path) -> Dict{String, Any}

The folders and `.jl` files of the folder `path` (the home folder when empty),
for the page's file picker: the browser cannot give the server a path, and a
remote session has no local file dialog anyway.
"""
function browse(path::AbstractString)
    dir = isempty(path) ? homedir() : abspath(expanduser(path))
    isdir(dir) || (dir = homedir())
    folders, files = String[], String[]
    names = try
        readdir(dir)
    catch err
        # a folder the user may not read lists as empty, its parent still reachable
        err isa Base.IOError || rethrow()
        String[]
    end
    for name in names
        startswith(name, '.') && continue
        full = joinpath(dir, name)
        # a dangling link or a folder the user may not read is left out
        isdir(full) ? push!(folders, name) : (endswith(name, ".jl") && isfile(full) && push!(files, name))
    end
    return Dict{String, Any}(
        "path" => dir, "parent" => dirname(dir), "folders" => folders, "files" => files,
        "places" => [Dict("label" => "Home", "path" => homedir()), Dict("label" => "Current folder", "path" => pwd())],
    )
end

"The library of published formulations, as `library/index.toml` lists it."
function library_entries()
    return get(TOML.parsefile(joinpath(LIBRARY_DIR, "index.toml")), "entry", Any[])
end

function library_listing()
    return [
        Dict{String, Any}(k => e[k] for k in ("id", "title", "designation", "process", "summary", "sources", "level") if haskey(e, k))
            for e in library_entries()
    ]
end

function post_library(st::StudioState, req::HTTP.Request)
    id = get(read_json(req), "id", nothing)
    k = findfirst(e -> get(e, "id", nothing) == id, library_entries())
    k === nothing && return bad_request("no library entry $(repr(id))")
    set_model!(st, model_from_dict(library_entries()[k]["model"]); path = nothing, bind = true)
    return json_response(200, state_payload(st))
end

"""
    chemistrylab_routes!(router, session)

Register the routes of ChemistryLab Studio.
"""
function chemistrylab_routes!(router::HTTP.Router, session::Session)
    st = session.state
    HTTP.register!(router, "GET", "/api/state", _ -> json_response(200, state_payload(st)))
    HTTP.register!(router, "POST", "/api/model", req -> post_model(st, req))
    HTTP.register!(
        router, "GET", "/api/catalog",
        _ -> json_response(200, something(lock(() -> st.catalog, st.lock), Dict("pending" => true))),
    )
    HTTP.register!(router, "POST", "/api/run", _ -> post_run(session, st))
    HTTP.register!(router, "POST", "/api/save", req -> post_save(st, req))
    HTTP.register!(router, "POST", "/api/open", req -> post_open(st, req))
    HTTP.register!(router, "GET", "/api/browse", req -> json_response(200, browse(get(query(req), "path", ""))))
    HTTP.register!(router, "GET", "/api/library", _ -> json_response(200, library_listing()))
    HTTP.register!(router, "POST", "/api/library", req -> post_library(st, req))
    return nothing
end

"""
    warm_up!(session)

Read the catalog of ChemistryLab, then run the default calculation once so that
its compilation, about two minutes, is paid before the user's first run.
"""
function warm_up!(session::Session)
    st, cp = session.state, session.compute
    set_state!(cp, :ready, "Reading the catalog of ChemistryLab.")
    catalog = call!(cp, :(Main.StudioWorker.catalog()))
    lock(st.lock) do
        st.catalog = catalog.ok ? catalog.value : Dict{String, Any}("error" => catalog.error)
    end
    catalog.stopped === :none || return nothing
    set_state!(cp, :ready, "Compiling the first calculation: about two minutes, once per session.")
    warm = call!(cp, :(Main.StudioWorker.run_compute($(generate(default_model()).compute))))
    warm.stopped === :none || return nothing
    set_state!(cp, :ready)
    return nothing
end
