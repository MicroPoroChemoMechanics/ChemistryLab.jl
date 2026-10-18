# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Memory the machine can spare, read the way the kernel enforces it. A cgroup
# quota (a per-user slice capped below the physical memory, say) binds before
# the machine runs out, and libuv's `Sys.free_memory` sees only the cgroup a
# process sits in, never its ancestors nor `memory.high`. A calculation is
# refused below a floor before it starts, rather than discovered too big during
# the solve, and the compute process is stopped if it outgrows what was left.

const GIB = 2^30

"The least room a calculation starts with, in bytes."
const RUN_FLOOR = 3 * GIB ÷ 2

"What the compute process must leave to the rest of the machine, in bytes."
const HOST_MARGIN = GIB ÷ 2

"""
    cgroup_value(dir, file) -> Int

A cgroup (v2) file read as an integer: `typemax(Int)` for `max`, or for a file
that does not exist (no such controller, not Linux).
"""
function cgroup_value(dir::AbstractString, file::AbstractString)
    path = joinpath(dir, file)
    isfile(path) || return typemax(Int)
    return something(tryparse(Int, strip(read(path, String))), typemax(Int))
end

"Page cache charged to a cgroup, in bytes: the kernel reclaims it before a limit bites."
function page_cache(dir::AbstractString)
    path = joinpath(dir, "memory.stat")
    isfile(path) || return 0
    cache = 0
    for line in eachline(path)
        key, value = split(line)
        key in ("active_file", "inactive_file") && (cache += parse(Int, value))
    end
    return cache
end

"""
    cgroup_budget() -> Union{Nothing, NamedTuple}

The memory budget of the cgroups (v2) this process runs in, from its own up to
the root: `limit`, the tightest `min(memory.high, memory.max)`, and `room`, the
least left under one of them, page cache counted as free. `nothing` when no
limit applies or none can be read (not Linux, cgroup v1).
"""
function cgroup_budget()
    Sys.islinux() || return nothing
    isfile("/proc/self/cgroup") || return nothing
    lines = readlines("/proc/self/cgroup")
    k = findfirst(startswith("0::"), lines)
    k === nothing && return nothing
    dir = rstrip("/sys/fs/cgroup" * lines[k][4:end], '/')
    limit, room = typemax(Int), typemax(Int)
    while startswith(dir, "/sys/fs/cgroup/")
        cap = min(cgroup_value(dir, "memory.high"), cgroup_value(dir, "memory.max"))
        if cap < typemax(Int)
            limit = min(limit, cap)
            room = min(room, cap - cgroup_value(dir, "memory.current") + page_cache(dir))
        end
        dir = dirname(dir)
    end
    return limit == typemax(Int) ? nothing : (; limit, room)
end

"""
    available_memory() -> Int

The bytes a new process may still take: the lower of the machine's free memory
and of the room the tightest enclosing cgroup leaves.
"""
function available_memory()
    machine = Int(Sys.free_physical_memory())
    budget = cgroup_budget()
    return budget === nothing ? min(machine, Int(Sys.free_memory())) : min(machine, budget.room)
end

"""
    process_rss(pid) -> Union{Nothing, Int}

The resident memory of process `pid`, in bytes; `nothing` where it cannot be
read (a process that has ended, or a system without `/proc`, where the
watchdog relies on the floor checked before each run).
"""
function process_rss(pid::Integer)
    lines = try
        readlines("/proc/$pid/status")
    catch err
        # no such process (it has ended), or no `/proc` on this system
        err isa SystemError || rethrow()
        return nothing
    end
    k = findfirst(startswith("VmRSS:"), lines)
    return k === nothing ? nothing : 1024 * parse(Int, split(lines[k])[2])
end

"""
    heap_size_hint(available) -> String

The `--heap-size-hint` given to the compute process: six tenths of what is
available, between 1 and 4 GiB. The garbage collector works harder near it,
which keeps the process's footprint near what was measured, about 1.7 GiB for a
cement paste.
"""
heap_size_hint(available::Integer) = string(clamp(6 * (available ÷ 10), GIB, 4GIB) ÷ 2^20, "M")

"Memory in gigabytes, for messages."
gigabytes(bytes::Integer) = round(bytes / GIB; digits = 1)
