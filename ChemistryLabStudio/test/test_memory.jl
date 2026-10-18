# The memory figures: positive and bounded by the machine; the cgroup readers
# on files written for the test; the heap hint within its bounds.

avail = CLS.available_memory()
@test 0 < avail <= Sys.total_physical_memory()
# read on every system the studio supports
@test CLS.process_rss(getpid()) > 0
@test CLS.process_rss(typemax(Int32)) === nothing

@test CLS.heap_size_hint(0) == "1024M"
@test CLS.heap_size_hint(100 * CLS.GIB) == string(60 * 1024, "M")     # no upper cap
@test CLS.heap_size_hint(3 * CLS.GIB) == string(6 * (3 * CLS.GIB ÷ 10) ÷ 2^20, "M")
@test CLS.gigabytes(3 * CLS.GIB ÷ 2) == 1.5

mktempdir() do dir
    write(joinpath(dir, "memory.max"), "max\n")
    write(joinpath(dir, "memory.high"), "1073741824\n")
    write(joinpath(dir, "memory.stat"), "anon 10\nactive_file 100\ninactive_file 23\nshmem 5\n")
    @test CLS.cgroup_value(dir, "memory.max") == typemax(Int)
    @test CLS.cgroup_value(dir, "memory.high") == CLS.GIB
    @test CLS.cgroup_value(dir, "absent") == typemax(Int)
    @test CLS.page_cache(dir) == 123
    @test CLS.page_cache(joinpath(dir, "nowhere")) == 0
end

b = CLS.cgroup_budget()
@test b === nothing || (b.limit > 0 && b.room <= b.limit)

# the compute process is measured in the coverage mode of this process
path = "/some/path"
GC.@preserve path begin
    tracked = pointer(path)
    @test CLS.coverage_flags(0, tracked) == String[]
    @test CLS.coverage_flags(1, tracked) == ["--code-coverage=user"]
    @test CLS.coverage_flags(2, tracked) == ["--code-coverage=all"]
    @test CLS.coverage_flags(3, tracked) == ["--code-coverage=@/some/path"]
end
