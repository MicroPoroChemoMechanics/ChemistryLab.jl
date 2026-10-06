# Delayed ettringite formation (docs/src/examples/delayed_ettringite.md): the
# SRPC of Lothenbach et al. (2008) at 5, 20, 50 and 80 °C and back, against the
# pore solutions of Lothenbach et al. (2007). The numbers the page states are
# pinned.

isdefined(@__MODULE__, :de_state) || include(joinpath(pkgdir(ChemistryLab), "scripts", "delayed_ettringite.jl"))

@testsection "durability: delayed ettringite formation" begin
    cs = l08t_system()
    paste = Dict(20.0 => de_state(20.0; cs))
    for T in (5.0, 50.0, 80.0)
        paste[T] = de_state(T; start = paste[20.0].state, cs)
    end
    o = Dict(T => de_observables(paste[T]) for T in keys(paste))
    ratio = [de_measured("SO4", T) / o[T].SO4 for T in (5.0, 20.0, 50.0)]
    @test minimum(ratio) ≈ 12.6 atol = 0.5
    @test maximum(ratio) ≈ 29.1 atol = 0.5
    @test o[50.0].SO4 / o[5.0].SO4 ≈ 35 atol = 1
    @test de_measured("SO4", 50.0) / de_measured("SO4", 5.0) ≈ 15.3 atol = 0.1
    @test o[20.0].OH / de_measured("OH-", 50.0) < 0.6
    @test o[50.0].monosulfate == 0
    @test o[80.0].monosulfate > 5
    @test o[80.0].ettringite / o[20.0].ettringite ≈ 0.68 atol = 0.01
    back = de_observables(de_state(20.0; start = paste[80.0].state, cs))
    @test back.ettringite ≈ o[20.0].ettringite rtol = 1.0e-6
    @test back.monosulfate < 1.0e-6
    @test temperature_range(L08T_DB["ettringite"]) == (273.15, 333.15)
end
