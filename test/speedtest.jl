using SOLPS2IMAS: SOLPS2IMAS

println("-----------------------------------------------------------------------------")
b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
b2output = "$(@__DIR__)/../samples/b2time.nc"
b2mn = "$(@__DIR__)/../samples/b2mn.dat"
fort = (
    "$(@__DIR__)/../samples/fort.33",
    "$(@__DIR__)/../samples/fort.34",
    "$(@__DIR__)/../samples/fort.35")
print("solps2imas(b2gmtry, b2output; b2mn, fort) time with compilation: ")
@time dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output; b2mn, fort);
print("solps2imas(b2gmtry, b2output; b2mn, fort) time (true runtime): ")
@time dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output; b2mn, fort);
