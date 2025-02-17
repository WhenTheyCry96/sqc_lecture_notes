# original code from https://github.com/kpobrien/JosephsonCircuits.jl
# 4-Wave Mixing

using JosephsonCircuits
using Plots

@variables R Cc Lj Cj
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("Lj1","2","0",Lj),
    ("C2","2","0",Cj)]

circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    R => 50.0)

ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
Ip = 0.00565e-6
sources = [(mode=(1,),port=1,current=Ip)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)

@time jpa = hbsolve(ws, wp, sources, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs)

plot(
    jpa.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpa.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    label="JosephsonCircuits.jl",
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
)