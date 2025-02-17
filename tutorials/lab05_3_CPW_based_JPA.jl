# original code from https://github.com/kpobrien/JosephsonCircuits.jl
# circuit parameters from T Yamamoto et al., "Flux-driven Josephson parametric amplifier," Appl. Phys. Lett., 93, 042510 (2008)
# find the online link: https://doi.org/10.1063/1.2964182
# 3-Wave Mixing

using JosephsonCircuits
using Plots

@variables R Cc Cj Lj Cr Lr Ll Ldc K Lg
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    # a very large inductor so the DC node flux of this node isn't floating
    ("L0","1","0",Lg), 
    ("C1","1","2",Cc),
    ("L1","2","3",Lr),
    ("C2","2","0",Cr),
    ("Lj1","3","0",Lj),
    ("Cj1","3","0",Cj),
    ("L2","3","4",Ll),
    ("Lj2","4","0",Lj),
    ("Cj2","4","0",Cj),
    ("L3","5","0",Ldc), 
    ("K1","L2","L3",K),
    # a port with a very large resistor so we can apply the bias across the port
    ("P2","5","0",2),
    ("R2","5","0",1000.0),
] 

circuitdefs = Dict(
    Lj =>219.63e-12,
    Lr =>0.4264e-9,
    Lg =>100.0e-9,
    Cc => 16.0e-15,
    Cj => 10.0e-15, 
    Cr => 0.4e-12,
    R => 50.0, 
    Ll => 34e-12, 
    K => 0.999, # the inverse inductance matrix for K=1.0 diverges, so set K<1.0
    Ldc => 0.74e-12,
)

ws = 2*pi*(9.7:0.0001:9.8)*1e9
wp = (2*pi*19.50*1e9,)
Ip = 0.74e-6
Idc = 140.3e-6

# add the DC bias and pump to port 2
sourcespumpon = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=Ip)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)
@time jpapumpon = hbsolve(ws, wp, sourcespumpon, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing

plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="JosephsonCircuits.jl",
)