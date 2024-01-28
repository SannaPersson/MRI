ENV["CUDA_VISIBLE_DEVICES"] = ""

using KomaMRI, MAT 


sys = Scanner(
    B0 = 3,                # Static magnetic field strength in Tesla
    Gmax = 22-3,            # Convert max_grad from mT/m to T/m by multiplying by 1e-3
    Smax = 160,              # Slew rate is already in T/m/s, no conversion needed
    ADC_Δt = 1e-7,           # ADC raster time in seconds (10 microseconds)
    seq_Δt = 1e-5,           # Block duration raster in seconds (100 microseconds)
    GR_Δt = 1e-5,            # Gradient raster time in seconds (100 microseconds)
    RF_Δt = 1e-6,            # RF raster time in seconds (1 microsecond)
    RF_ring_down_T = 30e-6,  # RF ringdown time in seconds
    RF_dead_time_T = 100e-6, # RF dead time in seconds
    ADC_dead_time_T = 10e-6  # ADC dead time in seconds
)

obj = brain_phantom2D() # a slice of a brain


seq_to_use = "four_shot.seq"
file_path = "/home/sanna/Desktop/Code/MRI/Project/" * seq_to_use
seq_file = joinpath(dirname(pathof(KomaMRI)), file_path)
seq = read_seq(seq_file)



simParams = KomaMRICore.default_sim_params()
simParams["return_type"] = "mat"
# no gpu
simParams["use_gpu"] = false


raw = simulate(obj, seq, sys; simParams, w=nothing)


# save raw data
#
# replace seq_to_use .mat for out_file
#
out_file = replace(seq_to_use, ".seq" => ".mat")
matwrite(out_file, Dict("raw" => raw))
matwrite("outputs/test.mat", Dict("raw" => raw))




