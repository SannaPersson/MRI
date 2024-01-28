using KomaMRI

# Define simulation parameters dictionary (sim)
# This should be filled with the translated parameters from your Pulseq sequence
sim_params = Dict(
    "fov" => 256e-3,
    "Nx" => 96,
    "Ny" => 96,
    # ... Add other translated parameters here
)

# Define reconstruction parameters dictionary (rec)
# Fill this with any specific reconstruction parameters you have
rec_params = Dict(
    :param1 => value1,
    :param2 => value2,
    # ... Add other reconstruction parameters here
)

# Launch KomaUI with specified settings
KomaUI(
    darkmode=true,           # Enable dark mode
    frame=true,              # Display the upper frame
    phantom_mode="2D",       # Load a 2D brain phantom
    sim=sim_params,          # Pass the simulation parameters
    rec=rec_params,          # Pass the reconstruction parameters
    return_window=false,     # Do not return the Blink window object
    show_window=true         # Display the Blink window
)
