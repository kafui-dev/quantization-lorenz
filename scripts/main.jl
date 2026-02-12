include("required-packages.jl")

# -------------------------------------------------------------------
# 1. Lorenz system parameters and derivative
# -------------------------------------------------------------------
σ = 10.0
ρ = 28.0
β = 8/3

function lorenz(u, p, t)
    x, y, z = u
    dx = σ * (y - x)
    dy = x * (ρ - z) - y
    dz = x * y - β * z
    return [dx, dy, dz]
end

# -------------------------------------------------------------------
# 2. Runge–Kutta 4 integrator
# -------------------------------------------------------------------
function rk4_step(f, u, t, dt, p=nothing)
    k1 = f(u, p, t)
    k2 = f(u .+ 0.5*dt*k1, p, t+0.5*dt)
    k3 = f(u .+ 0.5*dt*k2, p, t+0.5*dt)
    k4 = f(u .+ dt*k3, p, t+dt)
    return u .+ (dt/6)*(k1 .+ 2k2 .+ 2k3 .+ k4)
end

# -------------------------------------------------------------------
# 3. Simulation settings
# -------------------------------------------------------------------
dt = 0.01
nsteps = 10000
u0 = [1.0, 0.0, 0.0]          # initial condition (x, y, z)
t0 = 0.0

# Pre‑allocate trajectory array
traj = zeros(nsteps+1, 3)
traj[1, :] = u0

# -------------------------------------------------------------------
# 4. Time integration
# -------------------------------------------------------------------
u = copy(u0)
t = t0
for i in 1:nsteps
    u = rk4_step(lorenz, u, t, dt, nothing)
    traj[i+1, :] = u
    t += dt
end

# -------------------------------------------------------------------
# 5. Trajectory components and time vector
# -------------------------------------------------------------------
x, y, z = traj[:, 1], traj[:, 2], traj[:, 3]
time = (0:nsteps) * dt

# -------------------------------------------------------------------
# 6. Quantization of the trajectory
# -------------------------------------------------------------------
# 10^14 magnification of x, y and z
x_mag = round.(Int64, x .* 1e14)
y_mag = round.(Int64, y .* 1e14)
z_mag = round.(Int64, z .* 1e14)

# 8 LSBs extraction
x_lsb = x_mag .& 0xFF
y_lsb = y_mag .& 0xFF
z_lsb = z_mag .& 0xFF

# XOR combination → final quantized sequence (8‑bit values, 0–255)
quantized_sequence = x_lsb .⊻ y_lsb .⊻ z_lsb

# -------------------------------------------------------------------
# 7. Convert quantized integers to a binary bitstream
# -------------------------------------------------------------------
function bits_from_ints(ints; nbits=8)
    bits = Bool[]
    for val in ints
        for b in 0:nbits-1
            push!(bits, ((val >> b) & 1) == 1)
        end
    end
    return bits
end

bitstream = bits_from_ints(quantized_sequence)

# -------------------------------------------------------------------
# 8. PLOT 1 – Waveform Transformation (Time‑Series)
# -------------------------------------------------------------------
p1_raw = plot(time, x, 
              title = "Raw x(t) – Lorenz waveform",
              xlabel = "Time", ylabel = "x(t)",
              linecolor = :blue, label = "x(t)")

p1_quant = plot(time, quantized_sequence,
                title = "Quantized output – 8‑bit integer sequence",
                xlabel = "Time", ylabel = "Quantized value (0–255)",
                linecolor = :red, label = "Q(t)",
                ylims = (0, 255))

waveform = plot(p1_raw, p1_quant, layout = (2,1), size = (800, 600))

# -------------------------------------------------------------------
# 9. PLOT 2 – Phase Space vs. Bitstream Raster ("Snow" Test)
# -------------------------------------------------------------------
# 3D phase space (Lorenz attractor)
p2_phase = plot(x, y, z,
                title = "Lorenz attractor – Phase space",
                xlabel = "x", ylabel = "y", zlabel = "z",
                linewidth = 0.5, linecolor = :viridis,
                label = "", legend = false)

# Raster image from the bitstream (first 256×256 = 65536 bits)
raster_size = 256
total_bits_needed = raster_size * raster_size
bit_slice = bitstream[1:total_bits_needed]
raster_matrix = reshape(bit_slice, (raster_size, raster_size))'

p2_raster = heatmap(raster_matrix,
                    title = "Bitstream raster ($raster_size×$raster_size) – should be white noise",
                    xlabel = "Column", ylabel = "Row",
                    color = :grays, clim = (0, 1),
                    aspect_ratio = :equal, legend = false)

phases = plot(p2_phase, p2_raster, layout = (1,2), size = (1200, 500))

# -------------------------------------------------------------------
# 10. PLOT 3 – Autocorrelation Function (ACF) Transformation
# -------------------------------------------------------------------
maxlag = 200
# Raw x ACF
x_acf = autocor(x, 0:maxlag)
p3_raw = plot(0:maxlag, x_acf,
              title = "ACF of raw x(t) – slow decay, periodic humps",
              xlabel = "Lag", ylabel = "Autocorrelation",
              linecolor = :steelblue, marker = :circle, ms = 2,
              label = "x(t)")

# Bitstream ACF (binary)
bits_float = Float64.(bitstream)   # autocor expects real numbers
bits_acf = autocor(bits_float, 0:maxlag)
p3_bits = plot(0:maxlag, bits_acf,
               title = "ACF of bitstream – should be δ‑function (spike at 0)",
               xlabel = "Lag", ylabel = "Autocorrelation",
               linecolor = :crimson, marker = :circle, ms = 2,
               label = "bitstream")

acf = plot(p3_raw, p3_bits, layout = (2,1), size = (800, 600))

# -------------------------------------------------------------------
# 11. PLOT 4 – Return Maps (xₙ vs. xₙ₊₁)
# -------------------------------------------------------------------
# Raw return map
p4_raw = scatter(x[1:end-1], x[2:end],
                 title = "Return map – raw x",
                 xlabel = "xₙ", ylabel = "xₙ₊₁",
                 markersize = 0.8, markerstrokewidth = 0,
                 markercolor = :black, alpha = 0.5,
                 label = "", legend = false)

# Quantized return map (8‑bit integer blocks)
q = quantized_sequence
p4_quant = scatter(q[1:end-1], q[2:end],
                   title = "Return map – quantized (8‑bit blocks)",
                   xlabel = "Qₙ", ylabel = "Qₙ₊₁",
                   markersize = 0.8, markerstrokewidth = 0,
                   markercolor = :black, alpha = 0.3,
                   label = "", legend = false)

returnmap = plot(p4_raw, p4_quant, layout = (1,2), size = (1000, 400))

# -------------------------------------------------------------------
# 12. Save all figures
# -------------------------------------------------------------------
savefig(waveform, "output/plot1_waveform.png")
savefig(phases, "output/plot2_phase_raster.png")
savefig(acf, "output/plot3_acf.png")
savefig(returnmap, "output/plot4_returnmap.png")

println("All plots generated successfully.")