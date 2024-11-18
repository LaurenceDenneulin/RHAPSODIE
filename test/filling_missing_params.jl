using Rhapsodie
using CSV
using DataFrames
using DelimitedFiles
using EasyFITS

max_iter = 700
par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
# DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

mse_list = Vector{Float64}()
if length(ARGS) < 1
    println("Usage: julia script.jl <file_path>")
    exit(1)
end
file_path = ARGS[1]

df = CSV.read(file_path, DataFrame)
for row in eachrow(df)
    lambda = row[:λ]
    alpha = row[:α]
    contrast = row[:contrast]
    true_polar_map = Rhapsodie.read_and_fill_polar_map("mixed", "test_results/contrast_10e$(contrast)/TRUE.fits")
    Rhapsodie.load_data("test_results/contrast_10e$(contrast)/DATA.fits", "test_results/contrast_10e$(contrast)/WEIGHT.fits")

    PSF = readfits("data_for_demo/PSF_parametered_Airy.fits");
    A = set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);

    X0 = TPolarimetricMap("mixed", zeros(Rhapsodie.get_par().cols));
    regularisation_parameters = 10 .^[0,  -1. , -1, -0.66] #(in log10) star, disk
    regularisation_parameters[1] = 0
    # regularisation_parameter_list = [10^i for i in range(-3, 0, length=10)]
    regularisation_parameter_list = [10^lambda]
    regularisation_parameters[4] = regularisation_parameter_list[1]
    x = apply_rhapsodie(X0, A, Rhapsodie.dataset, regularisation_parameters, α=10^alpha,
                        maxeval=1000, maxiter=max_iter);
    crop!(x)
    res = Rhapsodie.MSE_object(x, x_true_polar_map)
    empty!(Rhapsodie.dataset)

    Iu_disk_mse = res[8]
    Ip_disk_mse = res[9]
    theta_mse = res[10]
    push!(mse_list, [lambda, alpha, k, Iu_disk_mse, Ip_disk_mse, theta_mse])
end
writedlm("test_results/missing_mse_list.csv", mse_list)
