using Rhapsodie
using DelimitedFiles
using EasyFITS

#include("test_separable_reconstruction.jl")

# contrast_list = [i for i in range(-1.5, 0, step=0.5)]
contrast_list = [-2.0]
max_iter = 700
# α=10^-5
α=2.044595476128458
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

parameter_type = "mixed"

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

# mse_list = Vector{Float64}()
for k in contrast_list
    println("------Contrast parameter: ", k, "------")
    true_polar_map = Rhapsodie.read_and_fill_polar_map(parameter_type, "test_results/contrast_10e$(k)/TRUE.fits")
    diff_fits = EasyFITS.readfits("test_results/contrast_10e$(k)/Results_Separable_DoubleDifference.fits")
    diff_polar_map = Rhapsodie.TPolarimetricMap("mixed", diff_fits[ :, :, 1]', diff_fits[ :, :, 5]', diff_fits[:, :, 2]', diff_fits[:, :, 3]')

    Rhapsodie.load_data("test_results/contrast_10e$(k)/DATA.fits", "test_results/contrast_10e$(k)/WEIGHT.fits")

    PSF = readfits("data_for_demo/PSF_parametered_Airy.fits");
    A = set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);
    # X0 = TPolarimetricMap(parameter_type, zeros(Rhapsodie.get_par().cols));
    X0 = diff_polar_map;
    regularisation_parameters = 10 .^[0,  -1. , -1, -0.66] #(in log10) star, disk
    regularisation_parameters[1] = 0
    regularisation_parameter_list = [-2.8123070745690395]

    regularisation_parameters[4] = regularisation_parameter_list[1]
    x = apply_rhapsodie(X0, A, Rhapsodie.dataset, regularisation_parameters, α=α,
                        maxeval=1000, maxiter=max_iter);
    crop!(x)
    write_polar_map(x, "test_results/contrast_10e$(k)/rhapsodie_method_results/max_iter_$(max_iter)/RHAPSODIE_opti_params_$(parameter_type)_Iu_Ip_separate_from_dd.fits", overwrite=true)
    # append!(mse_list, Rhapsodie.MSE_object(x, true_polar_map))
    empty!(Rhapsodie.dataset)
end
# writedlm("test_results/mse_list.txt", mse_list)
