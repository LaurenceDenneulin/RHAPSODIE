using Rhapsodie
using PRIMA
using DelimitedFiles
using EasyFITS

par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end
		
psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");
Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)


PSF = readfits("data_for_demo/PSF_parametered_Airy.fits");
A = set_fft_op(PSF[1:end÷2,:]'[:,:],psf_center[1:2]);
header = Vector{Vector{String}}()
push!(header, ["contrast", "λ", "α", "Iu_disk_mse", "Ip_disk_mse", "theta_mse"])

for k in range(-3, 0, step=0.5)
    println("Starting contrast: 10e$(k)")
    mse_list = Vector{Vector{Any}}()
    root_path = "test_results/contrast_10e$(k)/"
    dir_path = "test_results/prima/contrast_10e$(k)"
    if !isdir(dir_path)
        mkdir(dir_path)
    end
    io = open("test_results/prima/contrast_10e$(k)/mse.csv", "w")
    writedlm(io, header, ',')
    Rhapsodie.load_data("$(root_path)DATA.fits", "$(root_path)WEIGHT.fits")

        function calculate_MSE_for_prima(X::Vector{Float64})
            λ, α = X
            regularisation_parameters = 10 .^[0,  -1. , -1, λ]
            regularisation_parameters[1] = 0
            X0 = TPolarimetricMap("mixed", zeros(Rhapsodie.get_par().cols));
            x_est = apply_rhapsodie(X0, A, Rhapsodie.dataset, regularisation_parameters, α=10^α, maxeval=1000, maxiter=1000)
            true_polar_map = Rhapsodie.read_and_fill_polar_map("mixed", "$(root_path)TRUE.fits")
            crop!(x_est)
            write_polar_map(x_est, "test_results/prima/contrast_10e$(k)/RHAPSODIE_$(λ)_$(α).fits", overwrite=true)
            curr_mse = Rhapsodie.MSE_object(x_est, true_polar_map)
            mse_entry = [k, λ, α, curr_mse[8], curr_mse[9], curr_mse[10]]
            push!(mse_list, mse_entry)
            return sum(curr_mse[8:9])
        end

        optimal_hyperparams, info = PRIMA.newuoa(calculate_MSE_for_prima, [-0.66, -5.], rhobeg=4, rhoend=1e-2, maxfun=100)
        writedlm(io, mse_list, ',')
        close(io)
        println("Optimal hyperparameters lambda, alpha: ", optimal_hyperparams[1], optimal_hyperparams[2])
        println("Info: ", info)
    end