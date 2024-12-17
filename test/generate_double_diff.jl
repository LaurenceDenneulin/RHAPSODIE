using Rhapsodie
using DelimitedFiles
using EasyFITS
using InterpolationKernels

#include("test_data_simulation.jl")
#-----------------------------------------------
# Loading the dataset parameters
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

parameter_type = "intensities"

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end

psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

Rhapsodie.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot
, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)


    load_data("test_results/contrast_10e-2.0/DATA.fits", 
            "test_results/contrast_10e-2.0/WEIGHT.fits")
	Sdim=length(Rhapsodie.dataset)
    
    DATA=zeros(get_par().cols[1], get_par().cols[2], Sdim,2);
    WEIGHT=zeros(get_par().cols[1], get_par().cols[2], Sdim,2);

    ker= CatmullRomSpline(Float64, Flat)       
    input_size=(get_par().rows[1], get_par().rows[2]÷2);
    output_size= get_par().cols[1:2];

# Pre processing
    for i=1:Sdim
        T1=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(Rhapsodie.Star_Disk_Table[i][2]))
        T2=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(Rhapsodie.Star_Disk_Table[i][4]))
    
        I1=T1*Rhapsodie.dataset[i].data[:,1:end÷2]
        I2=T2*Rhapsodie.dataset[i].data[:,end÷2+1:end]
        W1=T1*Rhapsodie.dataset[i].weights[:,1:end÷2]
        W2=T2*Rhapsodie.dataset[i].weights[:,end÷2+1:end]

        DATA[:,:,i,1]=I1;
        DATA[:,:,i,2]=I2;
        WEIGHT[:,:,i,1]=W1;
        WEIGHT[:,:,i,2]=W2;
    end
#-----------------------------------------------------
# Double Difference   
DD=Double_Difference(DATA);
write_polar_map(DD, "test_results/contrast_10e-2.0/Results_Separable_DoubleDifference.fits", overwrite=true)