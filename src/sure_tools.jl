#
# sure_tools.jl
#
# Provides an implementation of the Generalized Stein 
# Umbiased Risk Estimator (GSURE) and of the prevision error.
#
# [Stein, 1981] Stein, C. M. (1981). Estimation of the mean of a  
# multivariate normal distribution. The annals of Statistics, 
# pages 1135–1151.
#
# [Eldar, 2008] Eldar, Y. C. (2008). Generalized sure for
# exponential families : Applications to regularization. 
# IEEE Transactions on Signal Processing, 57(2) :471–481.
#
# [Ramani et al., 2008] Ramani, S., Blu, T., and Unser, M. (2008). 
# Monte-Carlo Sure : A Black-Box Optimization of Regularization 
# Parameters for General Denoising Algorithms. IEEE Trans.
# on Image Process., 17(9) :1540–1554.
#
#
#------------------------------------------------------
#
# This file is part of Rhapsodie
#
# Copyright (c) 2017-2021 Laurence DENNEULIN (see Licence.md)
#

#------------------------------------------------


function MSE_data(x_est::Array{T,N}, x_true::Array{T,N}, d::Array{data_table,1}) where {T <: AbstractFloat,N}
    MSE=0.0;
    n=0;
    for  data in d
        res = MSE_data(x_est, x_true, data)
        MSE +=res[1];
        n += res[2];
    end
    return MSE,n
end

function MSE_data(x_est::Array{T,N}, x_true::Array{T,N}, d::data_table) where {T <: AbstractFloat,N}
        minimum(d.weights) >=0 || error("invalid weights");    

        res= d.H*(x_est - x_true);
        MSE = vdot( res, d.weights.*res);
        n= count(d.weights .>0);
    return MSE,n
end

function sure_crit(x::Array{T,N},
                   δx::Array{T,N}, 
                   d::Array{data_table,1}, 
                   δd::Array{data_table,1}; 
                   mask::Array{Array{Bool,2},1}=Array{Array{Bool,2},1}()) where {T <: AbstractFloat,N}

    data_fidelity=0.0;
    trace= 0.0;
    n=0;

    for k=1:length(d)
        if !isempty(mask)       
            res = sure_crit(x, δx, d[k], δd[k], mask=mask[k])
        else
            res = sure_crit(x, δx, d[k], δd[k])   
        end
        data_fidelity += res[1];
        trace += res[2];
        n += res[3];        
    end
    return [data_fidelity ; trace; n]

end

function sure_crit(x::Array{T,N},
                   δx::Array{T,N}, 
                   d::data_table, 
                   δd::data_table; 
                   mask::Array{Bool,2}=Array{Bool,2}(undef, 0,0)) where {T <: AbstractFloat,N}

        MAP = (d.weights .!=0); #Valid pixels map
        
        !isempty(mask) && (size(MAP)==size(mask)) && (MAP .*= mask);
                
        n=count(MAP);
        if n != 0;
            εb=(δd.data -d.data)[MAP]; #Perturbation of variance ε^2
            εd=d.data - d.H*x; #difference between data and model
            εx=(d.H*(δx - x))[MAP]; #difference between model estimated on perturbed data and on data

            data_fidelity = vdot(εd, d.weights.*εd);  
            
            trace = n*vdot(εb, εx)/vdot(εb, εb)
        else
            data_fidelity=0.
            trace=0.
        end      
    return [data_fidelity ; trace; n]

end

function SURE(solver::Function, 
              par::Array{Float64,1}, 
              d::Array{data_table,1}, 
              δd::Array{data_table,1},       
              x0::Array{Float64,3}; 
              mask::Array{Array{Bool,2},1}=Array{Array{Bool,2},1}())

    x  = solver(x0, d, par);
    δx = solver(x0, δd, par);
 
    return sure_crit(x,δx, d, δd, mask=mask), x
end

function do_data_perturbation(dset::Array{data_table,1},eps::Float64)
    p = data_table[];
    for d in dset
        push!(p,do_data_perturbation(d,eps));
    end
    return p
end                      

function do_data_perturbation(d::data_table, eps::Float64)
    return data_table(d.data + eps* randn(size(d.data)), d.weights, d.H);
end                      


function sure_optim(solver, x)#; linear=true)
        res=sure_tools.SURE(solver, x, grad_tools.dataset, data_perturbed, X0nl)[1]
    return res[1] + 2*res[2];
end
