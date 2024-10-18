#
# Rhapsodie.jl
#
# Package for the Reconstruction of High-contrAst Polarized
# SOurces and Deconvolution for cIrcumstellar Environments (Rhapsodie)
#
#----------------------------------------------------------
#
# 
# Copiyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

module Rhapsodie

    export
        PolarimetricPixel,
        PolarimetricMap,
        TPolarimetricPixel,
        TPolarimetricMap,
        write_polar_map,
        read_and_fill_polar_map,
        convert,
        load_parameters,
        get_par,
        load_data,
        fg!,
        SetCropOperator,
        crop,
        crop!,
        pad,
        set_fft_op,
        TwoDimensionalTransformInterpolator,
        FieldTransformOperator,
        TFieldTransformOperator,
        data_simulator,
        generate_model,
        data_generator,
        generate_parameters,
        Double_Difference,
        Double_Ratio,
        Linear_Method,
        NonLinear_Method,
        apply_rhapsodie,
        apply_edge_preserving_smoothing!,
        MSE_object

    import Base: +, -, *, /, ==, getindex, setindex!, read, write, convert

    using OptimPackNextGen
    import OptimPackNextGen: BraDi #va devenir BraDi avec un D majuscule
    using SpecialFunctions
    using TwoDimensional
    using FFTW
    using LinearInterpolators
    using Statistics
    using LinearAlgebra
    using LazyAlgebra
    import LazyAlgebra: Mapping, vcreate, vcopy, apply!
    using StaticArrays
    using FITSIO
    using EasyFITS
    using DelimitedFiles
    #using TiPi
    using Random

    include("Polarimetric_Parameters.jl")
    include("tpolarimetric_parameters.jl")
    include("grad_tools.jl")
    include("tgrad_tools.jl")
    include("loaders.jl")
    include("separable_methods.jl")
    include("sure_tools.jl")
    include("rhapsodie_methods.jl")
    include("trhapsodie_methods.jl")
    include("datasimul_tools.jl")
    include("tdatasimul_tools.jl")
    include("utils.jl")
end