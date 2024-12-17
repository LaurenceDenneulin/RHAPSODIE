using Rhapsodie
using EasyFITS

true_polar_map = Rhapsodie.read_and_fill_polar_map(parameter_type, "test_results/contrast_10e-2.0/TRUE.fits")

x = readfits("test_results/contrast_10e-2.0/.fits")

Rhapsodie.MSE_object(x, true_polar_map)

println("MSE: ", mse)
