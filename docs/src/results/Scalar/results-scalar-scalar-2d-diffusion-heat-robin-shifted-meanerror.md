# Scalar 2D Diffusion Heat Robin Shifted MeanError

Scalar 2D transient diffusion with Robin interface conditions. The circular
interface of radius `R` is translated within the square domain and the global
error at final time is averaged over the translation space
`(x₀, y₀) ∈ [-L/2, L/2]²`, truncated so the circle remains inside the box.
The integral `∫∫ e(x₀, y₀) dx₀ dy₀` is approximated by the uniform mean of the
sampled errors.

**CSV source:** `results/scalar/Scalar_2D_Diffusion_Heat_Robin_Shifted_MeanError.csv`

| h | mean_all_err | mean_full_err | mean_cut_err | mean_empty_err | samples |
| --- | --- | --- | --- | --- | --- |
| 1 | 0.196114 | 0.0801228 | 0.178276 | 0 | 25 |
| 0.5 | 0.110816 | 0.0694031 | 0.0858086 | 0 | 25 |
| 0.25 | 0.0345676 | 0.0276349 | 0.02071 | 0 | 25 |
| 0.125 | 0.00787116 | 0.00683871 | 0.00389232 | 0 | 25 |
| 0.0625 | 0.00225706 | 0.00206123 | 0.000918099 | 0 | 25 |
