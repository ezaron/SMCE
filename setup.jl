
using Pkg
Pkg.add.([
        "IJulia"
        "NetCDF",
        "Interpolations",
        "ColorSchemes",
        "Printf",
        "CairoMakie",
        "WGLMakie",
        "DelimitedFiles",
        "LinearAlgebra",
        "DSP",
        "FFTW"
])
println("Your Julia environment should now be complete.")
println("Click on File > New Launcher > Notebook > Julia to launch a Julia session.")
println("Alternately, open the .ipynb in this directory.")
run(`mkdir -p ./Figures`)
exit()
