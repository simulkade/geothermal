# geothermal energy
Files related to our geothermal energy analysis paper and talk, presented in ECOS 2019 conference in Wroklaw, Poland. A preliminary version of the work was also presented in the Energy Transition symposium in Delft, the Netherlands.

## How to use
First, download and install Julia for your platform from [here](https://julialang.org/downloads/)  
Then, install the following packages by opening Julia prompt and typing
```
]add package_name
```
where the package_name is `IJulia, DataFrames, PyPlot, Polynomials, GR`.  
You also need to install the latest versions of [CoolProp](https://github.com/CoolProp/CoolProp.jl) and [JFVM](https://github.com/simulkade/JFVM.jl) as instructed in the links.  

From here, you can go through the files with `.ipynb` extension and run them to get your own results.  

## Only the published data
If you want to look at the published data and its analysis, follow the above steps and finally open the file `data_analysis_results.ipynb` and run it. You can also see the analysis by clicking [here](https://nbviewer.jupyter.org/github/simulkade/geothermal/blob/master/data_analysis_results.ipynb)  

## Citation
The work is currently being cleaned up to be submitted to a journal. For the time being, you can cite it as
```
Ali Akbar Eftekhari, Hamidreza M. Nick, Hedzer van der Kooi, and Hans Bruining, Thermodynamic analysis of low-temperature fossil-fuel-dependent geothermal energy in the Netherlands, PROCEEDINGS OF ECOS 2019 - THE 32ND INTERNATIONAL CONFERENCE ON
EFFICIENCY, COST, OPTIMIZATION, SIMULATION AND ENVIRONMENTAL IMPACT OF ENERGY SYSTEMS
JUNE 23-28, 2019, WROCLAW, POLAND
```
