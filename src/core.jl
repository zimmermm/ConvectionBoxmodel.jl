#============================================
* Data structures
* System of ODE
* Solver
=============================================#

#===================================
Bathymetry
===================================#

# Interface describing a continuous bathymetry
struct ContinuousBathymetryInterface
	# quick access properties
	surface_area::Float64
	total_volume::Float64
	depth::Float64
	renewal_interface_volume::Float64
	renewal_interface_area::Float64
	# functions
	area::Function
	volume_above::Function
end

#===================================
Lake Model
===================================#

# Interface to user defined functions of the model
struct LakeModelInterface
	bathymetry::ContinuousBathymetryInterface

	initial_condition::Array{Float64,1}

	# static profiles
	lake_temperature::Function
	concentration_profile::Function
	
	# thickening rate of mixed layer
	dhdt::Function
	# heat flux
	dHdt::Function
	# atmospheric exchange
	Fatm::Function
	# Diffusive flux from hypolimnion
	F_diff::Function
	
	biomass_growth_term::Function

	starttime::Float64
	endtime::Float64
	model_callback::ContinuousCallback
end


######################
# date conversion
######################
function datetime_to_days(datetime::DateTime)
	(datetime-DateTime(1970)).value/1e3/3600/24.
end
precompile(datetime_to_days, (DateTime,))

function days_to_datetime(days)
	DateTime(1970)+Dates.Day(days)
end
precompile(days_to_datetime, (Float64,))
precompile(days_to_datetime, (Int64,))

#===================================
Boxmodel
===================================#

function boxmodel_ode(du,u,lakemodel,t)
	#h_mix, V, C, B, T, Qatm, Qmix, H, Qmox = u
	h_mix, V, C, B, T, Qatm, Qmix, Qmox, Qdiff = u

	# include morphology to support non-equidistant time steps
	du[1] = lakemodel.dhdt(u, t)*3600*24
	du[2] = du[1]*lakemodel.bathymetry.area(h_mix)

	# Flux to the atmosphere
	du[6] = -lakemodel.Fatm(u,t)*lakemodel.bathymetry.surface_area*3600*24
	# Flux from hypolimnion into mixed layer
	du[7] = du[2]*lakemodel.concentration_profile(h_mix)+lakemodel.F_diff(u,t)*lakemodel.bathymetry.area(h_mix)*3600*24
	# MOX
	du[8] = 0
	du[9] = lakemodel.F_diff(u,t)*lakemodel.bathymetry.area(h_mix)*3600*24

	# basic ODE for expanding box size
	du[3] = 1/V*(du[7]-du[2]*C-du[6])
	du[4] = -1/V*du[2]*B
	du[5] = 1/V*du[2]*(lakemodel.lake_temperature(h_mix)-T)+1/V*lakemodel.dHdt(u,t)*lakemodel.bathymetry.surface_area*3600*24

	mu = lakemodel.biomass_growth_term(u, t)*3600*24*B
	du[:] = du+mu
	du[8] = du[8]*V
	return
end
precompile(boxmodel_ode, (Array{Float64,1}, Array{Float64,1}, LakeModelInterface, Float64))

# solver
function solve_boxmodel(lakemodel)
	# assemble initial conditions
	h_mix0, C0, B0, T0 = lakemodel.initial_condition
	V0 = lakemodel.bathymetry.volume_above(h_mix0)
	# diagnostic variables
	Qatm0 = 0.
	Qmix0 = 0.
	Qmox0 = 0.
	H0 = 0.
	Qdiff0 = 0.
	u0 = [h_mix0, V0, C0, B0, T0, Qatm0, Qmix0, Qatm0, Qdiff0]

	# tspan
	tspan = (lakemodel.starttime,lakemodel.endtime)

	# ODEProblem
	prob = ODEProblem(boxmodel_ode,u0,tspan,lakemodel,callback=lakemodel.model_callback)

	# solve
	@time solve(prob, Rosenbrock23(), reltol=1e-3, abstol=1e-3, dtmax=1./24., force_dtmin=true)
end
precompile(solve_boxmodel, (LakeModelInterface,))
