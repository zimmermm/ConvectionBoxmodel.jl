__precompile__()
module ConvectionBoxmodel

#=========================================
This module contains the API to
model a mixing lake as a one-box-model
with expanding box size
=========================================#

#using DifferentialEquations
using OrdinaryDiffEq
using Interpolations: interpolate, Gridded, Linear

export	ContinuousBathymetryObject,
		LakeModelObject,
		solve_boxmodel,
		interp1d,
		datetime_to_days,
		days_to_datetime

function interp1d(x::Array{Float64,1}, y::Array{Float64,1})
	itp = interpolate((x,), y, Gridded(Linear()))
	function interpolated_at(xi)
		itp[xi]
	end
	return interpolated_at
end
precompile(interp1d, (Array{Float64,1}, Array{Float64,1}))

#===================================
Bathymetry
===================================#

# Object describing a continuous bathymetry
struct ContinuousBathymetryObject
	# quick access properties
	surface_area::Float64
	total_volume::Float64
	depth::Float64
	# functions
	area::Function
	volume_above::Function
end

#===================================
Lake Model
===================================#

# Object holding user defined functions
struct LakeModelObject
	bathymetry::ContinuousBathymetryObject

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
	
	biomass_growth_term::Function

	starttime::Float64
	endtime::Float64
	model_callback::ContinuousCallback
end

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

# Definition of the ODE System
function boxmodel_ode(du,u,lakemodel,t)
	h_mix, V, C, B, T, Qatm, Qmix, H = u

	# include morphology to support non-equidistant time steps
	dhdt = lakemodel.dhdt(u, t)
	dVdt = dhdt*lakemodel.bathymetry.area(h_mix)

	# Flux to the atmosphere
	dQatm = lakemodel.Fatm(u,t)*lakemodel.bathymetry.surface_area
	# Flux from hypolimnion into mixed layer
	dQmix = dVdt*lakemodel.concentration_profile(h_mix)

	dHdt = 1/V*lakemodel.dHdt(u,t)*lakemodel.bathymetry.surface_area

	# basic ODE for expanding box size
	dCdt = 1/V*(dQmix-dVdt*C)+1/V*dQatm
	dBdt = -1/V*dVdt*B
	dTdt = 1/V*dVdt*(lakemodel.lake_temperature(h_mix)-T)+1/V*lakemodel.dHdt(u,t)*lakemodel.bathymetry.surface_area

	mu = lakemodel.biomass_growth_term(u, t)
	du[:] = [dhdt, dVdt, dCdt, dBdt, dTdt, dQatm, dQmix, dHdt]+mu*B
end
precompile(boxmodel_ode, (Array{Float64,1}, Array{Float64,1}, LakeModelObject, Float64))


# solver
function solve_boxmodel(lakemodel)
	# assemble initial conditions
	h_mix0, C0, B0, T0 = lakemodel.initial_condition
	V0 = lakemodel.bathymetry.volume_above(h_mix0)
	# diagnostic variables
	Qatm0 = 0.
	Qmix0 = 0.
	H0 = 0.
	u0 = [h_mix0, V0, C0, B0, T0, Qatm0, Qmix0, H0]

	# tspan
	tspan = (lakemodel.starttime,lakemodel.endtime)
	
	# ODEProblem
	prob = ODEProblem(boxmodel_ode,u0,tspan,lakemodel,callback=lakemodel.model_callback)

	# solve
	@time solve(prob, Rodas5(), reltol=1e-5,abstol=1e-5)
end
precompile(solve_boxmodel, (LakeModelObject,))

include("ConvectionBoxmodelToolbox.jl")

end # module