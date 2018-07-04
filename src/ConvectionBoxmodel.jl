__precompile__()
module ConvectionBoxmodel

#=========================================
This module contains the API to
model a mixing lake as a one-box-model
with expanding box size
=========================================#

#using DifferentialEquations
using DifferentialEquations
using ODEInterfaceDiffEq
using Interpolations: interpolate, Gridded, Linear
using DataFrames
using Calculus: derivative
using DataFrames
using CSV

export	ContinuousBathymetryInterface,
		LakeModelInterface,
		solve_boxmodel,
		interp1d,
		interp1d!,
		trapz,
		datetime_to_days,
		days_to_datetime

export	ContinuousBathymetry,
		LakePhysicsInterface,
		ConvectionLakePhysicsInterface,
		ConstantRate,
		FirstOrderRate,
		Disable,
		FermiProfile,
		StepProfile,
		GompertzProfile,
		DataProfile,
		MOBGrowthModel,
		MOBGrowthModelInterface,
		NoBiomass,
		fully_mixed,
		mixed_to,
		forcing_from_dataset,
		inflow_from_dataset,
		initial_T_profile_from_simstrat_ic,
		initial_T_profile_from_simstrat_output

export @load_prototype, call_prototype

macro load_prototype(path)
	quote
		file = open($(esc(path)))
		s = readstring(file)
		close(file)
		eval(parse(s))
	end
end

function call_prototype(f, args...)
	Base.invokelatest(f, args...)
end

# generates a function f(x)->y that interpolates over the 1d data (x,y) provided as arguments
function interp1d(x::Array{Float64,1}, y::Array{Float64,1})
	itp = interpolate((x,), y, Gridded(Linear()))
	function interpolated_at(xi)
		itp[xi]
	end
	return interpolated_at
end
precompile(interp1d, (Array{Float64,1}, Array{Float64,1}))

function trapz(x::Array{Float64,1}, y::Array{Float64,1})
	dx = x[2:end]-x[1:end-1]
	yt = y[1:end-1]#(y[1:end-1]+y[2:end])/2.
	[0;cumsum(yt.*dx)]
	#itp = interp1d(x,y)
	#f=cumsum(Fun(itp, x[1]..x[end]))
	#[f(xi) for xi in x]
end
precompile(trapz, (Array{Float64,1}, Array{Float64,1}))

function interp1d!(x,y)
	return interp1d(convert(Array{Float64}, x), convert(Array{Float64}, y))
end

include("core.jl")
include("bathymetries.jl")
include("lakephysics.jl")
include("profileshapes.jl")
include("growthmodels.jl")
include("callbacks.jl")

end # module
