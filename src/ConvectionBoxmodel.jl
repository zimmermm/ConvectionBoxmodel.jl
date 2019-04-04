"""
	ConvectionBoxmodel

API to model a lake overturn as a one-box-model with expanding box size.
Author: Matthias Zimmermann
"""
module ConvectionBoxmodel

# Imports
# ====================================
using DifferentialEquations
using Interpolations: interpolate, Gridded, Linear
using DataFrames
using Calculus: derivative
using CSV
using Dates

# Exports
# ====================================

# ConvectionBoxmodel
export	days_to_datetime,
		datetime_to_days

export	interp1d,
		interp1d!,
		Interpolation,
		LinearInterpolation,
		top,
		bottom

# bathymetries
export	surface_area,
		volume_above,
		total_volume,
		area_at

# growthmodels
export	MOBGrowthModel,
		Nogrowth,
		MOBMonodModel,
		mox,
		μ,
		scale,
		montecarlo_shuffle

# profileshapes
export	FermiProfile,
		GompertzProfile,
		StepProfile,
		DataProfile

# callbacks
export	fully_mixed,
		mixed_to

# convectionlakephysics
export	LakePhysics,
		DefaultPhysics,
		ConvectionLakePhysics,
		ConvectionLakePhysicsScenario,
		MeteorologicalForcing,
		forcing_from_dataset,
		Inflow,
		inflow_from_dataset,
		@physicsfn,
		ρ,
		N2,
		dTdt,
		dhdt,
		ϵ, ϵ_u, ϵ_B,
		k, k_u, k_B,
		Fatm, Fatm_u, Fatm_B,
		F_diff,
		AML,
		buoyancy_flux

# core
export	LakeModel,
		boxmodel_ode,
		createODEProblem,
		solve_boxmodel,
		solve_boxmodel_montecarlo,
		boxmodel_prob_func

#========================
Date conversion
=========================#

function datetime_to_days(datetime::Dates.DateTime)
	# Converts datetime object to days since 01.01.1970 as Float64.
	(datetime-Dates.DateTime(1970)).value/1.0e3/3600.0/24.0
end
precompile(datetime_to_days, (Dates.DateTime,))

function days_to_datetime(days)
	#Converts days since 01.01.1970 to DateTime object.
	Dates.DateTime(1970)+Dates.Nanosecond(floor(days*24.0*3600.0*1.0e9))
end
precompile(days_to_datetime, (Float64,))
precompile(days_to_datetime, (Int64,))

#========================
Linear interpolation
=========================#

# generic interpolator datastructure
struct Interpolation{InterpolationMethod}
	x::Array{<:Real,1}  # 1d domain
	y::Array{<:Real,1}  # 1d data
	at::Function  # interpolator
end

top(i::Interpolation{<:Any}) = i.x[1]
bottom(i::Interpolation{<:Any}) = i.x[end]

# linear interpolator
interp1d(x::Array{T,1}, y::Array{T,1}) where T<:Number = begin
							itp = interpolate((x,), y, Gridded(Linear()))
							(at) -> itp(at)
						end
precompile(interp1d, (Array{Float64,1}, Array{Float64,1}))
interp1d!(x,y) = interp1d(convert(Array{Float64}, x), convert(Array{Float64}, y)) 

# implementation for linear interpolation
const LinearInterpolation = Interpolation{:Linear}
(::Type{Interpolation{:Linear}})(x::Array{<:Real,1}, y::Array{<:Real,1}) = Interpolation{:Linear}(x, y, interp1d(x, y))
(::Type{Interpolation{:Linear}})(x::Array{<:Float64,1}, y::Array{<:Int64,1}) = Interpolation{:Linear}(x, y, interp1d!(x, y))


# Includes
# =========================
include("profileshapes.jl")
include("callbacks.jl")
include("bathymetries.jl")
include("convectionlakephysics.jl")
include("growthmodels.jl")
include("core.jl")

end