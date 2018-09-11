__precompile__()
module ConvectionBoxmodel

#=========================================
This module contains the API to
model a mixing lake as a one-box-model
with expanding box size
=========================================#

#using DifferentialEquations
using DifferentialEquations
using Interpolations: interpolate, Gridded, Linear
using DataFrames
using Calculus: derivative
using CSV
using Dates

export	@load_prototype,
		call_prototype

export	interp1d,
		interp1d!,
		Interpolation,
		LinearInterpolation,
		top,
		bottom

export	surface_area,
		volume_above,
		total_volume

export	MOBGrowthModel,
		Nogrowth,
		MOBMonodModel,
		mox,
		μ


####################################
# Prototyping
####################################
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

####################################
# Linear Interpolation
####################################
interp1d(x::Array{T,1}, y::Array{T,1}) where T<:Number = begin
							itp = interpolate((x,), y, Gridded(Linear()))
							(at) -> itp[at]
						end
precompile(interp1d, (Array{Float64,1}, Array{Float64,1}))
interp1d!(x,y) = interp1d(convert(Array{Float64}, x), convert(Array{Float64}, y)) 


struct Interpolation{InterpolationMethod}
	x::Array{<:Real,1}
	y::Array{<:Real,1}
	at::Function
end

const LinearInterpolation = Interpolation{:Linear}
(::Type{Interpolation{:Linear}})(x::Array{<:Real,1}, y::Array{<:Real,1}) = Interpolation{:Linear}(x, y, interp1d(x, y))
(::Type{Interpolation{:Linear}})(x::Array{<:Float64,1}, y::Array{<:Int64,1}) = Interpolation{:Linear}(x, y, interp1d!(x, y))

top(i::Interpolation{<:Any}) = i.x[1]
bottom(i::Interpolation{<:Any}) = i.x[end]

#function trapz(x::Array{Float64,1}, y::Array{Float64,1})
#	dx = x[2:end]-x[1:end-1]
#	yt = y[1:end-1]
#	[0;cumsum(yt.*dx)]
#end
#precompile(trapz, (Array{Float64,1}, Array{Float64,1}))


#include("core.jl")
include("bathymetries.jl")
#include("lakephysics.jl")
#include("profileshapes.jl")
include("growthmodels.jl")
#include("callbacks.jl")

end # module
