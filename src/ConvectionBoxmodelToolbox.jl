__precompile__()
module ConvectionBoxmodelToolbox
using ConvectionBoxmodel
using OrdinaryDiffEq
#using ParameterizedFunctions: 
using Interpolations: interpolate, Gridded, Linear
using Calculus: derivative
using DataFrames
using CSV

export	ContinuousBathymetry,
		ConvectionLakePhysics,
		ConstantRate,
		FirstOrderRate,
		Disable,
		FermiProfile,
		StepProfile,
		GompertzProfile,
		MOBGrowthModel,
		NoBiomass,
		fully_mixed,
		mixed_to,
		forcing_from_dataset,
		initial_T_profile_from_simstrat_ic,
		initial_T_profile_from_simstrat_output

#===================================
Bathymetry
===================================#

# Constructor for a linearly interpolated bathymetry
function ContinuousBathymetry(depths::Array{Float64,1}, areas::Array{Float64,1})
	# function definitions
	area = interp1d(depths, areas)

	function equigrid(x0,xend,dh=0.01)
		n = round((xend-x0)/dh)
		dh = (xend-x0)/n
		collect(x0:dh:xend)
	end

	function volume_above(depth, dh=0.01)
		grid = equigrid(0.,depth,dh)
		areas_i = area(grid)
		sum((areas_i[1:end-1]+areas_i[2:end])/2.)*dh
	end

	min_depth = minimum(depths)
	max_depth = maximum(depths)

	# create object
	ContinuousBathymetryObject(area(min_depth), volume_above(max_depth), max_depth, area, volume_above)
end

#===================================
Forcing Data
===================================#

struct ForcingObject
	air_temperature::Function
	cloud_cover::Function
	vapour_pressure::Function
	global_radiation::Function
	wind_speed::Function
end

function forcing_from_dataset(path)
	# load simstrat forcing dataset
	forcing_df = CSV.read(path, delim="\t")

	##########################################
	# interpolators for individual data series
	##########################################
	air_temperature = interp1d(forcing_df[:t], forcing_df[Symbol("Tair (°C)")]+273.15)
	cloud_cover = interp1d(forcing_df[:t], forcing_df[Symbol("cloud coverage")])
	vapour_pressure = interp1d(forcing_df[:t], forcing_df[Symbol("vap (mbar)")])

	# wind speed
	u10 = forcing_df[Symbol("u (m/s)")]
	v10 = forcing_df[Symbol("v (m/s)")]
	wind_speed = interp1d(forcing_df[:t], sqrt.(u10.^2+v10.^2))

	# global shortwave irradiation
	G = forcing_df[Symbol("Fsol (W/m2)")]
	C = forcing_df[Symbol("cloud coverage")]

	Fdir = (1.-C)./((1.-C)+0.5*C)
	Fdiff = (0.5*C)./((1.-C)+0.5*C)
	Alb_dir = 0.2
	Alb_diff = 0.066
	Hs = G.*Fdir*(1-Alb_dir)+G.*Fdiff*(1-Alb_diff)
	#global_radiation = interp1d(forcing_df[:t], forcing_df[Symbol("Fsol (W/m2)")])
	global_radiation = interp1d(forcing_df[:t], Hs)

	ForcingObject(air_temperature, cloud_cover, vapour_pressure, global_radiation, wind_speed)
end
precompile(forcing_from_dataset, (String,))

function initial_T_profile_from_simstrat_ic(path)
	ic_df = CSV.read(path, delim="\t")
	interp1d(-ic_df[Symbol("depth (m)")], ic_df[Symbol("T (°C)")]+273.15)
end
precompile(initial_T_profile_from_simstrat_ic, (String,))

function initial_T_profile_from_simstrat_output(path)
	data, header = readdlm(path, header=true)
	df = DataFrame(data)
	names!(df, [Symbol("$name") for name in header[1,:]])
	df
end
precompile(initial_T_profile_from_simstrat_output, (String,))


#===================================
Lake Physics
===================================#

struct LakePhysicsObject
	dhdt::Function
	dHdt::Function
	Ba::Function
end


# Constructor for Convective Lake Physics
function ConvectionLakePhysics(temperature_profile, concentration_profile, forcing, A=0.2)
	# some constants
	a_T = 2.14e-4 # thermal expansion [K-1]
	g = 9.80665  # gravitational acceleration [m^2 s-1]
	rho = 998  # density of water [kg m-3]
	β = g*a_T
	Mw = 18.01528/1000  # molar mass of water [kg mol-1]
	Cp = 75.375/Mw  # specific heat capacity [J kg-1 K-1]
	σ = 5.670367e-8  # Stefan Boltzmann constant [W m-2 K-4]

	# empirical constants
	A_L = 0.03  # Longwave Alebedo [-]
	a = 1.09  # calibration constant
	B = 0.61  # Bowen coefficient [-]


	# function definitions
	# Brunt-väisälä frequency
	N2 = depth -> β*derivative(temperature_profile, depth)
	# Heat flux at the surface
	# emissivity
	Ea = (Ta, cloud_cover, vapour_pressure) -> a*(1+0.17*cloud_cover^2)*1.24*(vapour_pressure/Ta)^(1/7)
	# atmospheric longwave radiation
	Ha = (Ta, cloud_cover, vapour_pressure) -> (1-A_L)*Ea(Ta, cloud_cover, vapour_pressure)*σ*Ta^4
	# water longwave radiation
	Hw = (Tw) -> -0.972*σ*Tw^4
	# sensible heat
	#fu = (wind_scal2, Tw, Ta) -> 4.4+1.82*wind_scal2+
	#Hc = 

	Hflux = (t, Tw, Ta, cloud_cover, vapour_pressure) -> Ha(Ta, cloud_cover, vapour_pressure)+Hw(Tw)+forcing.global_radiation(t)
	# heat flux wrapper
	dHdt = (u, t) -> Hflux(t, u[5], forcing.air_temperature(t), forcing.cloud_cover(t), forcing.vapour_pressure(t))/rho/Cp*3600*24
	# Buoyancy flux
	Ba = (t, Tw) -> g*a_T*Hflux(t,Tw,forcing.air_temperature(t), forcing.cloud_cover(t), forcing.vapour_pressure(t))/(Cp*rho)

	

	# Convective thickening of the mixed layer
	dhdt = (u, t) -> (1+2*A)*Ba(t, u[5])/(N2(u[1])*u[1])*3600*24

	LakePhysicsObject(dhdt, dHdt, Ba)
end
precompile(ConvectionLakePhysics, (Function, Function, ForcingObject, Float64))

# Constructor for processes at constant rate
function ConstantRate(r)
	(u, t) -> r
end
precompile(ConstantRate, (Float64,))

function FirstOrderRate(r, idx)
	(u, t) -> r*u[idx]
end
precompile(FirstOrderRate, (Float64,Int64))

# Constructor for processes at constant rate
function Disable()
	(u, t) -> 0.
end




#===================================
Profile shapes library
===================================#

function FermiProfile(Qtop, Qbottom, z_interface, w_interface)
	depth -> Qbottom+(Qtop-Qbottom)./(1.+exp.((depth-z_interface)/w_interface))
end

function GompertzProfile(Qtop, Qbottom, z_interface, w_interface)
	depth -> Qtop-(Qtop-Qbottom).*exp.(-20/w_interface.*exp.(-w_interface.*(depth-z_interface))./depth)
end

function StepProfile(Qtop, Qbottom, z_bottom, z_interface, sharpness)
	slope = (Qbottom-Qtop)/log(1+exp(sharpness*(z_bottom-z_interface)))
	depth -> Qtop+slope.*log.(1+exp.(sharpness*(depth-z_interface)))
end


#===================================
Growth Model
===================================#

# constructor for growth model
function MOBGrowthModel(vmax, Km, y)
	function growth_term(u, t)
		h_mix, V, C, B, T = u
		mu=vmax*(C/(Km+C))
		[0, 0, -mu, y*mu, 0, 0, 0, 0]
	end
end
precompile(MOBGrowthModel, (Float64, Float64, Float64))

function NoBiomass()
	function growth_rate(u, t)
		[0, 0, 0, 0, 0, 0, 0, 0]
	end
end


# event definition
function fully_mixed()
	condition(u,t,integrator) = 16.-u[1]
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=true)
end

function mixed_to(termination_depth)
	condition(u,t,integrator) = termination_depth-u[1]
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=true)
end

end
