#===================================
Forcing Data
===================================#

struct MeteorologicalForcing
	air_temperature::Interpolation
	cloud_cover::Interpolation
	vapour_pressure::Interpolation
	global_radiation::Interpolation
	wind_speed::Interpolation
end

function forcing_from_dataset(path)
	# load simstrat forcing dataset
	forcing_df = CSV.read(path, delim="\t")

	##########################################
	# interpolators for individual data series
	##########################################
	timestamps = forcing_df[:t]
	air_temperature = LinearInterpolation(timestamps, forcing_df[Symbol("Tair (°C)")]+273.15)
	cloud_cover = LinearInterpolation(timestamps, forcing_df[Symbol("cloud coverage")])
	vapour_pressure = LinearInterpolation(timestamps, forcing_df[Symbol("vap (mbar)")])

	# wind speed
	u10 = forcing_df[Symbol("u (m/s)")]
	v10 = forcing_df[Symbol("v (m/s)")]
	wind_speed = LinearInterpolation(timestamps, sqrt.(u10.^2+v10.^2))

	# global shortwave irradiation
	G = forcing_df[Symbol("Fsol (W/m2)")]
	C = forcing_df[Symbol("cloud coverage")]

	## Wuest et al.
	#Fdir = (1.-C)./((1.-C)+0.5*C)
	#Fdiff = (0.5*C)./((1.-C)+0.5*C)
	#Alb_dir = 0.2
	#Alb_diff = 0.066
	#Hs = G.*Fdir*(1-Alb_dir)+G.*Fdiff*(1-Alb_diff)

	## Simstrat
	Hs = (1-0.08)*G
	global_radiation = LinearInterpolation(timestamps, Hs)

	MeteorologicalForcing(air_temperature, cloud_cover, vapour_pressure, global_radiation, wind_speed)
end
precompile(forcing_from_dataset, (String,))

#===================================
Inflow
===================================#
struct Inflow
	flow::Interpolation
	temperature::Interpolation
end

function inflow_from_dataset(path)
	# load simstrat forcing dataset
	inflow_df = CSV.read(path, delim=",")

	##########################################
	# interpolators for individual data series
	##########################################
	timestamps = inflow_df[:Timestamp]
	flow = LinearInterpolation(timestamps, inflow_df[:Flow])
	temperature = LinearInterpolation(timestamps, inflow_df[:Temp]+273.15)

	Inflow(flow, temperature)
end
precompile(inflow_from_dataset, (String,))

#===================================
LakePhysics
===================================#
# traits
abstract type DefaultPhysics end

# lake physics
abstract type LakePhysics end

struct ConvectionLakePhysics{Traits} <: LakePhysics
	# constants
	a_T::Float64
	g::Float64
	β::Float64
	Mw::Float64
	Cp::Float64
	Cp_air::Float64
	Lv::Float64
	σ::Float64
	κ::Float64
	ν::Float64
	ν_a::Float64

	T0::Float64
	p0::Float64
	R::Float64

	rho::Float64
	rho_air::Float64
	p_air::Float64

	γ::Float64

	# empirical constants
	A_L::Float64
	a::Float64
	η::Float64
	B::Float64
	c1::Float64

	f_wind::Float64

	temperature_profile::Interpolation
	salinity_profile::Interpolation
	concentration_profile::Interpolation
	forcing::MeteorologicalForcing
	inflow::Inflow
	bathymetry::Interpolation
	A::Float64
	cheat1::Float64
	cheat2::Float64
	cheat3::Float64
end

# convection lake physics constructor
ConvectionLakePhysics{Traits}(temperature_profile, salinity_profile, concentration_profile, forcing, inflow, bathymetry, A, cheat1, cheat2, cheat3) where Traits <: Any = begin
		a_T = 2.14e-4 # thermal expansion [K-1]
		g = 9.80665  # gravitational acceleration [m^2 s-1]
		β = g*a_T
		Mw = 18.01528/1000  # molar mass of water [kg mol-1]
		Cp = 75.375/Mw  # specific heat capacity [J kg-1 K-1]
		Cp_air = 1005  # specific heat capacity of air [J kg-1 K-1]
		Lv = 2.47e6  # heat vaporization of water [J kg-1]
		σ = 5.670367e-8  # Stefan Boltzmann constant [W m-2 K-4]
		κ = 0.4  # Von Karman constant
		ν = 1.5e-6  # Kinematic Viscosity of Water
		ν_a = 1.4e-5  # Kinematic Viscosity of Air [m2 s-1]

		T0 = 273.15  # Standard Temperature [K]
		p0 = 1e5	# Standard Pressure [Pa]
		R = 8.314  # universal gas constant

		rho = 998  # density of water [kg m-3]
		rho_air = 1.25  # density of air [kg m-3]
		p_air = 960  # air pressure [mbar]

		γ = (Cp_air*p_air)/(0.622*Lv)

		# empirical constants
		A_L = 0.03  # Longwave Alebedo [-]
		a = 1.09  # calibration constant
		η = 1/(2*pi)  # according to Lorke et al. #0.29 # calibration constant
		B = 0.62  # Bowen coefficient [mbar K-1]
		c1 = 8.6

		f_wind = 0.5

		ConvectionLakePhysics{Traits}(
			a_T,
			g,
			β,
			Mw,
			Cp,
			Cp_air,
			Lv,
			σ,
			κ,
			ν,
			ν_a,

			T0,
			p0,
			R,

			rho,
			rho_air,
			p_air,

			γ,

			# empirical constants
			A_L,
			a,
			η,
			B,
			c1,

			f_wind,

			temperature_profile,
			salinity_profile,
			concentration_profile,
			forcing,
			inflow,
			bathymetry,
			A,
			cheat1,
			cheat2,
			cheat3
			)
	end
ConvectionLakePhysics(temperature_profile, salinity_profile, concentration_profile, forcing, inflow, bathymetry, A, cheat1, cheat2, cheat3) = ConvectionLakePhysics{DefaultLakePhysics}(temperature_profile, salinity_profile, concentration_profile, forcing, inflow, bathymetry, A, cheat1, cheat2, cheat3)

# macro to define lake physics functions in short notation
macro physicsfn(fnexpr)
	fnrepr = repr(fnexpr)
	fnsignature = split(fnrepr[3:end-1], "=")[1]
	# first argument of function signature
	fnfirstarg = split(split(match(r"\((.*?)\)", fnsignature)[1], ",")[1],"::")
	fn_classname = fnfirstarg[1]
	if length(fnfirstarg) > 1
		# look for field names and prefix them with the class
		fn_class_fieldnames = fieldnames(eval(Meta.parse(fnfirstarg[2])))
		for field in fn_class_fieldnames
			fnrepr = replace(fnrepr, Regex(string("(?<![\\w_])(|\\d+)(", field, ")(?![\\w_])"))=>SubstitutionString(string("\\g<1>",fn_classname, ".", field)))
		end
	else
		error(string("You have to define the type of ", fn_classname ," in the function signature"))
	end
	eval(Meta.parse(fnrepr))
end

########################
# lake physics functions
########################
# toplevel functions used by the solver
# dhdt(u,t)
# dTdt(u,t)
# Fatm(u,t)
# F_diff(u,t)

# density
#-----------
const ai = [999.8395, 6.7914e-2, -9.0894e-3, 1.0171e-4, -1.2846e-6, 1.1592e-8, -5.0125e-11]
const bi = [0.8181, -3.85e-3, 4.96e-5]

ρ(T::Float64,S::Float64) =	begin
								Ti = [(T-273.15)^i for i in 0:6]
								dot(ai, Ti)+S*dot(bi, Ti[1:3])
							end

@physicsfn ρ(p::ConvectionLakePhysics{<:DefaultPhysics}, depth) = ρ(temperature_profile.at(depth), salinity_profile.at(depth))

# water column stability
#--------------------------
@physicsfn N2(p::ConvectionLakePhysics{<:DefaultPhysics}, depth::Float64, Tmix::Float64) =	begin
																								ρmix = density(Tmix, 0.1)  # density in mixed layer. assumption that salinity in mixed layer is constantly 0.1 PSU
																								ρhypo = density(p, depth)  # density of hypolimnion at mixed layer depth
																								ρav = (ρmix+ρhypo)/2  # average density
																								dρdz = (ρhypo-ρmix)/0.1  # density gradient. assumption that gradient is 0.1m thick
																								g/ρhypo*dρdz  # Brunt-Väisälä frequency
																							end

# Heat flux
#----------------

# Heat flux at the surface
##########################
# positive flux: warming
# negative flux: cooling
#*** Longwave
# emissivity
Ea(Ta, cloud_cover, vapour_pressure) = (1.0+0.17*cloud_cover^2)*1.24*(vapour_pressure/Ta)^(1.0/7.0)
# atmospheric longwave radiation
Ha_simstrat(Ta, cloud_cover, vapour_pressure, σ) = (1.0-A_L)*Ea(Ta, cloud_cover, vapour_pressure)*σ*Ta^4
# water longwave radiation
Hw_simstrat(Tw, σ) = -0.972*σ*Tw^4
#*** Shortwave
# is given by forcing.global_radiation(t) (already corrected for albedo)
#*** evaporation & condensation
fu_simstrat(U10, Tw, Ta) = 4.4+1.82*U10^2+0.26*(Tw-Ta)
e_s_simstrat(Tw, Ta)=6.107*10^(7.5*(Tw-273.15)/(237.3+(Tw-273.15)))
He_simstrat(U10, Tw, Ta, vapour_pressure) = -fu_simstrat(U10, Tw, Ta)*(e_s_simstrat(Tw, Ta)-vapour_pressure)
#*** sensible heat
Hc_simstrat(U10, Tw, Ta) = -B*fu_simstrat(U10, Tw, Ta)*(Tw-Ta)
#*** inflow/outflow
Hfl(Ta,Tw,flow,temp,surface_area) = rho*Cp*flow/surface_area*(temp-Tw)

# total heat balance [W m-2]
heat_flux_simstrat(U10, Tw, Ta, global_radiation, cloud_cover, vapour_pressure) = cheat1*(Hc_simstrat(U10, Tw, Ta) + He_simstrat(U10, Tw, Ta, vapour_pressure))+cheat2*Ha_simstrat(Ta, cloud_cover, vapour_pressure)+Hw_simstrat(Tw)+cheat3*global_radiation

@physicsfn heat_flux(p::ConvectionLakePhysics{<:DefaultPhysics}, u, t) = (heat_flux_simstrat(f_wind*forcing.wind_speed.at(t), u[5], forcing.air_temperature.at(t), forcing.global_radiation.at(t), forcing.cloud_cover.at(t), forcing.vapour_pressure.at(t)))#+Hfl(forcing.air_temperature(t), u[5], inflow.flow(t), inflow.temperature(t))
