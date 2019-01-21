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
	forcing_df = disallowmissing!(forcing_df[completecases(forcing_df),:])

	##########################################
	# interpolators for individual data series
	##########################################
	timestamps = forcing_df[:t]
	air_temperature = LinearInterpolation(timestamps, forcing_df[Symbol("Tair (°C)")].+273.15)
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
	inflow_df = disallowmissing!(inflow_df[completecases(inflow_df),:])

	##########################################
	# interpolators for individual data series
	##########################################
	timestamps = inflow_df[:Timestamp]
	flow = LinearInterpolation(timestamps, inflow_df[:Flow])
	temperature = LinearInterpolation(timestamps, inflow_df[:Temp].+273.15)

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
ConvectionLakePhysics{Traits}(temperature_profile, salinity_profile, concentration_profile, forcing, inflow, bathymetry, A, cheat1, cheat2, cheat3, f_wind) where Traits <: Any = begin
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
								ai'*Ti+S*bi'*(Ti[1:3])
							end

@physicsfn ρ(p::ConvectionLakePhysics{<:DefaultPhysics}, depth) = ρ(temperature_profile.at(depth), salinity_profile.at(depth))

# water column stability
#--------------------------
@physicsfn N2(p::ConvectionLakePhysics{<:DefaultPhysics}, depth::Float64, Tmix::Float64) =	begin
																								ρmix = ρ(Tmix, 0.1)  # density in mixed layer. assumption that salinity in mixed layer is constantly 0.1 PSU
																								ρhypo = ρ(p, depth)  # density of hypolimnion at mixed layer depth
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
@physicsfn Ea(p::ConvectionLakePhysics{<:DefaultPhysics}, Ta, cloud_cover, vapour_pressure) = (1.0+0.17*cloud_cover^2)*1.24*(vapour_pressure/Ta)^(1.0/7.0)
# atmospheric longwave radiation
@physicsfn Ha_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, Ta, cloud_cover, vapour_pressure) = (1.0-A_L)*Ea(p,Ta, cloud_cover, vapour_pressure)*σ*Ta^4
# water longwave radiation
@physicsfn Hw_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, Tw) = -0.972*σ*Tw^4
#*** Shortwave
# is given by forcing.global_radiation(t) (already corrected for albedo)
#*** evaporation & condensation
@physicsfn fu_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, U10, Tw, Ta) = 4.4+1.82*U10^2+0.26*(Tw-Ta)
@physicsfn e_s_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, Tw, Ta)=6.107*10^(7.5*(Tw-273.15)/(237.3+(Tw-273.15)))
@physicsfn He_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, U10, Tw, Ta, vapour_pressure) = -fu_simstrat(p, U10, Tw, Ta)*(e_s_simstrat(p, Tw, Ta)-vapour_pressure)
#*** sensible heat
@physicsfn Hc_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, U10, Tw, Ta) = -B*fu_simstrat(p, U10, Tw, Ta)*(Tw-Ta)
#*** inflow/outflow
@physicsfn Hfl(p::ConvectionLakePhysics{<:DefaultPhysics}, Ta,Tw,flow,temp,surface_area) = rho*Cp*flow/surface_area*(temp-Tw)

# total heat balance [W m-2]
@physicsfn heat_flux_simstrat(p::ConvectionLakePhysics{<:DefaultPhysics}, U10, Tw, Ta, global_radiation, cloud_cover, vapour_pressure) = cheat1*(Hc_simstrat(p, U10, Tw, Ta) + He_simstrat(p, U10, Tw, Ta, vapour_pressure))+cheat2*Ha_simstrat(p, Ta, cloud_cover, vapour_pressure)+Hw_simstrat(p, Tw)+cheat3*global_radiation

@physicsfn heat_flux(p::ConvectionLakePhysics{<:DefaultPhysics}, u, t) = (heat_flux_simstrat(p, f_wind*forcing.wind_speed.at(t), u[5], forcing.air_temperature.at(t), forcing.global_radiation.at(t), forcing.cloud_cover.at(t), forcing.vapour_pressure.at(t)))#+Hfl(p, forcing.air_temperature(t), u[5], inflow.flow(t), inflow.temperature(t))

@physicsfn dTdt(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = heat_flux(p, u,t)/(rho*Cp)

# Buoyancy foring: Convective thickening of the mixed layer
##########################################
# Buoyancy Flux [m2 s-3]
@physicsfn buoyancy_flux(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = -β*dTdt(p, u,t)
# thickening rate
# Zilitinkevich 1991
@physicsfn dhdt(p::ConvectionLakePhysics{<:DefaultPhysics}, u, t) = begin
					# no thermocline deepening when the bottom of the lake is reached
					# should be given as parameter in future versions!!
					if u[1] > 15.9
						return 0.0
					end
					B0 = buoyancy_flux(p,u,t)
					ϵ_h = ϵ_u(p, u, t, u[1])+ϵ_B(p,u,t)
					# no thermocline erosion during warming
					# no thermocline erosion if mixed layer is warmer than hypolimnion
					#if ϵ_h<0 | (u[5]>temperature_profile.at(u[1]))
					#	return 0.0
					#else
						#v=(1+2*A)*B0/(N2(p, u[1], u[5])*u[1])
						v=(1+2*A)*ϵ_h/(N2(p, u[1], u[5])*u[1])
						# don't allow rising of the thermocline
						if v < 0.0
							return 0.0
						else
							return v
						end
					#end
				end


# Wind forcing
##########################
# Simstrat 1.6
@physicsfn CD(p::ConvectionLakePhysics{<:DefaultPhysics}, U10) =	begin
																		if U10 <= 0.1
																			0.06215
																		elseif U10 <= 3.85
																			0.0044*U10^(-1.15)
																		else #Polynomial approximation of Charnock's law
																			-0.000000712*U10^2+0.00007387*U10+0.0006605
																		end
																	end
@physicsfn τ(p::ConvectionLakePhysics{<:DefaultPhysics}, U10) = rho_air*CD(p, U10)*U10^2
# wind shear (friction velocity)
@physicsfn u_w(p::ConvectionLakePhysics{<:DefaultPhysics}, U10) = sqrt.(τ(p, U10)/rho)
# thickness of diffusive boundary layer
@physicsfn δv(p::ConvectionLakePhysics{<:DefaultPhysics}, U10) = c1*ν/u_w(p, U10)


# Dissipation
###############

# Read et al.
@physicsfn ϵ_u(p::ConvectionLakePhysics,u,t,z)= 0.6*u_w(p, f_wind*forcing.wind_speed.at(t))^3/(κ*z)
@physicsfn ϵ_u(p::ConvectionLakePhysics,u,t)= ϵ_u(p,u,t,AML(p,u,t))
@physicsfn ϵ_B(p::ConvectionLakePhysics, u,t)=	begin
							B0 = buoyancy_flux(p,u,t)
							if B0 < 0.0
								0.0
							else
								0.5*B0
							end
						end

@physicsfn ϵ(p::ConvectionLakePhysics,u,t) = ϵ_u(p,u,t)+ϵ_B(p,u,t)
#begin
#	return 1e-9
#	ϵ_uB = ϵ_u(p,u,t)+ϵ_B(p,u,t)
#	if ϵ_uB > 1e-7
#		return 1e-7
#	else
#		return ϵ_uB
#	end
#end

# Air/Water transfer velocity
#############################

# Schmitt-Number
Sc_ch4_Wanninkhof(T) = 1909.4 - 120.78*(T-273.15) + 4.1555*((T-273.15)^2) - 0.080578*((T-273.15)^3) + 0.00065777*((T-273.15)^4)
@physicsfn Ceq_ch4(p::ConvectionLakePhysics{<:DefaultPhysics},T,S) = exp(-68.8862 + 101.4956*(100/T) + 28.7314*log(T/100) + S*(-0.076146 + 0.043970*(T/100) - 0.006872*(T/100)^2))*(p0/(R*T0))*1.8e-6
# Gas Transfer Velocity

@physicsfn k(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = η*(ϵ(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
@physicsfn k_u(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = η*(ϵ_u(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
@physicsfn k_B(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = η*(ϵ_B(p,u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)

@physicsfn Fatm(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = -k(p,u,t)*(u[3]-Ceq_ch4(p,u[5],0)*1e6)
@physicsfn Fatm_u(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = -k_u(p,u,t)*(u[3]-Ceq_ch4(p,u[5],0)*1e6)
@physicsfn Fatm_B(p::ConvectionLakePhysics{<:DefaultPhysics}, u,t) = -k_B(p,u,t)*(u[3]-Ceq_ch4(p,u[5],0)*1e6)

@physicsfn F_diff(p::ConvectionLakePhysics{<:DefaultPhysics},u,t) = max(4e-8, 2.11e-8*u[1]-1.57e-7)/12*1e12/bathymetry.at(u[1])
# Wind penetration depth: deprecated!
#########################################
# Monin-Obukhov length scale (Tedford 2014)
@physicsfn AML(p::ConvectionLakePhysics{<:DefaultPhysics},u,t) = δv(p,f_wind*forcing.wind_speed.at(t))
