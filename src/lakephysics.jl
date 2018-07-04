#===================================
Forcing Data
===================================#

struct ForcingInterface
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
	timestamps = forcing_df[:t]
	air_temperature = interp1d!(timestamps, forcing_df[Symbol("Tair (°C)")]+273.15)
	cloud_cover = interp1d!(timestamps, forcing_df[Symbol("cloud coverage")])
	vapour_pressure = interp1d!(timestamps, forcing_df[Symbol("vap (mbar)")])

	# wind speed
	u10 = forcing_df[Symbol("u (m/s)")]
	v10 = forcing_df[Symbol("v (m/s)")]
	wind_speed = interp1d!(timestamps, sqrt.(u10.^2+v10.^2))

	# global shortwave irradiation
	G = forcing_df[Symbol("Fsol (W/m2)")]
	C = forcing_df[Symbol("cloud coverage")]

	Fdir = (1.-C)./((1.-C)+0.5*C)
	Fdiff = (0.5*C)./((1.-C)+0.5*C)
	Alb_dir = 0.2
	Alb_diff = 0.066
	Hs = G.*Fdir*(1-Alb_dir)+G.*Fdiff*(1-Alb_diff)
	global_radiation = interp1d!(timestamps, Hs)

	ForcingInterface(air_temperature, cloud_cover, vapour_pressure, global_radiation, wind_speed)
end
precompile(forcing_from_dataset, (String,))

struct InflowInterface
	flow::Function
	temperature::Function
end

function inflow_from_dataset(path)
	# load simstrat forcing dataset
	inflow_df = CSV.read(path, delim="\t")

	##########################################
	# interpolators for individual data series
	##########################################
	timestamps = inflow_df[:t]
	flow = interp1d!(timestamps, inflow_df[Symbol("Flow")])
	temperature = interp1d!(timestamps, inflow_df[Symbol("Temp")]+273.15)

	InflowInterface(flow, temperature)
end
precompile(inflow_from_dataset, (String,))


#===================================
Lake Physics
===================================#

struct LakePhysicsInterface
	N2::Function
	dhdt::Function
	dHdt::Function
	Ba::Function
	ϵ::Function
	ϵ_B::Function
	ϵ_P::Function
	ϵ_Pz::Function
	L_MO::Function
	k::Function
	k_B::Function
	k_P::Function
	Fatm::Function
	Fatm_B::Function
	Fatm_P::Function
	F_diff::Function
end


# Constructor for Convective Lake Physics
function ConvectionLakePhysicsInterface(temperature_profile, salinity_profile, concentration_profile, forcing, inflow, A=0.2)
	###########################
	# Model Constants
	###########################

	# physical constants
	const a_T = 2.14e-4 # thermal expansion [K-1]
	const g = 9.80665  # gravitational acceleration [m^2 s-1]
	const β = g*a_T
	const Mw = 18.01528/1000  # molar mass of water [kg mol-1]
	const Cp = 75.375/Mw  # specific heat capacity [J kg-1 K-1]
	const Cp_air = 1005  # specific heat capacity of air [J kg-1 K-1]
	const Lv = 2.47e6  # heat vaporization of water [J kg-1]
	const σ = 5.670367e-8  # Stefan Boltzmann constant [W m-2 K-4]
	const κ = 0.4  # Von Karman constant
	const ν = 1.5e-6  # Kinematic Viscosity of Water
	const ν_a = 1.4e-5  # Kinematic Viscosity of Air [m2 s-1]

	const T0 = 273.15  # Standard Temperature [K]
	const p0 = 1e5	# Standard Pressure [Pa]
	const R = 8.314  # universal gas constant

	const rho = 998  # density of water [kg m-3]
	const rho_air = 1.25  # density of air [kg m-3]
	const p_air = 960  # air pressure [mbar]

	const γ = (Cp_air*p_air)/(0.622*Lv)

	# empirical constants
	const A_L = 0.03  # Longwave Alebedo [-]
	const a = 1.09  # calibration constant
	const η = 1/(2*pi)  # according to Lorke et al. #0.29 # calibration constant
	const B = 0.62  # Bowen coefficient [mbar K-1]
	const c1 = 8.6

	const f_wind = 0.5

	# coefficients for density calculation
	const ai = [999.8395, 6.7914e-2, -9.0894e-3, 1.0171e-4, -1.2846e-6, 1.1592e-8, -5.0125e-11]
	const bi = [0.8181, -3.85e-3, 4.96e-5]

	density(T,S) =	begin
		    			const Tgrad = T-273.15
						const Ti = [Tgrad^i for i in 0:6]
						
						# approximate density
						dot(ai, Ti)+S*dot(bi, Ti[1:3])
					end

	# Water column stability: Brunt-väisälä frequency
	#########################
	# N2 based on gradient of initial temperature profile
	density_profile(depth) = density(temperature_profile(depth),salinity_profile(depth))
	N2_ρ(depth) = g/density_profile(depth)*derivative(density_profile, depth)

	# Heat flux at the surface
	##########################
	# positive flux: warming
	# negative flux: cooling
	#*** Longwave
	# emissivity
	Ea(Ta, cloud_cover, vapour_pressure) = (1+0.17*cloud_cover^2)*1.24*(vapour_pressure/Ta)^(1/7)
	# atmospheric longwave radiation
	Ha_simstrat(Ta, cloud_cover, vapour_pressure) = a*(1-A_L)*Ea(Ta, cloud_cover, vapour_pressure)*σ*Ta^4
	# water longwave radiation
	Hw_simstrat(Tw) = -0.972*σ*Tw^4
	#*** Shortwave
	# is given by forcing.global_radiation(t) (already corrected for albedo)
	#*** evaporation & condensation
	fu_simstrat(U10, Tw, Ta) = 4.4+1.82*U10^2+0.26*(Tw-Ta)
	e_s_simstrat(Tw, Ta)=10^((0.7859+0.03477*(Tw-273.15))/(1+0.00412*(Tw-273.15)))*(1+1e-6*p_air*(4.5+0.00006*(Tw-273.15)^2))
	He_simstrat(U10, Tw, Ta, vapour_pressure) = -fu_simstrat(U10, Tw, Ta)*(e_s_simstrat(Tw, Ta)-vapour_pressure)
	#*** sensible heat
	Hc_simstrat(U10, Tw, Ta) = -B*fu_simstrat(U10, Tw, Ta)*(Tw-Ta)

	# total heat balance [W m-2]
	heat_flux_simstrat(U10, Tw, Ta, global_radiation, cloud_cover, vapour_pressure) = Hc_simstrat(U10, Tw, Ta) + He_simstrat(U10, Tw, Ta, vapour_pressure)+Ha_simstrat(Ta, cloud_cover, vapour_pressure)+Hw_simstrat(Tw)+global_radiation
	# wrapper
	heat_flux(u, t) = heat_flux_simstrat(f_wind*forcing.wind_speed(t), u[5], forcing.air_temperature(t), forcing.global_radiation(t), forcing.cloud_cover(t), forcing.vapour_pressure(t))
	dTdt(u,t) = heat_flux(u,t)/(rho*Cp)

	# Buoyancy foring: Convective thickening of the mixed layer
	##########################################
	# Buoyancy Flux [m2 s-3]
	buoyancy_flux(u,t) = -β*dTdt(u,t)
	# thickening rate
	# Zilitinkevich 1991
	dhdt(u, t) = begin
					B0 = buoyancy_flux(u,t)
					# no thermocline erosion during warming
					# (we also assume that the thermocline is not rising anymore)
					if B0<0 | (u[5]>temperature_profile(u[1]))
						return 0
					else
						#println(N2_ρ(u[1]))
						v=(1+2*A)*B0/(N2_ρ(u[1])*u[1])
						# don't allow rising of the thermocline
						if v < 0
							return 0
						else
							return v
						end
					end
				end

	# Buoyancy
	B_w(B0,hmix)=(B0*hmix)^(1/3)

	# Wind forcing 
	##########################
	# Simstrat 1.6
	CD(U10) =	begin
					if U10 <= 0.1
						0.06215
					elseif U10 <= 3.85
						0.0044*U10^(-1.15)
					else #Polynomial approximation of Charnock's law
						-0.000000712*U10^2+0.00007387*U10+0.0006605
					end
				end
	τ(U10) = rho_air*CD(U10)*U10^2
	# wind shear (friction velocity)
	u_w(U10) = sqrt.(τ(U10)/rho)
	# thickness of diffusive boundary layer
	δv(U10)=c1*ν/u_w(U10)

	# Dissipation
	###############

	# Read et al.
	ϵ_read_u_z(u,t,z)= u_w(f_wind*forcing.wind_speed(t))^3/(κ*z)
	ϵ_read_u(u,t)= ϵ_read_u_z(u,t,AML(u,t))
	ϵ_read_B(u,t)=	begin 
				B0 = buoyancy_flux(u,t)
				if B0 < 0.
					0.
				else
					0.5*B0
				end
			end

	ϵ_read(u,t) = ϵ_read_u(u,t)+ϵ_read_B(u,t)

	# Air/Water transfer velocity
	#############################

	# Schmitt-Number
	Sc_ch4_Wanninkhof(T) = 1909.4 - 120.78*(T-273.15) + 4.1555*((T-273.15)^2) - 0.080578*((T-273.15)^3) + 0.00065777*((T-273.15)^4)
	Ceq_ch4(T,S) = exp(-68.8862 + 101.4956*(100/T) + 28.7314*log(T/100) + S*(-0.076146 + 0.043970*(T/100) - 0.006872*(T/100)^2))*(p0/(R*T0))*1.8e-6
	# Gas Transfer Velocity
	k_ch4°(ϵ) = (u,t) -> η*(ϵ(u,t)*ν)^0.25*Sc_ch4_Wanninkhof(u[5])^(-0.5)
	Fatm°(k) = (u,t) -> -k(u,t)*(u[3]-Ceq_ch4(u[5],0)*1e6)

	
	# Wind penetration depth
	#########################################
	# Monin-Obukhov length scale (Tedford 2014)
	#L_MO(u,t) = u_w(f_wind*forcing.wind_speed(t))^3/(κ*buoyancy_flux(u,t))
	L_MO(u,t,thres) = u_w(f_wind*forcing.wind_speed(t))^3/(κ*thres)
	L_MO(u,t) = L_MO(u,t,1e-8)
	AML(u,t) = δv(f_wind*forcing.wind_speed(t))
	

	ϵ_model = ϵ_read
	ϵ_model_u = ϵ_read_u
	ϵ_model_u_z = ϵ_read_u_z
	ϵ_model_B = ϵ_read_B

	k_model = k_ch4°(ϵ_model)
	k_model_u = k_ch4°(ϵ_model_u)
	k_model_B = k_ch4°(ϵ_model_B)
	Fatm_model = Fatm°(k_model)
	Fatm_model_u = Fatm°(k_model_u)
	Fatm_model_B = Fatm°(k_model_B)
	#F_diff(u,t) = 1.e-6*1000000/(16-u[1])
	F_diff(u,t) = 2.89e-8*derivative(concentration_profile, u[1])

	LakePhysicsInterface(N2_ρ, dhdt, dTdt, buoyancy_flux, ϵ_model, ϵ_model_B, ϵ_model_u, ϵ_model_u_z, L_MO, k_model, k_model_B, k_model_u, Fatm_model, Fatm_model_B, Fatm_model_u, F_diff)
end
precompile(ConvectionLakePhysicsInterface, (Function, Function, Function, ForcingInterface, Float64))

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
