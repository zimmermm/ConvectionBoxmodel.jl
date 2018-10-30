#============================================
* Data structures
* System of ODE
* Solver

I do not follow Julias intended code style and use structures as kind of objects with properties and functions
In future versions I might change this.. 
=============================================#


#===================================
Lake Model
===================================#

# Interface to user defined functions of the model
struct LakeModel
	bathymetry::Interpolation{<:Any}

	initial_condition::Array{Float64,1}

	# static profiles
	lake_temperature::Interpolation{<:Any}
	concentration_profile::Interpolation{<:Any}

	lake_physics::LakePhysics

	growth_model::MOBGrowthModel

	starttime::Float64
	endtime::Float64
	model_callback::ContinuousCallback
end


#===================================
Boxmodel
===================================#

function boxmodel_ode(du,u,lakemodel,t)
	h_mix, V, C, B, T, Qatm, Qmix, Qmox, Qdiff = u

	# include morphology to support non-equidistant time steps
	du[1] = dhdt(lakemodel.lake_physics, u, t)*3600.0*24.0
	du[2] = du[1]*lakemodel.bathymetry.at(h_mix)

	# Flux to the atmosphere
	du[6] = -Fatm(lakemodel.lake_physics, u, t)*surface_area(lakemodel.bathymetry)*3600.0*24.0
	# Flux from hypolimnion into mixed layer
	du[7] = du[2]*lakemodel.concentration_profile.at(h_mix)+F_diff(lakemodel.lake_physics, u,t)*lakemodel.bathymetry.at(h_mix)*3600.0*24.0
	# MOX
	du[8] = 0.0
	du[9] = F_diff(lakemodel.lake_physics, u,t)*lakemodel.bathymetry.at(h_mix)*3600.0*24.0

	# basic ODE for expanding box size
	du[3] = 1.0/V*(du[7]-du[2]*C-du[6])
	du[4] = -1.0/V*du[2]*B
	du[5] = 1.0/V*du[2]*(lakemodel.lake_temperature.at(h_mix)-T)+1.0/V*dTdt(lakemodel.lake_physics, u,t)*surface_area(lakemodel.bathymetry)*3600.0*24.0

	mu = [0.0, 0.0, -mox(lakemodel.growth_model, C), μ(lakemodel.growth_model, C), 0.0, 0.0, 0.0, mox(lakemodel.growth_model, C), 0.0]*3600.0*24.0*B
	du[:] = du+mu
	du[8] = du[8]*V
	return
end
precompile(boxmodel_ode, (Array{Float64,1}, Array{Float64,1}, LakeModel, Float64))


#===================================
Solver
===================================#

function solve_boxmodel(lakemodel; saveat=[])
	# assemble initial conditions
	h_mix0, C0, B0, T0 = lakemodel.initial_condition
	V0 = volume_above(lakemodel.bathymetry, h_mix0)
	# state variables
	u0 = [h_mix0, V0, C0, B0, T0, 0.0, 0.0, 0.0, 0.0]

	# tspan
	tspan = (lakemodel.starttime,lakemodel.endtime)

	# ODEProblem
	prob = ODEProblem(boxmodel_ode,u0,tspan,lakemodel,callback=lakemodel.model_callback)

	# solve
	if isempty(saveat)
		@time solve(prob, Rosenbrock23(autodiff=false), reltol=1.0e-2, abstol=1.0e-2, dtmax=1.0/24.0)
	else
		@time solve(prob, Rosenbrock23(autodiff=false), reltol=1.0e-2, abstol=1.0e-2, dtmax=1.0/24.0, saveat=saveat)
	end
end
#precompile(solve_boxmodel, (LakeModel,))


function solve_boxmodel_montecarlo(lakemodel, num_monte)
	# assemble initial conditions
	h_mix0, C0, B0, T0 = lakemodel.initial_condition
	V0 = volume_above(lakemodel.bathymetry, h_mix0)
	# state variables
	u0 = [h_mix0, V0, C0, B0, T0, 0.0, 0.0, 0.0, 0.0]

	# tspan
	tspan = (lakemodel.starttime,lakemodel.endtime)

	# ODEProblem
	prob = ODEProblem(boxmodel_ode,u0,tspan,lakemodel,callback=lakemodel.model_callback)
	montecarlo = MonteCarloProblem(	prob,
									#output_func = (sol, i) -> (sol, false),
									prob_func = boxmodel_prob_func
									#reduction = (u, data, I) -> (append!(u,data),false),
									#u_init = u0
								)
	# solve
	@time solve(montecarlo, Rosenbrock23(autodiff=false), reltol=1.0e-2, abstol=1.0e-2, dtmax=1.0/24.0, num_monte=num_monte, parallel_type=:none)
end
precompile(solve_boxmodel_montecarlo, (LakeModel,))


function boxmodel_prob_func(prob, i, repeat)
	p = prob.p
	lakemodel = LakeModel(
								p.bathymetry,						# bathymetry
								p.initial_condition,						# initial_conditions
								p.lake_temperature,						# temperature profile
								p.concentration_profile,						# concentration profile
								p.lake_physics,
								montecarlo_shuffle(p.growth_model),
								p.starttime,
								p.endtime,
								p.model_callback
							)
	u0 = [
			(1.0+randn()/13.0)*prob.u0[1],  # hmix0
			volume_above(prob.p.bathymetry, prob.u0[1]),  # V0
			(1.0+randn()/7.0)*prob.u0[3],  # C0
			(1.0+randn()/8.0)*prob.u0[4],  # B0
			(1.0+randn()/244.0)*prob.u0[5],  # T0
			0.0, 0.0, 0.0, 0.0
		]
	ODEProblem(prob.f, u0, prob.tspan, lakemodel, callback=prob.callback)
end
precompile(boxmodel_prob_func, (ODEProblem, Int64, Int64))
