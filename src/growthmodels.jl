#===================================
Growth Model Implementations
===================================#

module MOBGrowthModelLibrary
	μ(C,vmax,Km)=vmax*(C/(Km+C))
end

# constructor for growth model
function MOBGrowthModelInterface(vmax, Km, y)
	function mutate(f)
		MOBGrowthModelInterface(f[1]*vmax, f[2]*Km, f[3]*y)
	end

	function montecarlo_shuffle()
		generate_growth_term((1+randn()/164)*vmax,
							(1+randn()/164)*Km,
							(1+randn()/164)*y
							)
	end

	function generate_growth_term(vmax, Km,y)
		function growth_term(u, t)
			h_mix, V, C, B, T, Qatm, Qmix, Qmox = u
			# potentia temperature dependence of growth rate
			#mu=MOBGrowthModelLibrary.μ(C,vmax*exp(-0.00244*u[5]),Km)
			mu=MOBGrowthModelLibrary.μ(C,vmax,Km)
			[0., 0., -mu, y*mu, 0., 0., 0., mu, 0.]
		end
		growth_term
	end
	MOBGrowthModel(mutate, montecarlo_shuffle, generate_growth_term(vmax, Km, y), generate_growth_term(vmax, Km, y))
end
precompile(MOBGrowthModelInterface, (Float64, Float64, Float64))

function NoBiomass()
	no_growth = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
	function growth_rate(u, t)
		no_growth
	end
	MOBGrowthModel(() -> growth_rate, growth_rate, growth_rate)
end
