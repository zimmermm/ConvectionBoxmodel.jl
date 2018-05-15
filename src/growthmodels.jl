
#===================================
Growth Model Implementations
===================================#

module MOBGrowthModel
	μ(C,vmax,Km)=vmax*(C/(Km+C))
end

# constructor for growth model
function MOBGrowthModelInterface(vmax, Km, y)
	function growth_term(u, t)
		h_mix, V, C, B, T, Qatm, Qmix, Qmox = u
		mu=MOBGrowthModel.μ(C,vmax,Km)
		[0., 0., -mu, y*mu, 0., 0., 0., mu]
	end
end
precompile(MOBGrowthModelInterface, (Float64, Float64, Float64))

function NoBiomass()
    no_growth = [0., 0., 0., 0., 0., 0., 0., 0.]
	function growth_rate(u, t)
		no_growth
	end
end

