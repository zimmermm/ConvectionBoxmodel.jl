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

function DataProfile(path)
	profile_data = CSV.read(path, types=[Array{Float64}, Array{Float64}])
	LinearInterpolation(profile_data[Symbol("Depth")], profile_data[Symbol("Value")])
end
precompile(DataProfile, (String,))


#function initial_T_profile_from_simstrat_ic(path)
#	ic_df = CSV.read(path, delim="\t")
#	interp1d(-ic_df[Symbol("depth (m)")], ic_df[Symbol("T (°C)")]+273.15)
#end
#precompile(initial_T_profile_from_simstrat_ic, (String,))
#
#function initial_T_profile_from_simstrat_output(path)
#	data, header = readdlm(path, header=true)
#	df = DataFrame(data)
#	names!(df, [Symbol("$name") for name in header[1,:]])
#	df
#end
#precompile(initial_T_profile_from_simstrat_output, (String,))
