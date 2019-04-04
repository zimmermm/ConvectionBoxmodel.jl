#===================================
Profile shapes library
===================================#

# predefined mathematical shapes
function FermiProfile(Qtop, Qbottom, z_interface, w_interface)
	depth -> Qbottom+(Qtop-Qbottom)./(1.0+exp.((depth-z_interface)/w_interface))
end

function GompertzProfile(Qtop, Qbottom, z_interface, w_interface)
	depth -> Qtop-(Qtop-Qbottom).*exp.(-20.0/w_interface.*exp.(-w_interface.*(depth-z_interface))./depth)
end

function StepProfile(Qtop, Qbottom, z_bottom, z_interface, sharpness)
	slope = (Qbottom-Qtop)/log(1.0+exp(sharpness*(z_bottom-z_interface)))
	depth -> Qtop+slope.*log.(1.0+exp.(sharpness*(depth-z_interface)))
end

# read from tabular data
function DataProfile(path)
	profile_data = CSV.read(path, types=[Float64, Float64], allowmissing=:none)
	LinearInterpolation(profile_data[Symbol("Depth")], profile_data[Symbol("Value")])
end
precompile(DataProfile, (String,))