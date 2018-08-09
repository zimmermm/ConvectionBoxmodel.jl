#===================================
Bathymetry Implementations
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
	ContinuousBathymetryInterface(area(min_depth), volume_above(max_depth), max_depth, area, volume_above)
end
