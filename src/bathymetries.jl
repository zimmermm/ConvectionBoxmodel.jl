#===================================
Bathymetry Implementations
===================================#
# Datastructure is a (depth, area) Interpolation

# additional functionality
surface_area(bathymetry::Interpolation{<:Any}) = bathymetry.y[1]
volume_above(bathymetry::Interpolation{<:Any},depth; dh=0.01) = begin
									grid = equigrid(top(bathymetry), depth, dh=dh)
									areas_i = bathymetry.at(grid)
									sum((areas_i[1:end-1]+areas_i[2:end])/2.)*dh
								end
total_volume(bathymetry::Interpolation{<:Any}; dh=0.01) = volume_above(bathymetry, bottom(bathymetry), dh=dh)


# helper
function equigrid(x0,xend; dh=0.01)
	n = round((xend-x0)/dh)
	dh = (xend-x0)/n
	collect(x0:dh:xend)
end
