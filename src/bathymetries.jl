#==================================================================
Bathymetry Implementations

The Datastructure of a Bathymetry is a (depth, area) Interpolation
===================================================================#

const Bathymetry{InterpolationMethod} = Interpolation{InterpolationMethod}

# additional functionality related to bathymetry
surface_area(bathymetry::Bathymetry{<:Any}) = bathymetry.y[1]
volume_above(bathymetry::Bathymetry{<:Any},depth; dh=0.01) = begin
									grid = equigrid(top(bathymetry), depth, dh=dh)
									areas_i = bathymetry.at(grid)
									sum((areas_i[1:end-1]+areas_i[2:end])/2.)*dh
								end
total_volume(bathymetry::Bathymetry{<:Any}; dh=0.01) = volume_above(bathymetry, bottom(bathymetry), dh=dh)
area_at(bathymetry::Bathymetry{<:Any}, depth) = bathymetry.at(depth)


# helper
function equigrid(x0,xend; dh=0.01)
	n = round((xend-x0)/dh)
	dh = (xend-x0)/n
	collect(x0:dh:xend)
end