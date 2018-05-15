
# event definition
function fully_mixed()
	condition(u,t,integrator) = 16.-u[1]
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=true)
end

function mixed_to(termination_depth)
	condition(u,t,integrator) = termination_depth-u[1]
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=true)
end
