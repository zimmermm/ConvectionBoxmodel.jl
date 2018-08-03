
# event definition
function fully_mixed()
	condition(u,t,integrator) = begin
									if u[1]>16.
										return 1000.
									else
										return u[3]-0.1
									end
								end
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=false)
end

function mixed_to(termination_depth)
	condition(u,t,integrator) = termination_depth-u[1]
	affect!(integrator) = terminate!(integrator)
	ContinuousCallback(condition,affect!,rootfind=true)
end
