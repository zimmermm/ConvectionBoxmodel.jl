#===================================
Growth Model Implementations
===================================#
abstract type MOBGrowthModel end

struct NoGrowth <: MOBGrowthModel end

struct MOBMonodModel <: MOBGrowthModel
	Vmax::Float64
	Km::Float64
	y::Float64
end

mox(model::MOBMonodModel, C::Float64) = model.Vmax*(C/(model.Km+C))
μ(model::MOBMonodModel, C::Float64) = model.y*mox(model)

mox(model::NoGrowth, C::Float64) = 0
μ(model::NoGrowth, C::Float64) = 0

scale(model::MOBMonodModel, f_vmax, f_Km, f_y) = MOBMonodModel(model.Vmax*f_vmax, model.Km*f_Km, model.y*f_y)
montecarlo_shuffle(model::MOBMonodModel) = scale(model, 1.+randn()/164., 1.+randn()/164., 1+randn()/164.)
