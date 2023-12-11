using ForwardDiff

function periogardo(prob, θ)
  cl2(x) = prob.model(x, t)
  grad_model!(h, x, t_) = begin

    global t = t_

    return ForwardDiff.gradient(h, x)
  end
  (m, n) = size(prob.data)
  F = zeros(m)
  J = zeros(m, Int(prob.dim))

  for i = 1:m
    J[i, :] = grad_model!(cl2, θ, prob.data[i, :])
    F[i] = cl2(θ)
  end

  return F, J

end
#grad_model!(cl2,θ, prob.data[1,i])