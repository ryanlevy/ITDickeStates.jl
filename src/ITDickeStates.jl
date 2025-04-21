module ITDickeStates

using ITensors: hasqns, ITensor, dag, Index, QN, hastags, hasqns
using ITensorMPS: MPS, orthogonalize

export dicke_state, dense_dicke_vector

function dense_dicke_vector(n::Int, k::Int, d::Int=2)

  # This is from @lkdvos for checking the answer
  @assert d == 2 "TBA"
  D = zeros(ntuple(Returns(d), n))

  p = ntuple(i -> i <= k ? 2 : 1, n)
  for σ_p in permutations(p)
    D[σ_p...] = 1
  end

  return normalize(D)
end

"""
eqn 3.11 from doi:10.1103/PhysRevA.110.052438 
"""
function gamma(n, k, i, j, m)
  if (k-j) ≤ n-i+1
    return sqrt(1-m+(-1)^m*(j-k)/(n-i+1))
  else
    return 0
  end
end

"""
    dicke_state(sites::Vector{<:Index}, k::Integer; ortho=true)


return a n site MPS with a uniform superposition of all states
with hamming weight k
|Dⁿₖ⟩ has bond dimension k+1
for example dicke(4,2) has 6 non-zero states of the 2^4 total states

citation: doi:10.1103/PhysRevA.110.052438 

# Optional Keyword Arguments
  - `ortho = true`: orthogonalize the constructed MPS. Recommended

# Examples
```julia
using ITensors, ITensorMPS, ITDickeStates

N = 4
s = siteinds("Qubit",N; conserve_number=true)
k = 2 # all hamming weight 2 states
psi = dicke_state(s,k)

```
"""
function dicke_state(sites::Vector{<:Index}, k::Integer; ortho=true)
  n = length(sites)
  χ = k+1

  psi = MPS(sites)
  spaces = if hasqns(sites)
    if hastags(sites, "S=1/2")
      [[QN("Sz", -2*i+1) => 1 for i in 0:k] for j in 1:(n - 1)]
    elseif hastags(sites, "Qubit")
      [[QN("Number", i) => 1 for i in 0:k] for j in 1:(n - 1)]
    end
  else
    [χ for j in 1:(n - 1)]
  end

  l = [Index(spaces[ii], "Link,l=$ii") for ii in 1:(n - 1)]

  for i in 1:n
    M = zeros(χ, 2, χ)
    # _ denotes julia index,
    for m_ in 1:2
      m = m_ - 1
      for j_ in 1:χ
        j = j_-1
        jp = j+m
        if (jp+1 > 0) && (jp+1<=χ)
          γ = gamma(n, k, i, j, m)
          M[j_, m_, jp + 1] = γ
        end
      end # j_
    end # m_

    if i==1
      psi[i] = ITensor(M[1, :, :], sites[i], dag(l[i]))
    elseif i==n
      psi[i] = ITensor(M[:, :, end], (l[i - 1]), sites[i])
    else
      psi[i] = ITensor(M, (l[i - 1]), sites[i], dag(l[i]))
    end
  end
  if ortho
    psi = orthogonalize(psi, 1)
  end
  return psi
end

end # module ITDickeStates
