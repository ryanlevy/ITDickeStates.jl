module ITDickeStates

using ITensors: hasqns, ITensor, dag, Index, QN, hastags, hasqns, @SiteType_str, SiteType
using ITensors.SiteTypes: _sitetypes
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

# This produces link indices with the proper QN labels
linear_linkinds(::SiteType, n::Int, k::Int) = nothing

function linear_linkinds(::SiteType"Qubit", n::Int, k::Int)
  return [[QN("Number", i) => 1 for i in 0:k] for j in 1:(n - 1)]
end

function linear_linkinds(::SiteType"S=1/2", n::Int, k::Int)
  return [[QN("Sz", -2*i+1) => 1 for i in 0:k] for j in 1:(n - 1)]
end

function linear_linkinds(sites::Vector{<:Index}, n::Int, k::Int)
  stypes = _sitetypes(first(sites))
  for st in stypes
    linds = linear_linkinds(st, n, k)
    !isnothing(linds) && return linds
  end
  error("unknown SiteType link inds construction, please define your own")
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


return a n site MPS representing the Dicke state |Dⁿₖ⟩, a uniform superposition of all states
with hamming weight k. The MPS has bond dimension k+1.

Example: |D⁴₂⟩ has 6 non-zero states of the 2^4 total states

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
    linear_linkinds(sites, n, k)
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
