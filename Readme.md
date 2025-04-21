# ITDickeStates.jl

Code to prepare the Dicke State $|D^n_k\rangle$ as a MPS with ITensors.jl and ITensorMPS.jl

## Installation

Because the package isn't registered, please install with
```julia
julia> using Pkg; Pkg.add(url="https://github.com/ryanlevy/ITDickeStates.jl")
```

## Background

The Dicke state $|D^n_k\rangle$ is defined by 

$$
|D^n_k\rangle \sim \sum_{w; H(w) = k} |w\rangle
$$

Or the uniform superposition of all configurations $w$ of length $n$ that have hamming weight $k$.

The MPS of this state has bond dimension $\chi = k+1$

For example[1],

$$ 
|D^4_2\rangle = \frac{1}{\sqrt{6}} (∣0011\rangle+ ∣0101\rangle + ∣0110\rangle + ∣1001\rangle + ∣1010\rangle + ∣1100\rangle)
$$

## Examples

```julia
using ITensors, ITensorMPS, ITDickeStates

N = 4
s = siteinds("Qubit",N; conserve_number=true)
k = 2 # all hamming weight 2 states
psi = dicke_state(s,k)
```

# Credits

Uses construction from
```
[1] Raveh, David, and Rafael I. Nepomechie. "Dicke States as Matrix Product States."
Phys. Rev. A 110, 052438 (2024)  
DOI: https://doi.org/10.1103/PhysRevA.110.052438   
arxiv: [2408.04729  ](https://arxiv.org/abs/2408.04729)
```

Please cite their paper if using this package

Also thanks to @lkdvos for useful discussion
