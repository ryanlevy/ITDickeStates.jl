# ITDickeStates.jl

Code to prepare the Dicke State $|D^n_k\rangle$ as a MPS with ITensors.jl and ITensorMPS.jl

Uses construction from

Raveh, David, and Rafael I. Nepomechie. "Dicke States as Matrix Product States."
Phys. Rev. A *110*, 052438 (2024) 

DOI: https://doi.org/10.1103/PhysRevA.110.052438
arxiv: 2408.04729

## Installation

## Examples

```julia
using ITensors,ITensorMPS, ITDickeStates

N = 4
s = siteinds("Qubit",N; conserve_number=true)
k = 2 # all hamming weight 2 states
psi = dicke_state(s,k)
```
