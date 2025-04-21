# ITDickeStates.jl

Code to prepare the Dicke State $|D^n_k\rangle$ as a MPS with ITensors.jl and ITensorMPS.jl

## Installation

## Examples

```julia
N = 4
s = siteinds("Qubit",N; conserve_number=true)
k = 2 # all hamming weight 2 states
psi = dicke(s,k)
```
