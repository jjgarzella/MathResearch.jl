# MathResearch.jl
Various scripts to import data/interface with one-off pure math research projects in Julia.

This is intended to be "quick and dirty", so that researchers 
(especially those using [Oscar.jl](https://docs.oscar-system.org/stable/)) can use the
data/software in experiments. 
As such, not much attention is paid to software development or performance.

If you'd like to add some dataset/library that isn't already available in Julia,
contributions are welcome. 

If you'd like to take one of these projects and more formally integrate
it with a computer algebra system or large-scale database project,
such a contribution is also very welcome.
Feel free to reach out to the author if you'd like any support,
and also feel free to use the code from this repo (it's MIT licensed).

## Toric Controlled Reduction

ToricControlledReduction calculates the zeta function of toric hypersurfaces
over finite fields. See the [original Github repository](https://github.com/edgarcosta/ToricControlledReduction)
for more info. This wrapper only supports projective hypersurfaces right now.

To set up ToricControlledReduction, make sure you have Julia install,
clone this repo, and go into the root directory of the clone. 
Then open a Julia REPL and run:

```
Pkg.add("UUIDs")

# run these two lines every time you open the REPL
using Oscar 
include("src/ToriccontrolledReduction.jl") 

setup_tcr() # this will download ToricControlledReduction, you only ever need to run it once
```

Now, you can run examples like this:
```
p = 11

R, (x,y,z,w) = polynomial_ring(GF(p),4)

f = x^4 + y^4 + z^4 + w^4

zeta_function_tcr(p,f)
```

The final array printed is the array of coefficients of the zeta function,
you can copy and paste it into the Julia REPL if you want.

## A census of cubic fourfolds

TODO: insert link to repo and sample code here

## Wish list

* census of quartic K3 surfaces in characteristic 2 (Keldaya/sutherland)
