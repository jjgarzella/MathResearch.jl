#using Oscar

using UUIDs

# this is duplicated in another module... fix this later
exponent_vectors(poly) = leading_exponent_vector.(terms(poly))

function projective_space_fan_string(n)
  
  projstring = "[[" * join(fill(-1,n-1)," ") * "]"

  for i in n-1:-1:1
    coordarray = zeros(Int,n-1)
    coordarray[i] = 1
    
    projstring = projstring * "[" * join(coordarray, " ") * "]"
  end

  projstring = projstring * "]"

  projstring
end#function
"""
outputs the exponent vector of f in the form used
by controlled reduction scripts

Note: use trim=true for ToricControlledReduction,
use trim=false for controlledreduction
"""
function exponent_vector_string(f;trim=true)
  ev_string = ""
  ev_string *= "["

  for exp_vec in exponent_vectors(f)
    trimfactor = trim ? 1 : 0
    ev_string *= "[" * join(exp_vec[1:end-trimfactor], " ") * "]"
  end

  ev_string *= "]"

  ev_string
end

coeffs_string(f) = "[" * join(coefficients(f)," ") * "]"

function tcr_string(p,f,name)
  tcrstring = name * ":"

  n = length(gens(parent(f)))

  # put the exponent vectors here
  
  tcrstring *= exponent_vector_string(f)

  tcrstring *= ":"

  # put the coefficients here
  
  tcrstring *= coeffs_string(f)
  tcrstring *= ":"

  # The fan for projective space

  tcrstring *=  projective_space_fan_string(n) * ":"

  # the degree here

  degarr = zeros(Int,n)
  deg = sum(exponent_vectors(f)[1]) # we assume homogeneous
  degarr[1] = deg
  
  tcrstring *= "[" * join(degarr," ") * "]"

  # the prime number

  tcrstring *= ":" * string(p)

  # println(tcrstring)

  tcrstring
end#function

function cr_string(p,f)
  crstring = "$p\n"
  crstring *= exponent_vector_string(f,trim=false) * "\n"
  crstring *= coeffs_string(f)

  crstring
end

#
#  # similar to tcr_string
#end

function output_tcr_file(ps,polys,names,outputfilename)
  open(outputfilename,"w") do outfile
    for i in 1:length(ps)
      println(outfile,tcr_string(ps[i],polys[i],names[i]))
    end
  end
end#function

"""
assumes julia is running from the project directory
"""
function setup_tcr()
  run(`./configure_tcr.bash`)
end

"""
assumes julia is running from the project directory
"""
function zeta_function_tcr(p,polynomial)
  zeta_function_tcr([p],[polynomial])
end

"""
assumes julia is running from the project directory

assumes all polynomials are homogeneous
"""
function zeta_function_tcr(ps::Array,fs::Array)
  unique_id = string(uuid1())
  l = length(ps)
  names = fill("",l)
  for i in 1:l
    d = degree(f,1) # assumes homogeneous
    n = length(gens(parent(fs[i])))
    names[i] = "example$(i)_dim$(n-1)_degree$(d)_" * unique_id
  end
  filename = "data/examples_" * unique_id
  output_tcr_file(ps,fs,names,filename)


  run(`./run_tcr.bash $filename`)
  run(`cat $filename.out`)
end
