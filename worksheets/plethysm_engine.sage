from sage.all import *

SymmetricFunctions(QQ).inject_shorthands(verbose=False)

def do_P_k(k, n, BR, nvars):
    CR = s[1].expand(nvars).parent()
    return s[n].expand(nvars).subs({CR.gens()[i] : BR('x'+str(i+1))**k for i in range(nvars)})

def do_P_lambda(la, n, BR, nvars):
    a = BR(expand(mul(do_P_k(p, n, BR, nvars) for p in la) * mul(BR('x'+str(i+1))-BR('x'+str(j+1)) for i in range(nvars) for j in range(i+1, nvars))))
    return sum(c * mul(BR('x'+str(i+1))^(v[i]-(nvars-i-1)) for i in range(nvars)) for (v,c) in a.dict().items() if all(v[i] > v[i+1] for i in range(nvars-1)))

def get_series_terms(partition_input, nvars, max_degree, current_den):
    BR = current_den.parent()
    z = BR.gen(nvars)
    
    terms_by_degree = {}
    
    for d in range(max_degree + 1):
        CC = sum(current_den.coefficient({z: d-r}) * do_P_lambda(Partition(partition_input), r, BR, nvars) for r in range(d+1))
        
        terms_by_degree[d] = []
        if CC:
            for coeff, monomial in list(CC):
                exponents = monomial.exponents()[0][:-1]
                terms_by_degree[d].append({'coeff': coeff, 'exponents': exponents})
                
    return terms_by_degree