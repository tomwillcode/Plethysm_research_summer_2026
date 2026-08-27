from sage.all import *
load('plethysm_engine.sage')

def detect_root(term_history, min_repeats=2):
    for k in range(1, 11):
        if len(term_history) < k * min_repeats + k:
            continue
        
        is_cycle = True
        stride_diffs = []
        
        for var_idx in range(len(term_history[0]['exponents'])):
            diffs = [
                term_history[i]['exponents'][var_idx] - term_history[i-k]['exponents'][var_idx] 
                for i in range(len(term_history)-1, k-1, -1)
            ]
            
            if not all(d >= 0 for d in diffs):
                is_cycle = False
                break
                
            if len(set(diffs[:min_repeats])) != 1:
                is_cycle = False
                break
                
            stride_diffs.append(diffs[0])
        
        if is_cycle:
            sign_match = term_history[-1]['coeff'] == term_history[-1-k]['coeff']
            sign = 1 if sign_match else -1
            return {'z_power': k, 'x_powers': stride_diffs, 'sign': sign}
            
    return None

def automated_denominator_search(partition_input, target_nvars):
    nvars = 2
    
    while nvars <= target_nvars:
        var_names = [f'x{i+1}' for i in range(nvars)] + ['z']
        R = PolynomialRing(QQ, var_names)
        z = R.gen(nvars)
        x = R.gens()[:nvars]
        
        current_den = R(1)
        current_max_deg = 20
        solved = False
        
        while not solved:
            terms_by_degree = get_series_terms(partition_input, nvars, current_max_deg, current_den)
            
            if terms_by_degree and all(len(terms) == 0 for terms in terms_by_degree.values()):
                print(f"Solved for {nvars} variables!")
                print(f"Denominator: {current_den}")
                solved = True
                break
                
            factor_found = False
            
            for pos in range(4):
                front_history = []
                back_history = []
                
                for deg in sorted(terms_by_degree.keys()):
                    deg_terms = terms_by_degree[deg]
                    if len(deg_terms) > pos:
                        front_history.append(deg_terms[pos])
                    if len(deg_terms) > pos:
                        back_history.append(deg_terms[-(pos+1)])
                
                for history in [front_history, back_history]:
                    if not history:
                        continue
                        
                    root = detect_root(history)
                    
                    if root:
                        term_prod = z**root['z_power']
                        for i, p in enumerate(root['x_powers']):
                            term_prod *= x[i]**p
                            
                        if root['sign'] == 1:
                            factor = 1 - term_prod
                        else:
                            factor = 1 + term_prod
                            
                        current_den *= factor
                        current_max_deg = 20
                        factor_found = True
                        break
                        
                if factor_found:
                    break
            
            if not factor_found:
                current_max_deg += 5
                
        nvars += 1

automated_denominator_search([2], 2)