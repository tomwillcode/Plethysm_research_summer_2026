import time
import sys
import os
import resource
import numpy as np
from sage.all import *

# ==========================================
# CONFIGURATION
# ==========================================
NVARS = 5
TARGET_PARTITION = [2, 1, 1, 1]
ESTIMATED_FINAL_DEGREE = 46 # We know the 5-var matches the 4-var numerator bounds

LOG_FILE = "progress_log_P2111_5var.txt"
OUTPUT_FILE = "FINAL_NUMERATOR_P2111_5var.txt"

# Define Ring
BR = QQ[",".join("x"+str(i) for i in range(1, NVARS+1))+",z"]

# EXPLICIT UNPACKING: This prevents NameErrors when running headless via nohup
x1, x2, x3, x4, x5, z = BR.gens()
BR.inject_variables()

# EXPLICIT DEFINITION: Prevents NameErrors for 's' in headless scripts
Sym = SymmetricFunctions(QQ)
s = Sym.s()

def get_ram_mb():
    # Returns the max resident set size (RAM usage) in MB
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

def log_print(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

def do_P_k(k, n):
    global BR, NVARS
    CR = s[1].expand(NVARS).parent()
    return s[n].expand(NVARS).subs({CR.gens()[i] : BR('x'+str(i+1))**k for i in range(NVARS)})

@cached_function
def do_P_lambda(la, n):
    global BR, NVARS
    a = BR(expand(mul(do_P_k(p, n) for p in la) * mul(BR('x'+str(i+1)) - BR('x'+str(j+1)) for i in range(NVARS) for j in range(i+1, NVARS))))
    return sum(c * mul(BR('x'+str(i+1))^(v[i] - (NVARS-i-1)) for i in range(NVARS)) \
               for (v,c) in a.dict().items() if all(v[i] > v[i+1] for i in range(NVARS-1)))

@cached_function
def den_guess():
    global BR, x1, x2, x3, x4, x5, z
    
    # =======================================================
    # PASTE YOUR EXACT FACTORED 5-VARIABLE DENOMINATOR BELOW
    # =======================================================
    m1 = z*x1**5
    m2 = z**2*x1**5*x2**5
    m3 = z*x1**4*x2 
    m4 = -z**2*x1**5*x2**5
    m5 = z*x1**3*x2**2
    m6 = z*x1**3*x2*x3
    m7 = -z*x1**3*x2*x3
    m8 = z**2*x1**4*x2**4*x3**2
    m9 = z*x1**2*x2**2*x3
    m10 = -z**2*x1**4*x2**3*x3**3
    m11 = z**2*x1**4*x2**3*x3**3
    m12 = -z**3*x1**5*x2**5*x3**5
    m13 = z**3*x1**5*x2**5*x3**5
    m14 = -z**4*x1**5*x2**5*x3**5*x4**5
    m15 = z**3*x1**4*x2**4*x3**4*x4**3
    m16 = -z*x1**2*x2*x3*x4
    m17 = z*x1**2*x2*x3*x4
    m18 = z**2*x1**3*x2**3*x3**2*x4**2
    m19 = -z**2*x1**3*x2**3*x3**2*x4**2
    m20 = -z**2*x1**4*x2**4*x3*x4
    m21 = z**2*x1**4*x2**4*x3*x4
    m22 = -z**2*x1**3*x2**3*x3**3*x4
    m23 = -z*x1*x2*x3*x4*x5
    
    # PASTE YOUR FULL COMBINED DENOMINATOR STRING HERE:
    # Example format: (1 - m1) * (1 - m2) ...
    denominator_string = (1-m1)*(1-m2)**2*(1-m3)*(1-m4)*(1-m5)*(1-m6)**2*(1-m7)*(1-m8)*(1-m9)*(1-m10)*(1-m11)*(1-m12)**2*(1-m13)*(1-m14)*(1-m15)*(1-m16)**2*(1-m17)*(1-m18)*(1-m19)**2*(1-m20)*(1-m21)*(1-m22)*(1-m23)
    
    return BR(denominator_string)

EXPANDED_DENOMINATOR = None

def get_den_expanded():
    global EXPANDED_DENOMINATOR
    if EXPANDED_DENOMINATOR is None:
        EXPANDED_DENOMINATOR = den_guess()
    return EXPANDED_DENOMINATOR

@cached_function
def den_coeff(d):
    global BR, z
    return get_den_expanded().coefficient({z: d})

def calc_num(la, d):
    return sum(den_coeff(d-r) * do_P_lambda(Partition(la), r) for r in range(d+1))

out = 0
time_history = []
degree_history = []

# Clear/create the log file safely
with open(LOG_FILE, "w") as f:
    f.write(f"--- Starting Server Run for P_{TARGET_PARTITION} ({NVARS} Variables) ---\n")

log_print("Warming up native denominator expansion... (takes a few seconds)")
try:
    get_den_expanded()
    log_print("Expansion complete! Starting degrees...\n")
except Exception as e:
    log_print(f"CRITICAL ERROR during denominator expansion: {e}")
    sys.exit(1)

for d in range(0, 100): 
    start_time = time.time()
    
    CC = calc_num(TARGET_PARTITION, d)
    
    elapsed = time.time() - start_time
    time_history.append(elapsed)
    degree_history.append(d)
    
    current_ram = get_ram_mb()
    
    eta_str = ""
    if d < ESTIMATED_FINAL_DEGREE and len(time_history) >= 8:
        try:
            mask = np.array(degree_history) > 6
            x_data = np.array(degree_history)[mask]
            y_data = np.array(time_history)[mask]
            
            coeffs = np.polyfit(x_data, np.log(y_data), deg=2)
            A, B, C = coeffs[0], coeffs[1], coeffs[2]
            
            remaining_degrees = np.arange(d + 1, ESTIMATED_FINAL_DEGREE + 1)
            predicted_times = np.exp(A * (remaining_degrees**2) + B * remaining_degrees + C)
            projected_remaining = np.sum(predicted_times)
            
            decay_factor = np.exp(2 * A)
            eta_str = f" | ETA to d={ESTIMATED_FINAL_DEGREE}: ~{projected_remaining/60:.1f} min (Decay: {decay_factor:.3f})"
        except Exception:
            eta_str = " | ETA: Calc..."
            
    if CC:
        CC_list = list(CC)
        if len(CC_list) > 6:
            front = CC_list[:3]
            back = CC_list[-3:]
            log_print(f"d={d:02d} | Terms: {len(CC_list):4d} | RAM: {current_ram:.0f} MB | Time: {elapsed:.2f}s{eta_str}\n  FRONT: {front}\n  BACK: {back}")
        else:
            log_print(f"d={d:02d} | Terms: {len(CC_list):4d} | RAM: {current_ram:.0f} MB | Time: {elapsed:.2f}s{eta_str}\n  {CC_list}")
        out += z**d * CC
    else:
        # We hit exactly 0! The finite numerator has terminated.
        log_print(f"d={d:02d} | Terms:    0 | RAM: {current_ram:.0f} MB | Time: {elapsed:.2f}s")
        log_print(f"\n[SUCCESS] Numerator terminated at degree {d-1}!")
        
        # Write the final result to the output file
        with open(OUTPUT_FILE, "w") as f:
            f.write(f"# Exact Numerator for P_{TARGET_PARTITION} ({NVARS}-Variable)\n")
            f.write(f"# Terminated at degree {d-1}\n\n")
            f.write(f"Numerator_5var = {out}\n")
            
        log_print(f"Saved complete numerator polynomial to {OUTPUT_FILE}.")
        break