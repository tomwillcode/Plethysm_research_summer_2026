import time
import sys
import resource
from sage.all import *

NVARS = 5
TARGET_PARTITION = [2, 2, 1]
ESTIMATED_FINAL_DEGREE = 46

LOG_FILE = "progress_log_P221_5var.txt"
OUTPUT_FILE = "FINAL_NUMERATOR_P221_5var.txt"

BR = QQ[",".join("x"+str(i) for i in range(1, NVARS+1))+",z"]
x1, x2, x3, x4, x5, z = BR.gens()
BR.inject_variables()

Sym = SymmetricFunctions(QQ)
s = Sym.s()

def get_ram_mb():
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
    m1 = z*x1**5
    m2 = z*x1**3*x2**2
    m3 = -z**2 * x1**5 * x2**5 
    m4 = z * x1**4 * x2
    m5 = z**2 * x1**5 * x2**5
    m6 = -z * x1**3 * x2 * x3
    m7 = -z**3 * x1**5 * x2**5 * x3**5
    m8 = z * x1**2 * x2**2 * x3
    m9 = z**3 * x1**5 * x2**5 * x3**5
    m10 = z**2 * x1**6 * x2**2 * x3**2
    m11 = -z**2 * x1**4 * x2**3 * x3**3
    m12 = -z * x1**2 * x2**2 * x3
    m13 = z**4 * x1**5 * x2**5 * x3**5 * x4**5
    m14 = -z**2 * x1**4 * x2**4 * x3 * x4
    m15 = z**3 * x1**4 * x2**4 * x3**4 * x4**3
    m16 = z * x1**2 * x2 * x3 * x4
    m17 = -z**2 * x1**3 * x2**3 * x3**2 * x4**2
    m18 = -z * x1**2 * x2 * x3 * x4
    m19 = 0
    m20 = 0
    m21 = z**2 * x1**3 * x2**3 * x3**2 * x4**2
    m22 = z**2 * x1**3 * x2**3 * x3**3 * x4
    denominator_string = (1 - m1)*(1 - m2)*(1 - m3)**2*(1 - m4)*(1 - m5)*(1 - m6)*(1 - m7)**2*(1-m8)**2*(1-m9)*(1-m10)*(1-m11)*(1-m12)*(1-m13)*(1-m14)*(1-m15)*(1 - m16)*(1 - m17)**2*(1-m18)**2*(1-m19)**2*(1-m20)*(1-m21)*(1-m22)*(1 - z*x1*x2*x3*x4*x5)
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

with open(LOG_FILE, "a") as f:
    f.write(f"\n--- Starting Server Run for P_{TARGET_PARTITION} from scratch ---\n")

try:
    get_den_expanded()
except Exception as e:
    sys.exit(1)

for d in range(0, 100): 
    start_time = time.time()
    
    CC = calc_num(TARGET_PARTITION, d)
    
    elapsed = time.time() - start_time
    current_ram = get_ram_mb()
    
    if CC:
        CC_list = list(CC)
        if len(CC_list) > 6:
            front = CC_list[:3]
            back = CC_list[-3:]
            log_print(f"d={d:02d} | Terms: {len(CC_list):4d} | RAM: {current_ram:.0f} MB | Time: {elapsed:.2f}s\n  FRONT: {front}\n  BACK: {back}")
        else:
            log_print(f"d={d:02d} | Terms: {len(CC_list):4d} | RAM: {current_ram:.0f} MB | Time: {elapsed:.2f}s\n  {CC_list}")
        out += z**d * CC
        
        if (len(CC_list) <= 1) and (d > 0):
            with open(OUTPUT_FILE, "w") as f:
                f.write(f"Numerator_5var = {out}\n")
            break
    else:
        with open(OUTPUT_FILE, "w") as f:
            f.write(f"Numerator_5var = {out}\n")
        break