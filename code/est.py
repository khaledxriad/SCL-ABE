"""Parameter solver + primal-uSVP core-SVP estimator (ADPS16-style heuristic)."""
import math

def delta_beta(b):
    return ((math.pi*b)**(1.0/b) * b/(2*math.pi*math.e))**(1.0/(2*(b-1)))

def primal_beta(n, q, sigma, mmax):
    """smallest beta for primal uSVP, optimizing the number of samples m<=mmax."""
    best = None
    lq = math.log(q)
    for m in range(max(n//2, 50), mmax+1, max(1, (mmax-n//2)//60 or 1)):
        d = m + n + 1
        # binary-search beta
        lo, hi = 50, d
        ok = False
        while lo < hi:
            b = (lo+hi)//2
            # log of RHS: (2b-d-1)*log(delta) + (m/d)*log q
            rhs = (2*b-d-1)*math.log(delta_beta(b)) + (m/d)*lq
            lhs = 0.5*math.log(b) + math.log(sigma)
            if lhs <= rhs: hi = b
            else: lo = b+1
        b = lo
        if best is None or b < best: best = b
    return best

def bits_classical(beta): return 0.292*beta
def bits_quantum(beta): return 0.265*beta

C_R = 12.0      # ABB10 spectral-norm constant for +-1 matrices (s1(R) <= 12 sqrt(2m))
LAMBDA_STAT = 128
ETA = math.sqrt(math.log(2*(1+2.0**LAMBDA_STAT))/math.pi)   # smoothing parameter of Z for eps=2^-lambda (~5.33)
TAU = math.sqrt(2*math.log(2*2.0**LAMBDA_STAT))              # sub-Gaussian tail multiplier, delta=2^-lambda (~13.4)

def derived(n, k, L, t_bits, safety=1.0):
    """All derived quantities for MP12 gadget trapdoor (base 2, q ~ 2^k)."""
    mbar = n*k + 256
    w = k
    m = mbar + n*w
    sigR = ETA
    s1R = sigR*(math.sqrt(mbar)+math.sqrt(n*w)+6)
    sigG = math.sqrt(5)*ETA
    s_mp = math.sqrt(sigG**2*(s1R**2+1) + ETA**2)               # MP12 SamplePre / SampleLeft requirement
    rho = (L+1)*C_R*math.sqrt(2*m)                              # worst-case s1(R_v), R_v = sum of L+1 +-1 matrices
    s_sr = math.sqrt(5)*ETA*rho                                  # SampleRight requirement used in the security proof
    s = 1.02*max(s_mp, s_sr)                                     # key Gaussian width
    # correctness: noise <= alpha*q*s*( C_R*sqrt2*(L+1)*m^1.5 + m )  must be < q/(2p*safety/ ... )
    p = 2**t_bits
    tau = TAU
    # |e2^T R'^T y| <= 1.1*tau*sqrt((L+1)/(2 pi)) * s*m*alpha*q ;  |e1^T y|+|x| negligible extra (+tau/sqrt(2pi) s sqrt(m))
    coef = tau*math.sqrt(L+1)*s*m + tau*s*math.sqrt(m) + tau
    alpha = 1.0/(2*p*safety*coef)
    return dict(n=n,k=k,m=m,mbar=mbar,s1R=s1R,sigG=sigG,s=s,alpha=alpha,coef=coef)

def solve(n, L, t_bits, kmin=20, kmax=200, alpha_q_min=None):
    """smallest k such that alpha*q >= max(2 sqrt(n), 8)"""
    for k in range(kmin, kmax):
        d = derived(n, k, L, t_bits)
        aq = d['alpha']*2.0**k
        need = alpha_q_min or max(2*math.sqrt(n), 8.0)
        if aq >= need:
            d['aq'] = aq
            return d
    return None

def choose(N, lam=128, t_bits=16, n_grid=None, quantum=False):
    L = int(math.ceil(math.log2(2*N-1)))
    loss_bits = math.log2(2*N-1)
    target = lam + loss_bits
    n_grid = n_grid or range(256, 4097, 32)
    for n in n_grid:
        d = solve(n, L, t_bits)
        if d is None: continue
        sig = d['alpha']*2.0**d['k']/math.sqrt(2*math.pi)
        beta = primal_beta(n, 2.0**d['k'], sig, mmax=min(d['m']+16, 6*n))
        bits = (bits_quantum if quantum else bits_classical)(beta)
        if bits >= target:
            d.update(L=L, beta=beta, bits=bits, loss_bits=loss_bits, target=target, sigma_e=sig)
            return d
    return None

if __name__ == "__main__":
    # sanity: Kyber-512-like LWE: n=512,q=3329,sigma=sqrt(1.5) -> expected beta ~ 400, 118 bits core-SVP
    b = primal_beta(512, 3329, math.sqrt(1.5), 512+0)
    print("Kyber512-like sanity: beta=%d classical bits=%.1f" % (b, bits_classical(b)))
    for N in [2**6, 2**10, 2**16]:
        d = choose(N)
        print(N, d)
