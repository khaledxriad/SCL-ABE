"""Independent re-derivation of the paper's numbers (does not import est.py/sizes.py/cover.py)."""
import json, math, itertools, random, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('SCL_DATA', os.path.join(HERE, '..', 'data'))
P = json.load(open(os.path.join(DATA, 'params.json')))
LAM = 128; T = 16; ELL = LAM // T
fails = []
def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + (' ' + info if info else ''))
    if not ok: fails.append(name)

# ---------- 1. cover: brute force against the definition ----------
def cov_def(h, S):
    """Cov = {v : Leaves(v) subset of S and (v root or Leaves(parent) not subset of S)}."""
    N = 1 << h; S = set(S); out = []
    def lv(v):
        d = v.bit_length() - 1; span = 1 << (h - d); lo = (v - (1 << d)) * span
        return set(range(lo, lo + span))
    for v in range(1, 2 * N):
        if lv(v) <= S and (v == 1 or not lv(v >> 1) <= S): out.append(v)
    return sorted(out), lv
def cov_fast(h, R):
    R = sorted(R); out = []
    def rec(v, lo, hi):
        import bisect
        i = bisect.bisect_left(R, lo)
        if i == len(R) or R[i] >= hi: out.append(v); return
        if hi - lo == 1: return
        m = (lo + hi) // 2; rec(2 * v, lo, m); rec(2 * v + 1, m, hi)
    rec(1, 0, 1 << h); return sorted(out)
random.seed(5)
ok = True
for h in (3, 4, 5):
    N = 1 << h
    for _ in range(300):
        r = random.randint(0, N - 1); R = random.sample(range(N), r); S = set(range(N)) - set(R)
        a, lv = cov_def(h, S); b = cov_fast(h, R)
        if a != b: ok = False
        # partition property
        if sorted(x for v in a for x in lv(v)) != sorted(S): ok = False
        if r >= 1 and len(a) > min(N - r, r * (math.log2(N / r) + 1)) + 1e-9: ok = False
check('cover == definition, exact partition, size bound (h=3,4,5, 900 trials)', ok)
ok = True
for h in range(2, 13):
    N = 1 << h
    for _ in range(200):
        a = random.randrange(N); b = random.randrange(a, N)          # contiguous S=[a,b]
        S = set(range(a, b + 1)); R = [x for x in range(N) if x not in S]
        if R and len(cov_fast(h, R)) > 2 * h: ok = False
        if not R and cov_fast(h, R) != [1]: ok = False
check('contiguous range needs <= 2h nodes (h=2..12)', ok)
check('r=1 gives h nodes', all(len(cov_fast(h, [random.randrange(1 << h)])) == h for h in range(1, 17)))

# ---------- 2. parameter constraints (P1)-(P6) ----------
def eta(eps_bits=LAM): return math.sqrt(math.log(2 * (1 + 2.0 ** eps_bits)) / math.pi)
def uSVP_beta(n, q, sigma):
    """smallest beta s.t. sigma*sqrt(beta) <= delta^(2beta-d-1) q^(m/d), optimised over m (GSA, primal)."""
    best = 10 ** 9
    for m in range(n, 6 * n, max(1, n // 40)):
        d = m + n + 1
        for b in range(60, d):
            dl = ((math.pi * b) ** (1 / b) * b / (2 * math.pi * math.e)) ** (1 / (2 * (b - 1)))
            if sigma * math.sqrt(b) <= dl ** (2 * b - d - 1) * q ** (m / d): best = min(best, b); break
    return best
for N_s, v in P.items():
    N = int(N_s); h = v['h']; n, k, m, mbar, s = v['n'], v['k'], v['m'], v['mbar'], v['s']
    L = h + 1
    q_bits = k; q = 2.0 ** k  # q is the prime just above 2^(k-1) .. 2^k; only its bit length enters the checks
    tag = 'N=2^%d ' % h
    check(tag + 'P1 mbar=nk+2lam, m=mbar+nk, L=ceil(log2(2N-1))', mbar == n * k + 2 * LAM and m == mbar + n * k and L == math.ceil(math.log2(2 * N - 1)))
    eta_ = eta(); sR = eta_ * (math.sqrt(mbar) + math.sqrt(n * k) + 6)
    sMP = math.sqrt(5 * eta_ ** 2 * (sR ** 2 + 1) + eta_ ** 2); sSR = math.sqrt(5) * eta_ * (L + 1) * 12 * math.sqrt(2 * m)
    check(tag + 'P2 s >= max(sMP,sSR)', s >= max(sMP, sSR) * (1 - 1e-12), 's=%.0f need=%.0f' % (s, max(sMP, sSR)))
    aq = 2.0 ** v['alpha_log2'] * q
    check(tag + 'P3 alpha*q >= 2 sqrt(n)', v['aq'] >= 2 * math.sqrt(n) - 1e-9, 'aq=%.1f' % v['aq'])
    tau = math.sqrt(2 * math.log(2 * 2.0 ** LAM)); coef = tau * math.sqrt(L + 1) * s * m + tau * s * math.sqrt(m) + tau
    alpha = 2.0 ** v['alpha_log2']; Delta = q / 2 ** T
    check(tag + 'P4 alpha*q*coef+1 <= Delta/2', alpha * q * coef + 1 <= Delta / 2 * (1 + 1e-9), 'ratio=%.3f' % ((alpha * q * coef + 1) / (Delta / 2)))
    need = LAM + math.log2(2 * N - 1)
    check(tag + 'P5 bits >= lam+log2(2N-1)', v['bits'] >= need, 'bits=%.1f need=%.1f' % (v['bits'], need))
    check(tag + 'P5 bits == 0.292*beta', abs(v['bits'] - 0.292 * v['beta']) < 0.01)
    check(tag + 'P6 nk >= k+log2(m(L+1))+3', n * k >= k + math.log2(m * (L + 1)) + 3)
    # independent estimator at the reported noise
    sigma = alpha * q / math.sqrt(2 * math.pi)
    if h in (6, 10, 16, 20):
        b2 = uSVP_beta(n, q, sigma)
        check(tag + 'independent uSVP estimate within 3%% of reported beta', abs(b2 - v['beta']) <= 0.03 * v['beta'], 'indep=%d reported=%d' % (b2, v['beta']))
    # ---------- 3. sizes ----------
    ct = ((2 * m + ELL) * k + L) / 8.0
    kn = 2 * m * ELL * math.ceil(math.log2(s) + 4) / 8.0
    pk = (n * m * (2 * L + 2) + n * ELL) * k / 8.0
    check(tag + 'sizes ct/node, key/node, pk', abs(ct - v['ct_node_B']) < 1 and abs(kn - v['key_node_B']) < 1 and abs(pk - v['pk_B']) < 1)
    check(tag + 'key/attribute = (h+1)*key/node', abs(v['key_attr_B'] - (h + 1) * kn) < 1)
    check(tag + 'linear baseline (N+2)*m*k/8', abs(v['bgg_ct_B'] - (N + 2) * m * k / 8.0) < 1)

# ---------- 4. calibration of the estimator ----------
b = uSVP_beta(512, 3329, math.sqrt(1.5))
check('estimator calibration (n=512,q=3329,cbd var 1.5): beta ~ 406, bits ~ 118.6', abs(b - 406) <= 12 and abs(0.292 * b - 118.6) <= 3.6, 'indep beta=%d bits=%.1f' % (b, 0.292 * b))
print('\nFAILED: %d' % len(fails) if fails else '\nALL CHECKS PASSED'); print(fails)
sys.exit(1 if fails else 0)
