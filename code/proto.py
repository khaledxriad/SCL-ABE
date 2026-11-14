"""Toy-parameter prototype of the subset-cover lattice ABE (SC-LABE): dual-Regev IBE with MP12 gadget trapdoor
(base 2, q = 2^62) + complete-subtree cover.  Gaussians are rounded continuous Gaussians: benchmarking only."""
import numpy as np, math, time, hashlib, json, sys
from cover import cover, path_nodes
K = 62; Q = 1 << K; MASK = np.uint64(Q - 1)
ETA = 3.0
U64 = np.uint64

def u(x):   # signed int64 -> uint64 two's complement
    return np.asarray(x, dtype=np.int64).view(np.uint64)

def mm(a, b): return (a @ b) & MASK     # uint64 wrap arithmetic, modulus 2^62

def gauss(rng, std, size):
    return np.rint(rng.normal(0.0, std, size)).astype(np.int64)

class IBE:
    def __init__(self, n, L, t_bits=16, kappa=128, seed=1):
        self.rng = np.random.default_rng(seed)
        self.n, self.L, self.t = n, L, t_bits
        self.ell = kappa // t_bits
        self.k = K
        self.mbar = n*K + 256; self.nw = n*K; self.m = self.mbar + self.nw
        m, rng = self.m, self.rng
        tau = math.sqrt(2*math.log(2*2.0**40))
        sigR = ETA
        self.sigG = math.sqrt(5)*ETA
        t0 = time.perf_counter()
        # --- TrapGen (MP12): A = [Abar | G - Abar R]
        self.Abar = rng.integers(0, Q, size=(n, self.mbar), dtype=np.uint64)
        self.R = gauss(rng, sigR, (self.mbar, self.nw))
        G = np.zeros((n, self.nw), dtype=np.uint64)
        for i in range(n):
            G[i, i*K:(i+1)*K] = np.array([1 << j for j in range(K)], dtype=np.uint64)
        AR = mm(self.Abar, u(self.R))
        self.A = np.concatenate([self.Abar, (G - AR) & MASK], axis=1)       # n x m
        # --- public matrices
        self.B0 = rng.integers(0, Q, size=(n, m), dtype=np.uint64)
        self.Bj = rng.integers(0, Q, size=(L, 2, n, m), dtype=np.uint64)
        self.U = rng.integers(0, Q, size=(n, self.ell), dtype=np.uint64)
        # --- preprocessing for SamplePre: Sigma_p = s^2 I - sigG^2 M M^T
        s1 = np.linalg.norm(self.R.astype(np.float64), 2)
        self.s1 = s1
        s_mp = math.sqrt(self.sigG**2*(s1**2+1) + ETA**2)
        s_sr = math.sqrt(5)*ETA*(L+1)*12.0*math.sqrt(2*m)
        self.s = 1.02*max(s_mp, s_sr)
        M = np.concatenate([self.R.astype(np.float64), np.eye(self.nw)], axis=0)    # m x nw
        S = -(self.sigG**2) * (M @ M.T)
        S[np.diag_indices(m)] += self.s**2
        self.Lp = np.linalg.cholesky(S)
        del S, M
        self.t_setup = time.perf_counter() - t0
        # noise parameter (same derivation as est.derived)
        coef = tau*math.sqrt(L+1)*self.s*m + tau*self.s*math.sqrt(m) + tau
        self.alpha = 1.0/(2*(2**t_bits)*coef)
        self.sigma_e = self.alpha*Q/math.sqrt(2*math.pi)
        self.coef = coef

    # ---- helpers
    def Bid(self, idx):
        acc = self.B0.copy()
        for j in range(self.L):
            acc = (acc + self.Bj[j, (idx >> j) & 1]) & MASK
        return acc

    def sample_g(self, v):          # v: int64 array (n, c) in [0,q) -> z int64 (nw, c)
        n, c = v.shape; rng = self.rng
        z = np.zeros((n, K, c), dtype=np.int64)
        cur = v.copy()
        for j in range(K):
            b = cur & 1
            g = rng.normal(0.0, self.sigG, size=(n, c))
            zj = (2*np.rint((g - b)/2.0)).astype(np.int64) + b
            z[:, j, :] = zj
            cur = (cur - zj) >> 1
        return z.reshape(n*K, c)

    def keygen(self, idx):
        """SampleLeft: returns E in Z^{2m x ell} with [A|B_id] E = U (mod q)."""
        m, ell, rng = self.m, self.ell, self.rng
        B = self.Bid(idx)
        e2 = gauss(rng, self.s, (m, ell))
        Y = (self.U - mm(B, u(e2))) & MASK
        p = np.rint(self.Lp @ rng.normal(size=(m, ell))).astype(np.int64)
        v = ((Y - mm(self.A, u(p))) & MASK).astype(np.int64)
        z = self.sample_g(v)
        x_top = p[:self.mbar] + self.R @ z
        x_bot = p[self.mbar:] + z
        e1 = np.concatenate([x_top, x_bot], axis=0)
        E = np.concatenate([e1, e2], axis=0)
        return E

    def check_key(self, idx, E):
        F = np.concatenate([self.A, self.Bid(idx)], axis=1)
        return bool(np.all(mm(F, u(E)) == self.U))

    def enc(self, idx, mu):
        """mu: int array (ell,) in Z_p -> (c0 (ell,), c1 (2m,)) """
        n, m, ell, rng = self.n, self.m, self.ell, self.rng
        s = rng.integers(0, Q, size=n, dtype=np.uint64)
        x = gauss(rng, self.sigma_e, ell)
        y = gauss(rng, self.sigma_e, m)
        # z = R'^T y with R' = sum of (L+1) fresh +-1 matrices (entries Binomial)
        z = np.zeros(m, dtype=np.float64); yf = y.astype(np.float64)
        ch = 512
        for r0 in range(0, m, ch):
            r1 = min(m, r0+ch)
            Rc = (rng.binomial(self.L+1, 0.5, size=(r1-r0, m))*2 - (self.L+1)).astype(np.float64)
            z += Rc.T @ yf[r0:r1]
        z = np.rint(z).astype(np.int64)
        F = np.concatenate([self.A, self.Bid(idx)], axis=1)
        c1 = (F.T @ s) & MASK
        c1 = (c1 + u(np.concatenate([y, z]))) & MASK
        delta = U64(Q >> self.t)
        c0 = ((self.U.T @ s) + u(x) + u(np.asarray(mu, dtype=np.int64)) * delta) & MASK
        return c0, c1

    def dec(self, E, c0, c1):
        w = (c0 - (E.T.astype(np.int64).view(np.uint64) @ c1)) & MASK
        half = U64(1 << (K - self.t - 1))
        return (((w + half) & MASK) >> U64(K - self.t)).astype(np.int64)

    def noise_margin(self, E, c0, c1, mu):
        w = (c0 - (u(E.T) @ c1)) & MASK
        err = (w.astype(np.int64) - (np.asarray(mu, dtype=np.int64) << (K-self.t))) 
        err = ((err + (Q//2)) % Q) - Q//2
        return float(np.max(np.abs(err))) / (Q >> (self.t+1))

def dem(key_bytes, msg):
    out = bytearray()
    ctr = 0
    while len(out) < len(msg):
        out += hashlib.sha256(key_bytes + ctr.to_bytes(4, 'big')).digest(); ctr += 1
    return bytes(a ^ b for a, b in zip(msg, out))

class ABE:
    def __init__(self, h, n, seed=1):
        self.h = h; self.N = 1 << h
        self.ibe = IBE(n, h+1, seed=seed)         # node ids are heap indices in [1, 2^(h+1)-1]: L = h+1 bits
    def keygen(self, a):
        return {v: self.ibe.keygen(v) for v in path_nodes(self.h, a)}
    def encrypt(self, revoked, msg):
        ib = self.ibe
        mu = ib.rng.integers(0, 1 << ib.t, size=ib.ell)
        kb = b"".join(int(x).to_bytes(2, 'big') for x in mu)
        cts = {v: ib.enc(v, mu) for v in cover(self.h, revoked)}
        return cts, dem(kb, msg)
    def decrypt(self, sk, a, ct):
        cts, body = ct
        for v in path_nodes(self.h, a):
            if v in cts:
                mu = self.ibe.dec(sk[v], *cts[v])
                kb = b"".join(int(x).to_bytes(2, 'big') for x in mu)
                return dem(kb, body)
        return None

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    t0 = time.time()
    abe = ABE(4, n)
    ib = abe.ibe
    print("n=%d m=%d setup %.1fs  s=%.0f s1R=%.0f alpha=%.3g sigma_e=%.3g" % (n, ib.m, ib.t_setup, ib.s, ib.s1, ib.alpha, ib.sigma_e))
    sk = {a: abe.keygen(a) for a in [0, 3, 5, 9]}
    print("keys ok:", all(ib.check_key(v, E) for a in sk for v, E in sk[a].items()))
    msg = b"hello sub-linear lattice ABE"
    ct = abe.encrypt([3, 9], msg)
    print("cover size", len(ct[0]), {a: abe.decrypt(sk[a], a, ct) for a in sk})
    # noise margin
    v = 1; E = ib.keygen(v); mu = ib.rng.integers(0, 1 << ib.t, size=ib.ell)
    c0, c1 = ib.enc(v, mu)
    print("dec ok", np.array_equal(ib.dec(E, c0, c1), mu), "noise/threshold = %.3g" % ib.noise_margin(E, c0, c1, mu))
    print("total %.1fs" % (time.time()-t0))
