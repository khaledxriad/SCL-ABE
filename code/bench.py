import numpy as np, time, json, math, statistics as st
from proto import *
from cover import cover
res = []
h = 6; N = 1 << h
for n in [8, 16, 24, 32]:
    abe = ABE(h, n, seed=n); ib = abe.ibe
    row = dict(n=n, m=ib.m, L=ib.L, ell=ib.ell, N=N, setup_s=ib.t_setup, alpha=ib.alpha, s=ib.s)
    def timeit(f, reps):
        ts = []
        for _ in range(reps):
            t0 = time.perf_counter(); out = f(); ts.append(time.perf_counter()-t0)
        return st.mean(ts), (st.pstdev(ts) if len(ts) > 1 else 0.0), out
    kg_m, kg_s, E = timeit(lambda: ib.keygen(5), 5)
    row['keygen_node_s'] = kg_m; row['keygen_node_sd'] = kg_s
    mu = ib.rng.integers(0, 1 << ib.t, size=ib.ell)
    en_m, en_s, (c0, c1) = timeit(lambda: ib.enc(5, mu), 3)
    row['enc_node_s'] = en_m; row['enc_node_sd'] = en_s
    de_m, de_s, dm = timeit(lambda: ib.dec(E, c0, c1), 50)
    row['dec_node_s'] = de_m; row['dec_node_sd'] = de_s
    row['dec_ok'] = bool(np.array_equal(dm, mu))
    row['noise_ratio'] = ib.noise_margin(E, c0, c1, mu)
    # end-to-end encryption vs cover size
    e2e = []
    rng = np.random.default_rng(99)
    for r in [1, 2, 4, 8, 16]:
        times = []; cs = []
        for _ in range(2):
            R = list(rng.choice(N, size=r, replace=False))
            t0 = time.perf_counter(); ct = abe.encrypt(R, b"x"*64); times.append(time.perf_counter()-t0); cs.append(len(ct[0]))
        e2e.append(dict(r=r, cover=st.mean(cs), enc_s=st.mean(times)))
    row['e2e'] = e2e
    # keygen for a full attribute (h+1 nodes)
    t0 = time.perf_counter(); sk = abe.keygen(7); row['keygen_attr_s'] = time.perf_counter()-t0
    # linear (BGG+-style) baseline encryption cost model: (N+1) matvecs F_i^T s + noise
    Fs = [ib.rng.integers(0, Q, size=(n, ib.m), dtype=np.uint64) for _ in range(N+1)]
    def lin():
        s = ib.rng.integers(0, Q, size=n, dtype=np.uint64)
        return [((F.T @ s) + u(gauss(ib.rng, ib.sigma_e, ib.m))) & MASK for F in Fs]
    lm, ls, _ = timeit(lin, 3)
    row['linear_enc_s'] = lm; row['linear_N'] = N
    del Fs
    # sizes
    row['ct_node_bytes'] = (2*ib.m + ib.ell)*K/8
    row['key_node_bytes'] = 2*ib.m*ib.ell*8        # int64 storage in prototype
    res.append(row)
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != 'e2e'}), flush=True)
    json.dump(res, open('../data/bench.json', 'w'), indent=1)
print("done")
