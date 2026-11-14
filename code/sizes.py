import math, json, csv, random
from est import choose, derived
from cover import cover
random.seed(11)
T_BITS, KAPPA = 16, 128
ELL = KAPPA // T_BITS

def bits_per_key_entry(s): return math.ceil(math.log2(s) + 4.0)

out = {}
for h in [6, 8, 10, 12, 14, 16, 20]:
    N = 1 << h
    d = choose(N, lam=128, t_bits=T_BITS)
    L, n, k, m, s = d['L'], d['n'], d['k'], d['m'], d['s']
    ct_node = ((2*m + ELL)*k + L)/8.0
    kb = bits_per_key_entry(s)
    key_node = 2*m*ELL*kb/8.0
    pk = (n*m*(2*L+2) + n*ELL)*k/8.0
    msk = d['mbar']*n*k*5/8.0
    bgg_ct = (N+2)*m*k/8.0            # (N+1) attribute encodings + c0, all m-vectors over Z_q (same lattice parameters)
    # measured average cover sizes by Monte-Carlo (random revocation) for this N
    cs = {}
    for r in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        if r > N//2: continue
        tr = 100 if N <= (1 << 14) else 30
        cs[r] = sum(len(cover(h, random.sample(range(N), r))) for _ in range(tr))/tr
    out[N] = dict(h=h, L=L, n=n, k=k, m=m, mbar=d['mbar'], s=s, alpha_log2=math.log2(d['alpha']), aq=d['aq'], beta=d['beta'],
                  bits=d['bits'], loss_bits=d['loss_bits'], sigma_e=d['sigma_e'],
                  ct_node_B=ct_node, key_node_B=key_node, key_attr_B=(h+1)*key_node, pk_B=pk, msk_B=msk,
                  bgg_ct_B=bgg_ct, key_entry_bits=kb, cover_avg=cs)
    print(N, {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in out[N].items() if kk != 'cover_avg'}, flush=True)
json.dump(out, open('../data/params.json', 'w'), indent=1)
