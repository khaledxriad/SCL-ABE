"""Complete-subtree 'maximal fully-authorized subtree' cover over a heap-indexed binary tree."""
import bisect, random, math, csv, json, sys

def cover(h, revoked):
    """Return list of heap node indices (root=1, leaves 2^h..2^(h+1)-1) whose subtrees
    are maximal subtrees containing no revoked leaf. Exact partition of [N]\\R."""
    R = sorted(revoked)
    out = []
    def rec(node, lo, hi):            # leaves [lo,hi)
        i = bisect.bisect_left(R, lo)
        if i == len(R) or R[i] >= hi:
            out.append(node); return
        if hi - lo == 1:
            return
        mid = (lo+hi)//2
        rec(2*node, lo, mid); rec(2*node+1, mid, hi)
    rec(1, 0, 1 << h)
    return out

def leaves_of(h, node):
    d = node.bit_length()-1
    span = 1 << (h-d)
    lo = (node - (1 << d)) * span
    return lo, lo+span

def path_nodes(h, leaf):
    v = (1 << h) + leaf
    out = []
    while v >= 1:
        out.append(v); v >>= 1
    return out            # leaf ... root, h+1 nodes

def check(h, revoked):
    cv = cover(h, revoked)
    N = 1 << h
    seen = [0]*N
    for v in cv:
        lo, hi = leaves_of(h, v)
        for x in range(lo, hi): seen[x] += 1
    R = set(revoked)
    for x in range(N):
        assert seen[x] == (0 if x in R else 1), (x, seen[x])
    # key-path property: a leaf u is in cover subtree v iff v on path(u)
    return len(cv)

def bound(N, r):
    return r*(math.log2(N/r)+1) if r >= 1 else 1

if __name__ == "__main__":
    random.seed(7)
    # exhaustive correctness test, h<=6 (all subsets for h<=4, random for larger)
    for h in range(1, 5):
        N = 1 << h
        for mask in range(1 << N):
            check(h, [i for i in range(N) if mask >> i & 1])
    for h in range(5, 9):
        N = 1 << h
        for _ in range(300):
            r = random.randint(0, N)
            check(h, random.sample(range(N), r))
    print("cover correctness tests passed")
    rows = []
    for h in [10, 12, 14, 16, 18, 20]:
        N = 1 << h
        for r in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384]:
            if r > N//2: continue
            trials = 300 if N <= (1<<16) else 100
            rand_sizes, blk_sizes = [], []
            for _ in range(trials):
                R = random.sample(range(N), r)
                rand_sizes.append(len(cover(h, R)))
                s0 = random.randrange(0, N-r+1)
                blk_sizes.append(len(cover(h, range(s0, s0+r))))
            rows.append(dict(h=h, N=N, r=r,
                avg_rand=sum(rand_sizes)/trials, max_rand=max(rand_sizes), min_rand=min(rand_sizes),
                avg_block=sum(blk_sizes)/trials, max_block=max(blk_sizes),
                bound=bound(N, r), worst_complete=min(N-r, bound(N, r))))
    with open('../data/cover_stats.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    bad = [x for x in rows if x['max_rand'] > x['bound'] + 1e-9 or x['max_block'] > x['bound'] + 1e-9]
    print("rows:", len(rows), "bound violations:", len(bad))
    for x in rows:
        if x['h'] in (10, 16) and x['r'] in (1, 8, 64, 512):
            print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in x.items()})
