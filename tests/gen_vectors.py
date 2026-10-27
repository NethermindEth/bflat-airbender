"""Builds tests/vectors/ec.txt, the BN254, BLS12-381 and KZG cases of run.py, as lines
"<op> <hex>... <expected>".

Sources: the precompile vectors Nethermind carries (EIP-2537 JSONs from the EIP and geth, the
bn256 JSONs of the precompile benchmarks, the c-kzg point-evaluation vectors) and cases generated
here with py_ecc (pip install py_ecc), the Ethereum Foundation's reference implementation of
both curves and of RFC 9380 hash-to-curve. Inputs are reshaped into the accelerator layout:
EIP-2537 field elements lose their 16 zero bytes, and inputs the precompile rejects before
calling the accelerator (bad padding, bad length, a versioned-hash mismatch) are left out.

Usage: gen_vectors.py <nethermind repo root> <kzg_trusted_setup.txt>
"""
import glob
import hashlib
import json
import os
import random
import sys

from py_ecc import optimized_bls12_381 as B
from py_ecc import optimized_bn128 as N
from py_ecc.bls.hash_to_curve import clear_cofactor_G1, clear_cofactor_G2, map_to_curve_G1, map_to_curve_G2
from py_ecc.bls.point_compression import compress_G1, decompress_G1
from py_ecc.fields import optimized_bls12_381_FQ as BFQ, optimized_bls12_381_FQ2 as BFQ2

rng = random.Random(4844)
root, setup_path = sys.argv[1], sys.argv[2]
out = []


def add(op, args, want):
    out.append(f"{op} {' '.join(a.hex() if a else '-' for a in args)} {want}")


def h(b):
    return b.hex()


# --- vectors from the repository -----------------------------------------------------------------
def unhex(s):
    return bytes.fromhex(s[2:] if s.startswith("0x") else s)


def trim(data, fields):
    """EIP-2537 padded field elements to the accelerator's 48-byte ones; None on bad padding."""
    res = b""
    for i in range(fields):
        fe = data[64 * i: 64 * i + 64]
        if any(fe[:16]):
            return None
        res += fe[16:]
    return res


bls_dir = os.path.join(root, "src/Nethermind/Nethermind.Evm.Test/PrecompileVectors/Bls")
shapes = {  # file stem -> (op, input item size, fields per item, scalar bytes per item, output fields)
    "add_G1_bls": ("blsg1add", 256, 4, 0, 2),
    "add_G2_bls": ("blsg2add", 512, 8, 0, 4),
    "mul_G1_bls": ("blsg1msm", 160, 2, 32, 2),
    "mul_G2_bls": ("blsg2msm", 288, 4, 32, 4),
    "multiexp_G1_bls": ("blsg1msm", 160, 2, 32, 2),
    "multiexp_G2_bls": ("blsg2msm", 288, 4, 32, 4),
    "pairing_check_bls": ("blspair", 384, 6, 0, None),
    "map_fp_to_G1_bls": ("blsmapg1", 64, 1, 0, 2),
    "map_fp2_to_G2_bls": ("blsmapg2", 128, 2, 0, 4),
}
repo_counts = {}
for path in sorted(glob.glob(os.path.join(bls_dir, "*.json"))):
    stem = os.path.basename(path)[:-5]
    fail = stem.startswith("fail-")
    op, size, fields, scalar, ofields = shapes[stem.removeprefix("fail-")]
    for v in json.load(open(path)):
        data = unhex(v["Input"])
        if not data or len(data) % size or (op in ("blsg1add", "blsg2add", "blsmapg1", "blsmapg2") and len(data) != size):
            continue
        acc = b""
        for off in range(0, len(data), size):
            item = data[off: off + size]
            t = trim(item[: size - scalar], fields)
            if t is None:
                break
            acc += t + item[size - scalar:]
        else:
            if fail:
                want = "fail"
            elif op == "blspair":
                want = str(unhex(v["Expected"])[-1])
            else:
                want = trim(unhex(v["Expected"]), ofields).hex()
            add(op, [acc], want)
            repo_counts[stem] = repo_counts.get(stem, 0) + 1

bench = os.path.join(root, "src/bench_precompiles/vectors")
for name, op, size in (("bnadd", "bnadd", 128), ("bnmul", "bnmul", 96), ("bnpair", "bnpair", None)):
    for v in json.load(open(os.path.join(bench, name, f"common_{name}.json"))):
        data = unhex(v["Input"])
        if size:
            data = data[:size].ljust(size, b"\0")
        want = unhex(v["Expected"])
        add(op, [data], str(want[-1]) if op == "bnpair" else want.hex())
        repo_counts[name] = repo_counts.get(name, 0) + 1

BLS_R = B.curve_order
BLS_R_BYTES = BLS_R.to_bytes(32, "big")


def add_kzg(data, valid, source):
    if len(data) != 192:
        return
    vh, z, y, c, proof = data[:32], data[32:64], data[64:96], data[96:144], data[144:]
    if b"" + hashlib.sha256(c).digest()[1:] != vh:
        return
    add("kzg", [c, z, y, proof], "1" if valid else "0")
    repo_counts[source] = repo_counts.get(source, 0) + 1


# The c-kzg verify_kzg_proof vectors (mainnet setup) as EEST carries them; their names tell the
# outcome. The benchmark's own success vector does not verify under the mainnet setup (py_ecc
# agrees), so only its failure vectors are taken.
eest = os.path.join(root, "src/tests/LegacyTests/Cancun/GeneralStateTests/Pyspecs/cancun/eip4844_blobs/"
                          "point_evaluation_precompile_external_vectors.json")
for name, t in json.load(open(eest)).items():
    for data in t["transaction"]["data"]:
        add_kzg(unhex(data), "case_correct_proof" in name, "eest_kzg_external")
for v in json.load(open(os.path.join(bench, "point_evaluation/current/fail-pointEvaluation.json"))):
    add_kzg(unhex(v["Input"]), False, "fail-pointEvaluation")

# --- BN254 cases generated with py_ecc -------------------------------------------------------------
NP = N.field_modulus
NR = N.curve_order


def bn_g1(P):
    if N.is_inf(P):
        return bytes(64)
    x, y = N.normalize(P)
    return x.n.to_bytes(32, "big") + y.n.to_bytes(32, "big")


def bn_g2(Q):
    if N.is_inf(Q):
        return bytes(128)
    x, y = N.normalize(Q)
    be = lambda v: int(v).to_bytes(32, "big")
    return be(x.coeffs[1]) + be(x.coeffs[0]) + be(y.coeffs[1]) + be(y.coeffs[0])


def bn_rand_g1():
    return N.multiply(N.G1, rng.randrange(1, NR))


for _ in range(10):
    P, Q = bn_rand_g1(), bn_rand_g1()
    add("bnadd", [bn_g1(P) + bn_g1(Q)], h(bn_g1(N.add(P, Q))))
    add("bnadd", [bn_g1(P) + bn_g1(P)], h(bn_g1(N.double(P))))
    add("bnadd", [bn_g1(P) + bn_g1(N.neg(P))], h(bytes(64)))
    add("bnadd", [bytes(64) + bn_g1(Q)], h(bn_g1(Q)))
    k = rng.randrange(0, 2**256)
    add("bnmul", [bn_g1(P) + k.to_bytes(32, "big")], h(bn_g1(N.multiply(P, k % NR))))
for k in (0, 1, 2, NR - 1, NR, NR + 1, 2**256 - 1):
    add("bnmul", [bn_g1(N.G1) + k.to_bytes(32, "big")], h(bn_g1(N.multiply(N.G1, k % NR))))
add("bnmul", [bytes(64) + (5).to_bytes(32, "big")], h(bytes(64)))
G = bn_g1(N.G1)
bad_y = G[:32] + (int.from_bytes(G[32:], "big") ^ 1).to_bytes(32, "big")
add("bnadd", [bad_y + G], "fail")
add("bnadd", [G + NP.to_bytes(32, "big") + bytes(32)], "fail")
add("bnadd", [G + (1).to_bytes(32, "big") + (2 + NP).to_bytes(32, "big")], "fail")
add("bnmul", [bad_y + (1).to_bytes(32, "big")], "fail")

for n in (2, 3):
    pairs, acc = [], 0
    for i in range(n - 1):
        a, b = rng.randrange(1, NR), rng.randrange(1, NR)
        pairs.append(bn_g1(N.multiply(N.G1, a)) + bn_g2(N.multiply(N.G2, b)))
        acc += a * b
    pairs.append(bn_g1(N.neg(N.multiply(N.G1, acc % NR))) + bn_g2(N.G2))
    good = b"".join(pairs)
    add("bnpair", [good], "1")
    add("bnpair", [good[:-192] + bn_g1(N.double(N.G1)) + good[-128:]], "0")
add("bnpair", [bn_g1(N.G1) + bn_g2(N.G2)], "0")
add("bnpair", [bytes(64) + bn_g2(N.G2)], "1")
add("bnpair", [bn_g1(N.G1) + bytes(128)], "1")
add("bnpair", [], "1")
# A twist point outside G2: the first x = k + u with x^3 + b' a square.
FQ2 = N.FQ2
twist_pts = []
for xi in range(1, 400):
    X = FQ2([xi, 1])
    rhs = X**3 + N.b2
    # rhs is a square in Fp2 iff its norm is a square in Fp
    c0, c1 = int(rhs.coeffs[0]), int(rhs.coeffs[1])
    if pow((c0 * c0 + c1 * c1) % NP, (NP - 1) // 2, NP) != 1:
        continue
    # sqrt in Fp2 (p = 3 mod 4): Algorithm 9 of Adj and Rodriguez-Henriquez
    a1 = rhs ** ((NP - 3) // 4)
    alpha = a1 * a1 * rhs
    x0 = a1 * rhs
    Y = FQ2([0, 1]) * x0 if alpha == FQ2([-1, 0]) else (alpha + FQ2.one()) ** ((NP - 1) // 2) * x0
    assert Y * Y == rhs
    P = (X, Y, FQ2.one())
    assert not N.is_inf(N.multiply(P, NR))
    twist_pts.append(P)
    if len(twist_pts) == 4:
        break
for P in twist_pts:
    add("bnpair", [bn_g1(N.G1) + bn_g2(P)], "fail")
    add("bnpair", [bytes(64) + bn_g2(P)], "fail")
add("bnpair", [bad_y + bn_g2(N.G2)], "fail")
g2b = bytearray(bn_g2(N.G2))
g2b[127] ^= 1
add("bnpair", [bn_g1(N.G1) + bytes(g2b)], "fail")
add("bnpair", [bn_g1(N.G1) + NP.to_bytes(32, "big") + bn_g2(N.G2)[32:]], "fail")

# --- BLS12-381 cases generated with py_ecc -----------------------------------------------------------
BP = B.field_modulus


def fe(v):
    return int(v).to_bytes(48, "big")


def g1(P):
    if B.is_inf(P):
        return bytes(96)
    x, y = B.normalize(P)
    return fe(x.n) + fe(y.n)


def g2(Q):
    if B.is_inf(Q):
        return bytes(192)
    x, y = B.normalize(Q)
    return fe(x.coeffs[0]) + fe(x.coeffs[1]) + fe(y.coeffs[0]) + fe(y.coeffs[1])


def rand_g1():
    return B.multiply(B.G1, rng.randrange(1, BLS_R))


def rand_g2():
    return B.multiply(B.G2, rng.randrange(1, BLS_R))


def off_g1():
    """On the curve, outside G1: the SWU map before cofactor clearing."""
    return map_to_curve_G1(BFQ(rng.randrange(BP)))


def off_g2():
    return map_to_curve_G2(BFQ2([rng.randrange(BP), rng.randrange(BP)]))


for _ in range(4):
    for enc, rnd, off, op, double, neg in ((g1, rand_g1, off_g1, "blsg1add", B.double, B.neg),
                                           (g2, rand_g2, off_g2, "blsg2add", B.double, B.neg)):
        P, Q, O = rnd(), rnd(), off()
        assert not B.is_inf(B.multiply(O, BLS_R))
        add(op, [enc(P) + enc(Q)], h(enc(B.add(P, Q))))
        add(op, [enc(P) + enc(P)], h(enc(double(P))))
        add(op, [enc(P) + enc(neg(P))], h(enc(B.Z1 if enc is g1 else B.Z2)))
        add(op, [enc(O) + enc(Q)], h(enc(B.add(O, Q))))  # no subgroup check on addition
        add(op, [enc(O) + enc(O)], h(enc(double(O))))
        add(op, [bytes(len(enc(P))) + enc(O)], h(enc(O)))
bad = bytearray(g1(B.G1))
bad[-1] ^= 1
add("blsg1add", [bytes(bad) + g1(B.G1)], "fail")
add("blsg1add", [fe(BP) + fe(0) + g1(B.G1)], "fail")
bad2 = bytearray(g2(B.G2))
bad2[-1] ^= 1
add("blsg2add", [bytes(bad2) + g2(B.G2)], "fail")

for op, enc, rnd, off, mul, add_ in (("blsg1msm", g1, rand_g1, off_g1, B.multiply, B.add),
                                     ("blsg2msm", g2, rand_g2, off_g2, B.multiply, B.add)):
    zero = B.Z1 if enc is g1 else B.Z2
    for n in (1, 1, 2, 3, 5):
        pts = [rnd() for _ in range(n)]
        ks = [rng.choice([0, 1, BLS_R - 1, BLS_R, 2**256 - 1, rng.randrange(2**256), rng.randrange(2**256)])
              for _ in range(n)]
        acc = zero
        for P, k in zip(pts, ks):
            acc = add_(acc, mul(P, k % BLS_R))
        add(op, [b"".join(enc(P) + k.to_bytes(32, "big") for P, k in zip(pts, ks))], h(enc(acc)))
    P = rnd()
    add(op, [enc(zero) + (7).to_bytes(32, "big") + enc(P) + (3).to_bytes(32, "big")], h(enc(mul(P, 3))))
    add(op, [enc(P) + bytes(32)], h(enc(zero)))
    add(op, [enc(off()) + (1).to_bytes(32, "big")], "fail")
    add(op, [enc(P) + (1).to_bytes(32, "big") + enc(off()) + bytes(32)], "fail")

for n in (2, 3):
    pairs, acc = [], 0
    for i in range(n - 1):
        a, b = rng.randrange(1, BLS_R), rng.randrange(1, BLS_R)
        pairs.append(g1(B.multiply(B.G1, a)) + g2(B.multiply(B.G2, b)))
        acc += a * b
    pairs.append(g1(B.neg(B.multiply(B.G1, acc % BLS_R))) + g2(B.G2))
    good = b"".join(pairs)
    add("blspair", [good], "1")
    add("blspair", [good[:-288] + g1(B.double(B.G1)) + good[-192:]], "0")
add("blspair", [g1(B.G1) + g2(B.G2)], "0")
add("blspair", [g1(B.G1) + g2(B.Z2)], "1")
add("blspair", [], "1")
add("blspair", [g1(off_g1()) + g2(B.G2)], "fail")
add("blspair", [g1(B.G1) + g2(off_g2())], "fail")
add("blspair", [g1(B.Z1) + g2(off_g2())], "fail")

for _ in range(8):
    u = rng.randrange(BP)
    add("blsmapg1", [fe(u)], h(g1(clear_cofactor_G1(map_to_curve_G1(BFQ(u))))))
    u0, u1 = rng.randrange(BP), rng.randrange(BP)
    add("blsmapg2", [fe(u0) + fe(u1)], h(g2(clear_cofactor_G2(map_to_curve_G2(BFQ2([u0, u1]))))))
for u in (0, 1, BP - 1):
    add("blsmapg1", [fe(u)], h(g1(clear_cofactor_G1(map_to_curve_G1(BFQ(u))))))
    add("blsmapg2", [fe(u) + fe(0)], h(g2(clear_cofactor_G2(map_to_curve_G2(BFQ2([u, 0]))))))
    add("blsmapg2", [fe(0) + fe(u)], h(g2(clear_cofactor_G2(map_to_curve_G2(BFQ2([0, u]))))))
add("blsmapg1", [fe(BP)], "fail")
add("blsmapg2", [fe(0) + fe(BP)], "fail")

# --- KZG proofs over the mainnet setup ------------------------------------------------------------
lines = open(setup_path).read().split()
n1, n2 = int(lines[0]), int(lines[1])
monomial = [decompress_G1(int(s, 16)) for s in lines[2 + n1 + n2: 2 + n1 + n2 + 5]]


def commit(coeffs):
    acc = B.Z1
    for c, P in zip(coeffs, monomial):
        acc = B.add(acc, B.multiply(P, c % BLS_R))
    return compress_G1(acc).to_bytes(48, "big")


for _ in range(4):
    coeffs = [rng.randrange(BLS_R) for _ in range(4)]
    z = rng.randrange(BLS_R)
    y = sum(c * pow(z, i, BLS_R) for i, c in enumerate(coeffs)) % BLS_R
    q2 = coeffs[3]  # (p(X) - y) / (X - z) by synthetic division
    q1 = (coeffs[2] + z * q2) % BLS_R
    q0 = (coeffs[1] + z * q1) % BLS_R
    q = [q0, q1, q2]
    C, P = commit(coeffs), commit(q)
    zb, yb = z.to_bytes(32, "big"), y.to_bytes(32, "big")
    add("kzg", [C, zb, yb, P], "1")
    add("kzg", [C, zb, ((y + 1) % BLS_R).to_bytes(32, "big"), P], "0")
    add("kzg", [C, ((z + 1) % BLS_R).to_bytes(32, "big"), yb, P], "0")
    add("kzg", [C, (z + BLS_R).to_bytes(32, "big") if z + BLS_R < 2**256 else BLS_R_BYTES, yb, P], "0")
    add("kzg", [C, zb, BLS_R_BYTES, P], "0")
    add("kzg", [bytes([C[0] & 0x7F]) + C[1:], zb, yb, P], "0")
    add("kzg", [bytes([C[0] ^ 0x20]) + C[1:], zb, yb, P], "0")
inf = bytes([0xC0]) + bytes(47)
add("kzg", [inf, (5).to_bytes(32, "big"), bytes(32), inf], "1")
add("kzg", [inf, (5).to_bytes(32, "big"), (1).to_bytes(32, "big"), inf], "0")
add("kzg", [inf, BLS_R_BYTES, bytes(32), inf], "0")
add("kzg", [bytes([0xC0]) + bytes(46) + b"\x01", bytes(32), bytes(32), inf], "0")
add("kzg", [bytes([0xE0]) + bytes(47), bytes(32), bytes(32), inf], "0")
off = compress_G1(off_g1()).to_bytes(48, "big")
add("kzg", [off, bytes(32), bytes(32), inf], "0")
xp = bytearray(BP.to_bytes(48, "big"))
xp[0] |= 0x80
add("kzg", [bytes(xp), bytes(32), bytes(32), inf], "0")

os.makedirs(os.path.join(os.path.dirname(__file__), "vectors"), exist_ok=True)
with open(os.path.join(os.path.dirname(__file__), "vectors", "ec.txt"), "w", newline="\n") as f:
    f.write("\n".join(out) + "\n")
print(f"{len(out)} cases; from the repository: {repo_counts}", file=sys.stderr)
