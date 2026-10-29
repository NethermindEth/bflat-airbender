"""Checks the accelerators against pure-Python references. Usage: run.py <harness-command...>"""
import hashlib
import os
import random
import subprocess
import sys

rng = random.Random(1234)


# --- Keccak-256 reference -------------------------------------------------------------------
def keccak_f(a):
    rc = [0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
          0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
          0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
          0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
          0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
          0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008]
    rot = [[0, 36, 3, 41, 18], [1, 44, 10, 45, 2], [62, 6, 43, 15, 61], [28, 55, 25, 21, 56], [27, 20, 39, 8, 14]]
    m = (1 << 64) - 1
    rl = lambda x, n: ((x << n) | (x >> (64 - n))) & m if n else x
    for r in range(24):
        c = [a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20] for x in range(5)]
        d = [c[(x - 1) % 5] ^ rl(c[(x + 1) % 5], 1) for x in range(5)]
        a = [a[i] ^ d[i % 5] for i in range(25)]
        b = [0] * 25
        for x in range(5):
            for y in range(5):
                b[y + 5 * ((2 * x + 3 * y) % 5)] = rl(a[x + 5 * y], rot[x][y])
        a = [b[i] ^ (~b[(i % 5 + 1) % 5 + 5 * (i // 5)] & b[(i % 5 + 2) % 5 + 5 * (i // 5)]) for i in range(25)]
        a[0] ^= rc[r]
    return a


def keccak256(data):
    rate = 136
    data = bytearray(data) + b"\x01"
    data += b"\x00" * (-len(data) % rate)
    data[-1] |= 0x80
    st = [0] * 25
    for off in range(0, len(data), rate):
        for i in range(rate // 8):
            st[i] ^= int.from_bytes(data[off + 8 * i: off + 8 * i + 8], "little")
        st = keccak_f(st)
    return b"".join(x.to_bytes(8, "little") for x in st[:4])


# --- curves -----------------------------------------------------------------------------------
class Curve:
    def __init__(self, p, n, a, b, gx, gy):
        self.p, self.n, self.a, self.b, self.g = p, n, a, b, (gx, gy)

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0:
                return None
            lam = (3 * P[0] * P[0] + self.a) * pow(2 * P[1], -1, p) % p
        else:
            lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
        x = (lam * lam - P[0] - Q[0]) % p
        return x, (lam * (P[0] - x) - P[1]) % p

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R


K1 = Curve(2**256 - 2**32 - 977, 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141, 0, 7,
           0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
           0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
R1 = Curve(0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
           0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
           0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC,
           0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
           0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
           0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5)


def sign(cv, d, h):
    while True:
        k = rng.randrange(1, cv.n)
        R = cv.mul(k, cv.g)
        r = R[0] % cv.n
        s = pow(k, -1, cv.n) * (h + r * d) % cv.n
        if r and s:
            recid = (R[1] & 1) | (2 if R[0] >= cv.n else 0)
            return r, s, recid


h32 = lambda v: v.to_bytes(32, "big").hex()

cases = []  # (line, expected)

for data in [b"", b"abc"] + [rng.randbytes(rng.randrange(0, 400)) for _ in range(40)] + [rng.randbytes(n) for n in (55, 56, 63, 64, 64, 65, 119, 120, 128)] + [bytes(64), bytes(136), bytes(272)]:
    arg = data.hex() or "-"
    cases.append((f"keccak {arg}", keccak256(data).hex()))
    cases.append((f"sha256 {arg}", hashlib.sha256(data).hexdigest()))
try:
    for data in [b"", b"abc"] + [rng.randbytes(rng.randrange(0, 300)) for _ in range(20)]:
        cases.append((f"ripemd {data.hex() or '-'}", "00" * 12 + hashlib.new("ripemd160", data).hexdigest()))
except ValueError:
    cases.append(("ripemd -", "00" * 12 + "9c1185a5c5e9fc54612808977ee8f548b2258d31"))
    cases.append(("ripemd 616263", "00" * 12 + "8eb208f7e05d987a9b044a8e98c6b087f15a0bfc"))

# blake2b("abc") is one F call with the final flag, rounds 12, over the parameter-xored IV.
iv = [0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
      0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179]
for msg in [b"abc", b"", rng.randbytes(100)]:
    h = iv[:]
    h[0] ^= 0x01010040
    inp = (12).to_bytes(4, "big") + b"".join(x.to_bytes(8, "little") for x in h)
    inp += msg.ljust(128, b"\x00") + len(msg).to_bytes(8, "little") + bytes(8) + b"\x01"
    cases.append((f"blake2f {inp.hex()}", hashlib.blake2b(msg).hexdigest()))

# Random keys and hashes, then the edge keys and hashes: the GLV split and the wNAF tables see
# every shape of scalar, including u1 = 0 (hash 0 or n).
k1_cases = [(rng.randrange(1, K1.n), rng.randrange(0, 2**256)) for _ in range(200)]
k1_cases += [(d, rng.randrange(0, 2**256)) for d in (1, 2, 3, K1.n - 1, K1.n - 2)]
k1_cases += [(rng.randrange(1, K1.n), h) for h in (0, K1.n, 2**256 - 1)]
for d, h in k1_cases:
    r, s, recid = sign(K1, d, h % K1.n)
    Q = K1.mul(d, K1.g)
    cases.append((f"ecrecover {h32(h)} {h32(r)}{h32(s)} {recid:02x}", h32(Q[0]) + h32(Q[1])))
    # The high-s twin recovers the same key with the other parity.
    cases.append((f"ecrecover {h32(h)} {h32(r)}{h32(K1.n - s)} {recid ^ 1:02x}", h32(Q[0]) + h32(Q[1])))
cases.append((f"ecrecover {h32(1)} {h32(0)}{h32(1)} 00", "fail"))
cases.append((f"ecrecover {h32(1)} {h32(1)}{h32(K1.n)} 00", "fail"))
cases.append((f"ecrecover {h32(1)} {h32(1)}{h32(1)} 04", "fail"))
# x = 5 is not on secp256k1 (5^3 + 7 = 132 is a non-residue mod p).
cases.append((f"ecrecover {h32(1)} {h32(5)}{h32(1)} 00",
              "fail" if pow(132, (K1.p - 1) // 2, K1.p) != 1 else None))


# The same recoveries in batches, the failing signatures among them: each item gives what its own call gives.
# 600 items cross the native side's 512-signature chunk.
def ecrecover_batch(items):
    msgs = "".join(m for m, _, _, _ in items)
    sigs = "".join(sg for _, sg, _, _ in items)
    recids = "".join(rc for _, _, rc, _ in items)
    return (f"ecrecoverbatch {msgs} {sigs} {recids}", ",".join(e for _, _, _, e in items))


k1_items = [tuple(line.split()[1:4]) + (expected,) for line, expected in cases if line.startswith("ecrecover ") and expected]
bad = [item for item in k1_items if item[3] == "fail"]
good = [item for item in k1_items if item[3] != "fail"]
cases.append(ecrecover_batch(good[:1]))
cases.append(ecrecover_batch(bad[:1] + good[:1]))
cases.append(ecrecover_batch(bad[:2]))
mixed = []
for i, item in enumerate(good):
    mixed.append(item)
    if i % 50 == 7:
        mixed.append(bad[i % len(bad)])
cases.append(ecrecover_batch(mixed))
cases.append(ecrecover_batch((mixed * 2)[:600]))


# Batch verification of the recovered keys: R's y from r and the recovery id, as a host gives it. Any key, y,
# signature or hash out of place fails the whole batch; 600 items cross the 512-signature chunk.
def r_y(sig, recid):
    x = int(sig[:64], 16) + (K1.n if recid & 2 else 0)
    y = pow((x ** 3 + 7) % K1.p, (K1.p + 1) // 4, K1.p)
    return y if (y & 1) == (recid & 1) else K1.p - y


def verify_batch(items, keys=None, ys=None):
    msgs = "".join(m for m, _, _, _ in items)
    sigs = "".join(sg for _, sg, _, _ in items)
    recids = "".join(rc for _, _, rc, _ in items)
    keys = keys or [e for _, _, _, e in items]
    ys = ys or [h32(r_y(sg, int(rc, 16))) for _, sg, rc, _ in items]
    weights = rng.randbytes(16 * len(items)).hex()
    return f"verifybatch {msgs} {sigs} {recids} {''.join(keys)} {''.join(ys)} {weights}"


for items in (good[:1], good[:2], good[:40], good, (good * 3)[:600]):
    cases.append((verify_batch(items), "1"))
    keys = [e for _, _, _, e in items]
    ys = [h32(r_y(sg, int(rc, 16))) for _, sg, rc, _ in items]
    i = rng.randrange(len(items))
    # A key that differs (the high-s twins share theirs), the key negated, then R's y negated: each a point on the curve.
    other = next((k for k in keys if k != keys[i]), keys[i][:64] + h32(K1.p - int(keys[i][64:], 16)))
    cases.append((verify_batch(items, keys[:i] + [other] + keys[i + 1:], ys), "0"))
    cases.append((verify_batch(items, keys[:i] + [keys[i][:64] + h32(K1.p - int(keys[i][64:], 16))] + keys[i + 1:], ys), "0"))
    cases.append((verify_batch(items, keys, ys[:i] + [h32(K1.p - int(ys[i], 16))] + ys[i + 1:]), "0"))
    tampered = list(items)
    m, sg, rc, e = tampered[i]
    tampered[i] = (h32(int(m, 16) ^ 1), sg, rc, e)
    cases.append((verify_batch(tampered, keys, ys), "0"))
cases.append((verify_batch(good[:3] + bad[:1], [e for _, _, _, e in good[:3]] + [good[0][3]],
                           [h32(r_y(sg, int(rc, 16))) for _, sg, rc, _ in good[:3]] + [h32(1)]), "0"))

for _ in range(6):
    d = rng.randrange(1, R1.n)
    h = rng.randrange(0, 2**256)
    r, s, _ = sign(R1, d, h % R1.n)
    Q = R1.mul(d, R1.g)
    pk = h32(Q[0]) + h32(Q[1])
    cases.append((f"p256 {h32(h)} {h32(r)}{h32(s)} {pk}", "1"))
    cases.append((f"p256 {h32(h ^ 1)} {h32(r)}{h32(s)} {pk}", "0"))
    cases.append((f"p256 {h32(h)} {h32(r)}{h32(s)} {h32(Q[0])}{h32(Q[1] ^ 1)}", "0"))

for _ in range(60):
    bl, el, ml = rng.randrange(0, 70), rng.randrange(0, 40), rng.randrange(0, 70)
    b, e, m = rng.randbytes(bl), rng.randbytes(el), rng.randbytes(ml)
    if ml and rng.random() < 0.3:
        m = b"\x00" * (ml // 2) + m[ml // 2:]
    mi = int.from_bytes(m, "big")
    exp = "-" if ml == 0 else (pow(int.from_bytes(b, "big"), int.from_bytes(e, "big"), mi) if mi else 0).to_bytes(ml, "big").hex()
    cases.append((f"modexp {b.hex() or '-'} {e.hex() or '-'} {m.hex() or '-'}", exp))
for ml in (128, 256, 512):
    b, e, m = rng.randbytes(ml), rng.randbytes(8), rng.randbytes(ml)
    exp = pow(int.from_bytes(b, "big"), int.from_bytes(e, "big"), int.from_bytes(m, "big")).to_bytes(ml, "big").hex()
    cases.append((f"modexp {b.hex()} {e.hex()} {m.hex()}", exp))
# Odd moduli of up to 1024 bits take the Montgomery path: every length, exponents with leading zeros and of
# one bit, bases above the modulus.
for mbytes in list(range(1, 33)) + [32, 32, 33, 48, 63, 64, 64, 65, 96, 127, 128, 128]:
    m = rng.randrange(1, 256 ** mbytes) | 1
    if m == 1:
        m = 3
    for b, e in [(rng.randrange(0, 256 ** 40), rng.randrange(0, 256 ** 32)), (rng.randrange(0, m), 1),
                 (rng.randrange(0, m), 0), (m + 5, rng.randrange(0, 256 ** 3))]:
        bb = b.to_bytes(max(1, (b.bit_length() + 7) // 8), "big")
        eb = e.to_bytes(32, "big")
        mb = m.to_bytes(mbytes, "big")
        cases.append((f"modexp {bb.hex()} {eb.hex()} {mb.hex()}", pow(b, e, m).to_bytes(mbytes, "big").hex()))

# BN254, BLS12-381 and KZG: precompile vectors and py_ecc cases built by gen_vectors.py.
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vectors", "ec.txt")) as f:
    for line in f:
        *args, want = line.split()
        cases.append((" ".join(args), want))

cases = [c for c in cases if c[1] is not None]
out = subprocess.run(sys.argv[1:], input="\n".join(c[0] for c in cases) + "\n", capture_output=True, text=True).stdout.split("\n")
bad = 0
for (line, want), got in zip(cases, out):
    if got != want:
        bad += 1
        if bad <= 10:
            print(f"MISMATCH {line.split()[0]}: got {got[:80]} want {want[:80]}")
print(f"{len(cases) - bad}/{len(cases)} passed")
sys.exit(1 if bad else 0)
