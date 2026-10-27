"""Prints the constants of src/bn254.c and src/bls12.c: the Montgomery-form field and curve
tables, and the BLS12-381 generators and [tau]G2 of the KZG trusted setup. Needs py_ecc
(pip install py_ecc) for the RFC 9380 isogeny-map coefficients and point decompression.

Usage: gen_ec_constants.py <kzg_trusted_setup.txt>"""
import sys
from py_ecc.optimized_bls12_381 import constants as K
from py_ecc import optimized_bls12_381 as B
from py_ecc.bls.point_compression import decompress_G1, decompress_G2


def limbs(v, n=12):
    return [(v >> (32 * i)) & 0xFFFFFFFF for i in range(n)]


def fp_c(v, n):
    ls = limbs(v, n)
    return "{{" + ", ".join(f"0x{x:08x}" for x in ls) + "}}"


class Ctx:
    def __init__(self, p, n):
        self.p, self.n = p, n
        self.R = 1 << (32 * n)

    def mont(self, v):
        return v * self.R % self.p

    def fp(self, v):
        return fp_c(self.mont(v % self.p), self.n)

    def fp2(self, a):
        return "{{" + self.fp(a[0]) + ", " + self.fp(a[1]) + "}}"

    # Fp2 = Fp[u]/(u^2 + 1)
    def mul2(self, a, b):
        p = self.p
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def pow2(self, a, e):
        r = (1, 0)
        while e:
            if e & 1:
                r = self.mul2(r, a)
            a = self.mul2(a, a)
            e >>= 1
        return r

    def inv2(self, a):
        p = self.p
        t = pow(a[0] * a[0] + a[1] * a[1], -1, p)
        return (a[0] * t % p, -a[1] * t % p)


def field(name, ctx):
    p, n = ctx.p, ctx.n
    inv = (-pow(p, -1, 1 << 32)) % (1 << 32)
    print(f"/* {name} */")
    print(f"    .n = {n},")
    print(f"    .inv = 0x{inv:08x},")
    print(f"    .p = {fp_c(p, n)},")
    print(f"    .one = {ctx.fp(1)},")
    print(f"    .r2 = {fp_c(ctx.R * ctx.R % p, n)},")


# ---------------------------------------------------------------- BN254
x = 4965661367192848881
p = 36 * x**4 + 36 * x**3 + 24 * x**2 + 6 * x + 1
r = 36 * x**4 + 36 * x**3 + 18 * x**2 + 6 * x + 1
c = Ctx(p, 8)
xi = (9, 1)
field("BN254", c)
print(f"    .b1 = {c.fp2((3, 0))},")
print(f"    .b2 = {c.fp2(c.mul2((3, 0), c.inv2(xi)))},")
print(f"    .gamma = {c.fp2(c.pow2(xi, (p - 1) // 6))},")
print(f"    loop 6x+2 = 0x{6 * x + 2:x}, r = 0x{r:064x}")

# ---------------------------------------------------------------- BLS12-381
x = 0xD201000000010000
p = B.field_modulus
r = B.curve_order
c = Ctx(p, 12)
xi = (1, 1)
field("BLS12-381", c)
print(f"    .b1 = {c.fp2((4, 0))},")
print(f"    .b2 = {c.fp2((4, 4))},")
print(f"    .gamma = {c.fp2(c.pow2(xi, (p - 1) // 6))},")
print(f"    r = 0x{r:064x}")

# The cube root of unity beta for which phi(P) = (beta x, y) acts on G1 as -x^2.
g = 2
while True:
    w = pow(g, (p - 1) // 3, p)
    if w != 1:
        break
    g += 1
G = B.normalize(B.G1)
for beta in (w, w * w % p):
    phi = (B.FQ(beta * G[0].n % p), G[1])
    m = B.normalize(B.neg(B.multiply(B.G1, x * x)))
    if (phi[0], phi[1]) == (m[0], m[1]):
        print(f"beta = {c.fp(beta)}")


def fq_list(name, coeffs, deg):
    items = []
    for k in coeffs:
        if deg == 1:
            items.append(c.fp(k.n))
        else:
            items.append(c.fp(k.coeffs[0].n if hasattr(k.coeffs[0], "n") else k.coeffs[0]))
            items.append(c.fp(k.coeffs[1].n if hasattr(k.coeffs[1], "n") else k.coeffs[1]))
    print(f"/* {name}: {len(coeffs)} coefficients */")
    print("    " + ",\n    ".join(items) + ",")


print(f"g1 A = {c.fp(K.ISO_11_A.n)}")
print(f"g1 B = {c.fp(K.ISO_11_B.n)}")
print(f"g1 Z = {c.fp(11)}")
for name, k in zip(["xnum", "xden", "ynum", "yden"], K.ISO_11_MAP_COEFFICIENTS):
    fq_list("g1 " + name, list(k), 1)


def c2(v):
    return tuple(int(t.n) if hasattr(t, "n") else int(t) for t in v.coeffs)


print(f"g2 A = {c.fp2(c2(K.ISO_3_A))}")
print(f"g2 B = {c.fp2(c2(K.ISO_3_B))}")
print(f"g2 Z = {c.fp2(c2(K.ISO_3_Z))}")
for name, k in zip(["xnum", "xden", "ynum", "yden"], K.ISO_3_MAP_COEFFICIENTS):
    k = list(k)
    while c2(k[-1]) == (0, 0):
        k.pop()
    fq_list("g2 " + name, k, 2)


# ---------------------------------------------------------------- KZG setup points
def arr(b, ind="    "):
    s = [f"0x{x:02x}" for x in b]
    return "\n".join(ind + ", ".join(s[i:i + 16]) + "," for i in range(0, len(s), 16))


def be(v, n=48):
    return int(v).to_bytes(n, "big")


def g2bytes(P):
    return be(P[0].coeffs[0]) + be(P[0].coeffs[1]) + be(P[1].coeffs[0]) + be(P[1].coeffs[1])


lines = open(sys.argv[1]).read().split()
n1, n2 = int(lines[0]), int(lines[1])
g2m = [decompress_G2((int(s[:96], 16), int(s[96:], 16))) for s in lines[2 + n1: 2 + n1 + 2]]
tau1 = decompress_G1(int(lines[2 + n1 + n2 + 1], 16))
assert B.normalize(g2m[0]) == B.normalize(B.G2)
# [tau]G1 of the same file pairs consistently with [tau]G2.
assert B.final_exponentiate(B.pairing(B.G2, tau1)) == B.final_exponentiate(B.pairing(g2m[1], B.G1))
g1 = B.normalize(B.G1)
print("G1Gen")
print(arr(be(g1[0]) + be(g1[1])))
print("G2Gen")
print(arr(g2bytes(B.normalize(B.G2))))
print("TauG2")
print(arr(g2bytes(B.normalize(g2m[1]))))
