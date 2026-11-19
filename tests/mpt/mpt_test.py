# Differential test of src/mpt.c against a trie built from scratch. Usage: mpt_test.py <lib.so> [rounds] [seed]
import ctypes, random, sys
L = ctypes.CDLL(sys.argv[1]); rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 200; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
L.airbender_mpt_open.restype = ctypes.c_void_p
for f in (L.airbender_mpt_get_raw, L.airbender_mpt_set_raw, L.airbender_mpt_root, L.airbender_mpt_clear): f.argtypes = None
def keccak(b):
    out = ctypes.create_string_buffer(32); L.t_keccak(b, len(b), out); return out.raw
assert keccak(b"").hex() == "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470"
EMPTY = keccak(b"\x80")
assert EMPTY.hex() == "56e81f171bcc55a6ff8345e692c0f86e5b48e01b996cadc001622fb5e363b421"
def enc_len(n, base):
    if n < 56: return bytes([base + n])
    b = n.to_bytes((n.bit_length() + 7) // 8, "big"); return bytes([base + 55 + len(b)]) + b
def enc_str(b): return b if len(b) == 1 and b[0] < 0x80 else enc_len(len(b), 0x80) + b
def enc_list(items): p = b"".join(items); return enc_len(len(p), 0xc0) + p
def hp(nibs, leaf):
    f = 2 if leaf else 0
    if len(nibs) % 2: nibs = [f + 1] + list(nibs)
    else: nibs = [f, 0] + list(nibs)
    return bytes(nibs[i] << 4 | nibs[i + 1] for i in range(0, len(nibs), 2))
def nibs_of(k): return [x for b in k for x in (b >> 4, b & 15)]
def build(items, depth, nodes, root=False):
    """items: sorted list of (nibbles, value). Returns the encoding of the node, b'' for none."""
    if not items: return b""
    if len(items) == 1:
        e = enc_list([enc_str(hp(items[0][0][depth:], True)), enc_str(items[0][1])])
    else:
        cp = 0; a, z = items[0][0], items[-1][0]
        while a[depth + cp] == z[depth + cp]: cp += 1
        if cp:
            e = enc_list([enc_str(hp(a[depth:depth + cp], False)), ref(build(items, depth + cp, nodes))])
        else:
            groups = [[] for _ in range(16)]
            for it in items: groups[it[0][depth]].append(it)
            e = enc_list([ref(build(g, depth + 1, nodes)) for g in groups] + [b"\x80"])
    if len(e) >= 32 or root: nodes[keccak(e)] = e
    return e
def ref(e): return b"\x80" if not e else e if len(e) < 32 else b"\xa0" + keccak(e)
def root_of(state, nodes):
    items = sorted((nibs_of(k), v) for k, v in state.items())
    e = build(items, 0, nodes, True)
    return keccak(e) if e else EMPTY
rnd = random.Random(seed)
def rkey(pool):
    m = rnd.random()
    if pool and m < 0.5:
        base = bytearray(rnd.choice(pool)); n = rnd.choice([1, 1, 2, 3, 8, 20, 31, 32])
        for i in range(32 - n, 32): base[i] = rnd.randrange(256)
        if rnd.random() < 0.3: base[31] = (base[31] & 0xf0) | rnd.randrange(16)
        return bytes(base)
    return bytes(rnd.randrange(256) for _ in range(32))
def rval():
    m = rnd.random()
    if m < 0.25: return bytes([rnd.randrange(1, 128)])
    if m < 0.5: return bytes(rnd.randrange(256) for _ in range(rnd.randrange(1, 6)))
    if m < 0.95: return bytes(rnd.randrange(256) for _ in range(rnd.randrange(1, 110)))
    return bytes(rnd.randrange(256) for _ in range(rnd.randrange(110, 400)))
def get(t, k):
    p = ctypes.c_void_p(); n = ctypes.c_uint32()
    e = L.airbender_mpt_get_raw(ctypes.c_void_p(t), k, ctypes.byref(p), ctypes.byref(n))
    if e == 1: return None
    assert e == 0, ("get error", e)
    return ctypes.string_at(p.value, n.value)
def check(t, state, extra, what):
    out = ctypes.create_string_buffer(32)
    e = L.airbender_mpt_root(ctypes.c_void_p(t), out); assert e == 0, (what, "root error", e)
    want = root_of(state, {})
    assert out.raw == want, (what, "root", out.raw.hex(), want.hex(), len(state))
    for k in list(state)[:60] + extra:
        assert get(t, k) == state.get(k), (what, "get", k.hex())
ops = 0
for r in range(rounds):
    state = {}; pool = []
    for _ in range(rnd.choice([0, 1, 2, 3, 5, 17, 40, 200, 600])):
        k = rkey(pool); pool.append(k); state[k] = rval()
    nodes = {}; root = root_of(state, nodes)
    L.t_reset()
    for e in nodes.values(): L.t_add_node(e, len(e))
    t = L.airbender_mpt_open(root)
    check(t, state, [rkey(pool) for _ in range(20)], "open")
    for batch in range(rnd.choice([1, 2, 4])):
        for _ in range(rnd.choice([1, 3, 10, 60, 300])):
            m = rnd.random(); ops += 1
            if m < 0.35 and state:
                k = rnd.choice(list(state)); v = rval(); state[k] = v
                assert L.airbender_mpt_set_raw(ctypes.c_void_p(t), k, v, len(v)) == 0
            elif m < 0.65 and state:
                k = rnd.choice(list(state)); del state[k]
                assert L.airbender_mpt_set_raw(ctypes.c_void_p(t), k, b"", 0) == 0
            elif m < 0.72:
                k = rkey(pool)
                if k in state: del state[k]
                assert L.airbender_mpt_set_raw(ctypes.c_void_p(t), k, b"", 0) == 0
            else:
                k = rkey(pool); pool.append(k); v = rval(); state[k] = v
                assert L.airbender_mpt_set_raw(ctypes.c_void_p(t), k, v, len(v)) == 0
            if rnd.random() < 0.05: check(t, state, [], "mid")
        check(t, state, [rkey(pool) for _ in range(20)], "batch")
    if rnd.random() < 0.2:
        for k in list(state):
            assert L.airbender_mpt_set_raw(ctypes.c_void_p(t), k, b"", 0) == 0
        state = {}; check(t, state, [], "emptied")
# The account calls: fields in, the encoding and the fields out.
L.airbender_mpt_open_storage.restype = ctypes.c_void_p
t = L.airbender_mpt_open(EMPTY); accounts = {}
for i in range(400):
    addr = bytes(rnd.randrange(256) for _ in range(20))
    nonce = rnd.choice([0, 1, 0x7f, 0x80, rnd.randrange(1 << 64), rnd.randrange(1 << 20)])
    bal = rnd.choice([0, 1, 0x7f, 0x80, rnd.randrange(1 << 256), rnd.randrange(1 << 70), (1 << 256) - 1])
    sr = rnd.choice([EMPTY, bytes(rnd.randrange(256) for _ in range(32))]); ch = rnd.choice([keccak(b""), bytes(rnd.randrange(256) for _ in range(32))])
    n2 = (ctypes.c_uint32 * 2)(nonce & 0xffffffff, nonce >> 32); b8 = (ctypes.c_uint32 * 8)(*[(bal >> (32 * j)) & 0xffffffff for j in range(8)])
    assert L.airbender_mpt_set_account_fields(ctypes.c_void_p(t), addr, n2, b8, sr, ch) == 0
    def num(x): return enc_str(x.to_bytes((x.bit_length() + 7) // 8, "big")) if x else b"\x80"
    accounts[addr] = (nonce, bal, sr, ch, enc_list([num(nonce), num(bal), enc_str(sr), enc_str(ch)]))
state = {keccak(a): v[4] for a, v in accounts.items()}
out = ctypes.create_string_buffer(32); assert L.airbender_mpt_root(ctypes.c_void_p(t), out) == 0 and out.raw == root_of(state, {})
for addr, (nonce, bal, sr, ch, enc) in accounts.items():
    f = (ctypes.c_uint32 * 27)(); p = ctypes.c_void_p(); n = ctypes.c_uint32()
    assert L.airbender_mpt_get_account_fields(ctypes.c_void_p(t), addr, f, ctypes.byref(p), ctypes.byref(n)) == 0
    assert ctypes.string_at(p.value, n.value) == enc
    assert f[0] | f[1] << 32 == nonce and sum(f[2 + j] << (32 * j) for j in range(8)) == bal
    assert bytes(f)[40:72] == sr and bytes(f)[72:104] == ch and f[26] == (sr == EMPTY) | (ch == keccak(b"")) << 1
    st = ctypes.c_int(); s2 = L.airbender_mpt_open_storage(ctypes.c_void_p(t), addr, ctypes.byref(st)); assert st.value == 0
    assert L.airbender_mpt_root(ctypes.c_void_p(s2), out) == 0 and out.raw == sr
f = (ctypes.c_uint32 * 27)(); p = ctypes.c_void_p(); n = ctypes.c_uint32()
assert L.airbender_mpt_get_account_fields(ctypes.c_void_p(t), bytes(20), f, ctypes.byref(p), ctypes.byref(n)) == 1
# A deletion from a two-child branch needs the other leaf to fold its path.
keys = [bytes([n]) + bytes(31) for n in (0x10, 0x20)]
nodes = {}
root = root_of({k: k for k in keys}, nodes)
L.t_reset()
for node in nodes.values():
    if not node.endswith(keys[1]):
        L.t_add_node(node, len(node))
t = ctypes.c_void_p(L.airbender_mpt_open(root))
assert get(t.value, keys[0]) == keys[0]
assert L.airbender_mpt_set_raw(t, keys[0], b"", 0) == -1
assert L.airbender_mpt_root(t, out) == 0 and out.raw == root
assert get(t.value, keys[0]) == keys[0]

addr = bytes(20)
t = ctypes.c_void_p(L.airbender_mpt_open(EMPTY))
hashes = [enc_str(EMPTY), enc_str(keccak(b""))]
for fields in ([b"\x81\x01", b"\x80"], [b"\x80", b"\x81\x7f"], [b"\x80"] * 3):
    account = enc_list(fields + hashes)
    assert L.airbender_mpt_set_account(t, addr, account, len(account)) == 0
    assert L.airbender_mpt_get_account_fields(t, addr, f, ctypes.byref(p), ctypes.byref(n)) == 2
    status = ctypes.c_int()
    assert not L.airbender_mpt_open_storage(t, addr, ctypes.byref(status)) and status.value == 2

slot = (ctypes.c_uint32 * 8)(1)
slot_key = keccak((1).to_bytes(32, "big"))
value = (ctypes.c_uint32 * 8)()
for number in (1, 127, 128, (1 << 256) - 1, 0):
    words = (ctypes.c_uint32 * 8)(*[(number >> (32 * j)) & 0xffffffff for j in range(8)])
    assert L.airbender_mpt_set_slot(t, slot, words) == 0
    assert L.airbender_mpt_get_slot(t, slot, value) == 0 and list(value) == list(words)
    assert get(t.value, slot_key) == (enc_str(number.to_bytes((number.bit_length() + 7) // 8, "big")) if number else None)
assert L.airbender_mpt_set_raw(t, slot_key, b"\x81\x7f", 2) == 0
assert L.airbender_mpt_get_slot(t, slot, value) == -2
print("ok: rounds", rounds, "ops", ops, "edited fast path calls", L.t_edited_calls())
