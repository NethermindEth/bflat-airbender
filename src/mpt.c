/*
 * The state and storage tries of the guest's world state, over the witness nodes of the trusted input.
 *
 * A trie starts as the hash of its root. An operation builds each node it first reaches from the witness encoding
 * its hash finds, once: later walks follow pointers. A branch built from an encoding lists only the children an
 * operation reached or changed; the rest stay where they lie in the encoding. A write changes nodes as a Merkle
 * Patricia trie does: a leaf splits into a branch, a branch left with one child folds into it. The root is hashed
 * last, each changed node once: a branch whose children changed but kept their shape is its witness encoding with
 * the new hashes written over the old, and a full one resumes the sponge from the state kept after its last
 * unchanged block.
 *
 * Every witness node is checked when an operation first reaches it, as SecuredTrieNodeValidator of Nethermind
 * does after the zkEVM spec's incremental_mpt.py: a malformed node off the block's path does not fail the block.
 *
 * Included by a platform file, which defines before it:
 *   uint32_t mpt_w_count(void), const uint32_t *mpt_w_hash(uint32_t k), const uint8_t *mpt_w_rlp(uint32_t k,
 *   uint32_t *len): the witness nodes, their digests as eight aligned words;
 *   void mpt_keccak(const uint8_t *data, uint32_t len, uint32_t out[8]);
 *   void mpt_hash_edited(uint32_t k, int block, const uint8_t *original, const uint8_t *edited, uint32_t out[8]):
 *   the digest of full branch k as the 532 bytes at `edited`, which are valid from rate block `block` on and
 *   equal `original` before it;
 *   void mpt_copy(uint8_t *to, const uint8_t *from, uint32_t len), which may copy up to 31 bytes more on either
 *   side, for `to` and `from` equal modulo 32;
 *   void mpt_keccak_key(const uint32_t *words, uint32_t len, uint32_t out[8]): the keccak of the 20 or 32 bytes
 *   in aligned words;
 *   MPT_UNALIGNED_WORDS where a word may be read as the aligned words that cover it.
 */

/* Its own, though the platform file has them: read alone, as a code scanner reads it, the file knows its types. */
#include <stddef.h>
#include <stdint.h>

enum
{
    MptOk = 0,
    MptAbsent = 1,
    MptMissing = -1,
    MptInvalid = -2,
    MptTooLong = -3
};

enum
{
    TEmpty = 0,
    THash = 1,
    TEmbed = 2,
    TNode = 3
};

enum
{
    KBranch = 1,
    KExt = 2,
    KLeaf = 3
};

enum
{
    MptKeyNibbles = 64,
    MptFullBranch = 532,
    MptBufferBytes = 2048
};

/* A child: nothing, the 32 bytes of a hash, an embedded encoding with its length above the tag, or a built node. */
typedef struct
{
    uintptr_t p;
    uint32_t t;
} mref;

typedef struct mnode
{
    union
    {
        struct
        {
            uint8_t kind;
            uint8_t dirty;   /* its encoding is not computed */
            uint8_t changed; /* it is not the witness node it was built from */
            uint8_t hashed;  /* its parent refers to it by hash; else by the elen bytes at eraw */
        };
        uint32_t state;
    };
    union
    {
        struct
        {
            uint8_t plen, poff;
            uint8_t full; /* src is a full branch of hashes */
            uint8_t pad;
        };
        uint32_t shape;
    };
    uint32_t mask;  /* a branch's children that are in child[]; the others lie in src */
    uint32_t dmask; /* those of them a write changed */
    uint32_t elen;
    uint32_t slen;        /* the length of src */
    const uint8_t *src;   /* the encoding it was built from */
    const uint8_t *eraw;  /* its encoding, where a parent embeds it */
    const uint32_t *hash; /* its digest, where a parent refers to it by hash */
    const uint8_t *pk;    /* the path: nibble i is nibble poff + i of the bytes at pk */
    const uint8_t *val;
    uint32_t vlen;
    uint32_t widx; /* the witness node a full branch was built from, for its kept states */
    mref child[16]; /* a branch's; an extension has child[0], a leaf none */
} mnode;

typedef struct
{
    mref root;
    uint32_t rooth[8];
} mtrie;

static const uint32_t MptEmptyRoot[8] = {0x171fe856, 0xa655cc1b, 0xe64583ff, 0x6ef8c092,
                                         0x1be0485b, 0xc0ad6c99, 0xb52f6201, 0x21b463e3};

static int MptErr;
static int MptChanged;

/* ---- memory -------------------------------------------------------------------------------------- */

static uint8_t *MptArenaAt, *MptArenaEnd;

static void *mpt_alloc(uint32_t size)
{
    size = (size + 7) & ~7u;
    if ((uint32_t)(MptArenaEnd - MptArenaAt) < size)
    {
        uint32_t chunk = size > (1u << 22) ? size : (1u << 22);
        MptArenaAt = malloc(chunk);
        MptArenaEnd = MptArenaAt + chunk;
    }
    void *p = MptArenaAt;
    MptArenaAt += size;
    return p;
}

static const uint8_t *mpt_keep(const uint8_t *data, uint32_t len)
{
    uint8_t *p = mpt_alloc(len ? len : 1);
    if (((uintptr_t)data & 3) == 0)
    {
        /* Whole words: the arena rounds every piece up to eight bytes, and the source is read to the end of its
           last word. */
        const uint32_t *from = (const uint32_t *)data;
        uint32_t *to = (uint32_t *)p;
        for (uint32_t i = 0; i < (len + 3) / 4; i++)
            to[i] = from[i];
    }
    else
    {
        memcpy(p, data, len);
    }
    return p;
}

/* ---- bytes at any alignment ---------------------------------------------------------------------- */

static inline uint32_t mpt_nib(const uint8_t *b, uint32_t pos)
{
    uint32_t x = b[pos >> 1];
    return (pos & 1) ? (x & 15) : (x >> 4);
}

/* The 32 bytes at p as eight words. */
static inline __attribute__((always_inline)) void mpt_load32(uint32_t h[8], const uint8_t *p)
{
#ifdef MPT_UNALIGNED_WORDS
    uint32_t m = (uintptr_t)p & 3;
    const uint32_t *w = (const uint32_t *)(p - m);
    if (m == 0)
    {
#pragma clang loop unroll(full)
        for (int i = 0; i < 8; i++)
            h[i] = w[i];
        return;
    }
    uint32_t lo = 8 * m, hi = 32 - lo, cur = w[0];
#pragma clang loop unroll(full)
    for (int i = 0; i < 8; i++)
    {
        uint32_t next = w[i + 1];
        h[i] = (cur >> lo) | (next << hi);
        cur = next;
    }
#else
    memcpy(h, p, 32);
#endif
}

/* Eight words to the 32 bytes at p; the bytes around them stay. */
static inline __attribute__((always_inline)) void mpt_store32(uint8_t *p, const uint32_t h[8])
{
#ifdef MPT_UNALIGNED_WORDS
    uint32_t m = (uintptr_t)p & 3;
    uint32_t *w = (uint32_t *)(p - m);
    if (m == 0)
    {
#pragma clang loop unroll(full)
        for (int i = 0; i < 8; i++)
            w[i] = h[i];
        return;
    }
    uint32_t lo = 8 * m, hi = 32 - lo, mask = (1u << lo) - 1;
    w[0] = (w[0] & mask) | (h[0] << lo);
#pragma clang loop unroll(full)
    for (int i = 1; i < 8; i++)
        w[i] = (h[i - 1] >> hi) | (h[i] << lo);
    w[8] = (w[8] & ~mask) | (h[7] >> hi);
#else
    memcpy(p, h, 32);
#endif
}

/* ---- the witness nodes by hash ------------------------------------------------------------------- */

static uint32_t *MptSlots;
static uint32_t MptSlotMask;
/* 0 until the node is checked; then 1 for a branch, or 1 + where the items of a leaf or an extension start. */
static uint8_t *MptInfo;

/* The witness node the next first use is expected to be of: a writer that lists the nodes in the order a block
   first reaches them saves the guest the table below. */
static uint32_t MptNext;

/* A node was not the expected one: k it was, with `next` expected. For whoever orders a witness by tracing. */
__attribute__((noinline)) void mpt_node_out_of_order(uint32_t k, uint32_t next)
{
    __asm__ volatile("" : : "r"(k), "r"(next));
}

static void mpt_info(void)
{
    uint32_t n = mpt_w_count();
    MptInfo = malloc(n ? n : 1);
    memset(MptInfo, 0, n ? n : 1);
}

static void mpt_index(void)
{
    uint32_t n = mpt_w_count(), cap = 16;
    while (cap < 2 * n)
        cap <<= 1;
    MptSlots = malloc((size_t)cap * 4);
    memset(MptSlots, 0, (size_t)cap * 4);
    MptSlotMask = cap - 1;
    for (uint32_t k = 0; k < n; k++)
    {
        uint32_t i = mpt_w_hash(k)[0] & MptSlotMask;
        while (MptSlots[i])
            i = (i + 1) & MptSlotMask;
        MptSlots[i] = k + 1;
    }
}

static __attribute__((noinline)) int mpt_lookup_table(const uint32_t h[8])
{
    if (!MptSlots)
        mpt_index();
    uint32_t mask = MptSlotMask, i = h[0] & mask;
    const uint32_t *slots = MptSlots;
    for (;;)
    {
        uint32_t s = slots[i];
        if (!s)
            return -1;
        const uint32_t *c = mpt_w_hash(s - 1);
        if (c[0] == h[0]
            && ((c[1] ^ h[1]) | (c[2] ^ h[2]) | (c[3] ^ h[3]) | (c[4] ^ h[4]) | (c[5] ^ h[5]) | (c[6] ^ h[6])
                | (c[7] ^ h[7])) == 0)
            return (int)(s - 1);
        i = (i + 1) & mask;
    }
}

/* The witness node with digest h, or -1: the expected one, or else whichever the table finds. `use` moves the
   expectation on. */
static inline __attribute__((always_inline)) int mpt_lookup_words(const uint32_t h[8], int use)
{
    if (__builtin_expect(!MptInfo, 0))
        mpt_info();
    uint32_t next = MptNext;
    if (next < mpt_w_count())
    {
        const uint32_t *c = mpt_w_hash(next);
        if (((c[0] ^ h[0]) | (c[1] ^ h[1]) | (c[2] ^ h[2]) | (c[3] ^ h[3]) | (c[4] ^ h[4]) | (c[5] ^ h[5])
             | (c[6] ^ h[6]) | (c[7] ^ h[7])) == 0)
        {
            if (use)
                MptNext = next + 1;
            return (int)next;
        }
    }
    int k = mpt_lookup_table(h);
    if (use && k >= 0)
        mpt_node_out_of_order((uint32_t)k, next);
    return k;
}

/* ---- the check of a witness node ----------------------------------------------------------------- */

static int mpt_item_end(const uint8_t *n, int pos, int end);

static int mpt_are_items(const uint8_t *n, int pos, int end)
{
    while (pos < end)
    {
        pos = mpt_item_end(n, pos, end);
        if (pos < 0)
            return 0;
    }
    return 1;
}

/* Where the canonical RLP item at pos, below end, ends, or -1. */
static int mpt_item_end(const uint8_t *n, int pos, int end)
{
    int prefix = n[pos];
    if (prefix < 0x80)
        return pos + 1;
    int is_list = prefix >= 0xc0;
    int sl = prefix & 0x3f;
    int cs = pos + 1, ce;
    if (sl < 56)
    {
        ce = cs + sl;
        if (ce > end || (sl == 1 && !is_list && n[cs] < 0x80))
            return -1;
    }
    else
    {
        int ll = sl - 55;
        if (ll > 2 || cs + ll > end)
            return -1;
        int length = n[cs];
        if (length == 0)
            return -1;
        if (ll == 2)
            length = length << 8 | n[cs + 1];
        cs += ll;
        ce = cs + length;
        if (length < 56 || ce > end)
            return -1;
    }
    return is_list && !mpt_are_items(n, cs, ce) ? -1 : ce;
}

static inline int mpt_is_embedded(int prefix)
{
    return (unsigned)(prefix - 0xc0) < 31;
}

static inline int mpt_skip_hash_or_empty(const uint8_t *n, int *pos, int length, int *children)
{
    if (*pos >= length)
        return 0;
    int prefix = n[*pos];
    if (prefix == 0x80)
    {
        (*pos)++;
        return 1;
    }
    if (prefix != 0xa0)
        return 0;
    (*children)++;
    *pos += 33;
    return 1;
}

/* -1 unless the node may sit `remaining` nibbles above the end of a key; then 0 for a branch, and for a leaf or an
   extension where its items start. hashed is whether a hash refers to it below the root. */
static int mpt_validate(const uint8_t *node, int length, int remaining, int hashed)
{
    if (length == MptFullBranch)
    {
        uint32_t bad = (node[0] ^ 0xf9) | (node[1] ^ 0x02) | (node[2] ^ 0x11) | (node[MptFullBranch - 1] ^ 0x80);
#pragma clang loop unroll(full)
        for (int i = 0; i < 16; i++)
            bad |= node[3 + 33 * i] ^ 0xa0;
        if (bad == 0)
            return remaining <= 0 ? -1 : 0;
    }
    if (length < (hashed ? 32 : 1))
        return -1;
    int header = node[0];
    int items_start = header <= 0xf7 ? 1 : header - 0xf6;
    int list_end;
    if (header < 0xc0)
        list_end = -1;
    else if (header <= 0xf7)
        list_end = 1 + header - 0xc0;
    else if (header == 0xf8 && length > 1 && node[1] >= 56)
        list_end = 2 + node[1];
    else if (header == 0xf9 && length > 2 && node[1] != 0)
        list_end = 3 + (node[1] << 8 | node[2]);
    else
        list_end = -1;
    if (list_end != length || items_start >= length)
        return -1;
    int pos = items_start, children = 0, slot = 0;
    while (slot < 16 && mpt_skip_hash_or_empty(node, &pos, length, &children))
        slot++;
    if (slot == 16 && pos + 1 == length && node[pos] == 0x80)
        return remaining <= 0 || children < 2 ? -1 : 0;
    /* Not a branch of hashes: a branch with an embedded child, a leaf or an extension. */
    int mismatch = pos;
    int may_be_branch = mismatch < length && mpt_is_embedded(node[mismatch]);
    pos = items_start;
    children = 0;
    slot = 0;
    while (may_be_branch && slot < 16 && pos < length)
    {
        if (!mpt_skip_hash_or_empty(node, &pos, length, &children))
        {
            int prefix = node[pos];
            pos = mpt_is_embedded(prefix) ? mpt_item_end(node, pos, length) : -1;
            if (pos < 0)
                break;
            children++;
        }
        slot++;
    }
    if (slot == 16 && pos + 1 == length && node[pos] == 0x80)
        return remaining <= 0 || children < 2 ? -1 : 0;
    int pp = node[items_start];
    if (pp == 0x80 || pp > 0xb7)
        return -1;
    int second = mpt_item_end(node, items_start, length);
    if (second < 0 || second >= length || node[second] >= 0xc0 + 31 || mpt_item_end(node, second, length) != length)
        return -1;
    int path_bytes = second - items_start - (pp < 0x80 ? 0 : 1);
    int flags = node[second - path_bytes];
    int nibbles = 2 * path_bytes - 2 + ((flags >> 4) & 1);
    int child = node[second];
    int ok = (flags & 0x20) ? (nibbles == remaining && child < 0xc0)
                            : (nibbles > 0 && nibbles < remaining && (child == 0xa0 || mpt_is_embedded(child)));
    return ok ? items_start : -1;
}

/* ---- references ---------------------------------------------------------------------------------- */

static inline mref mpt_ref_hash(const uint8_t *p)
{
    mref r = {(uintptr_t)p, THash};
    return r;
}

static inline mref mpt_ref_node(mnode *n)
{
    mref r = {(uintptr_t)n, TNode};
    return r;
}

/* The reference an item of a branch or the second item of an extension holds. */
static inline mref mpt_ref_item(const uint8_t *item)
{
    mref r;
    uint32_t p = item[0];
    if (p == 0x80)
    {
        r.p = 0;
        r.t = TEmpty;
    }
    else if (p == 0xa0)
    {
        r.p = (uintptr_t)(item + 1);
        r.t = THash;
    }
    else
    {
        r.p = (uintptr_t)item;
        r.t = TEmbed | (1 + p - 0xc0) << 8;
    }
    return r;
}

static inline uint32_t mpt_item_size(uint32_t prefix)
{
    return prefix == 0x80 ? 1 : prefix == 0xa0 ? 33 : 1 + prefix - 0xc0;
}

static inline const uint8_t *mpt_items(const uint8_t *rlp)
{
    return rlp + (rlp[0] <= 0xf7 ? 1 : rlp[0] - 0xf6);
}

/* Child i of a branch. */
static inline __attribute__((always_inline)) mref mpt_child(const mnode *n, uint32_t i)
{
    if ((n->mask >> i) & 1)
        return n->child[i];
    if (n->full)
        return mpt_ref_hash(n->src + 4 + 33 * i);
    const uint8_t *at = mpt_items(n->src);
    for (uint32_t j = 0; j < i; j++)
        at += mpt_item_size(at[0]);
    return mpt_ref_item(at);
}

/* The sixteen children of a branch. */
static void mpt_children(const mnode *n, mref out[16])
{
    uint32_t mask = n->mask;
    if (mask == 0xffff)
    {
        for (int i = 0; i < 16; i++)
            out[i] = n->child[i];
        return;
    }
    const uint8_t *at = n->full ? n->src + 3 : mpt_items(n->src);
    for (uint32_t i = 0; i < 16; i++)
    {
        out[i] = (mask >> i) & 1 ? n->child[i] : mpt_ref_item(at);
        at += mpt_item_size(at[0]);
    }
}

/* The content of the string item at `item`. */
static inline const uint8_t *mpt_string(const uint8_t *item, uint32_t *len)
{
    uint32_t v = item[0];
    if (v < 0x80)
    {
        *len = 1;
        return item;
    }
    if (v <= 0xb7)
    {
        *len = v - 0x80;
        return item + 1;
    }
    uint32_t ll = v - 0xb7, length = 0;
    for (uint32_t i = 0; i < ll; i++)
        length = length << 8 | item[1 + i];
    *len = length;
    return item + 1 + ll;
}

/* How many of the `count` nibbles from nibble `pos` of path equal those of key from nibble `depth`. */
static uint32_t mpt_common(const uint8_t *path, uint32_t pos, uint32_t count, const uint8_t *key, uint32_t depth)
{
    uint32_t i = 0;
    if (((pos ^ depth) & 1) == 0)
    {
        if ((pos & 1) && count)
        {
            if (mpt_nib(path, pos) != mpt_nib(key, depth))
                return 0;
            i = 1;
        }
        const uint8_t *a = path + ((pos + i) >> 1), *b = key + ((depth + i) >> 1);
#ifdef MPT_UNALIGNED_WORDS
        /* Four bytes at a time, each side read as the two aligned words that hold them. */
        if (i + 8 <= count)
        {
            uint32_t ma = (uintptr_t)a & 3, mb = (uintptr_t)b & 3, la = 8 * ma, lb = 8 * mb;
            const uint32_t *wa = (const uint32_t *)(a - ma), *wb = (const uint32_t *)(b - mb);
            uint32_t ca = wa[0], cb = wb[0];
            do
            {
                uint32_t na = wa[1], nb = wb[1];
                if (((ca >> la) | ((na << 1) << (31 - la))) != ((cb >> lb) | ((nb << 1) << (31 - lb))))
                    break;
                ca = na;
                cb = nb;
                wa++;
                wb++;
                a += 4;
                b += 4;
                i += 8;
            } while (i + 8 <= count);
        }
#endif
        while (i + 2 <= count && *a == *b)
        {
            a++;
            b++;
            i += 2;
        }
    }
    while (i < count && mpt_nib(path, pos + i) == mpt_nib(key, depth + i))
        i++;
    return i;
}

/* ---- nodes --------------------------------------------------------------------------------------- */

static mnode *mpt_new(int kind, int children)
{
    mnode *n = mpt_alloc((uint32_t)(offsetof(mnode, child) + (size_t)children * sizeof(mref)));
    n->state = (uint32_t)kind | 1u << 8 | 1u << 16; /* dirty and changed */
    n->shape = 0;
    n->mask = 0xffff;
    n->dmask = 0xffff;
    n->elen = 0;
    n->slen = 0;
    n->src = 0;
    n->eraw = 0;
    n->hash = 0;
    n->pk = 0;
    n->val = 0;
    n->vlen = 0;
    n->widx = ~0u;
    for (int i = 0; i < children; i++)
    {
        n->child[i].p = 0;
        n->child[i].t = TEmpty;
    }
    return n;
}

static inline __attribute__((always_inline)) int mpt_is_full_branch(const uint8_t *node)
{
    uint32_t bad = (node[0] ^ 0xf9) | (node[1] ^ 0x02) | (node[2] ^ 0x11) | (node[MptFullBranch - 1] ^ 0x80);
#pragma clang loop unroll(full)
    for (int i = 0; i < 16; i++)
        bad |= node[3 + 33 * i] ^ 0xa0;
    return bad == 0;
}

/* Whether node is a branch whose children are hashes and empty slots only, two or more, and that holds no value:
   most of the nodes that are not full branches. Anything else is left to mpt_validate. */
static inline __attribute__((always_inline)) int mpt_is_plain_branch(const uint8_t *node, uint32_t len)
{
    uint32_t header = node[0], items;
    if (header == 0xf8)
    {
        if (len < 58 || node[1] != len - 2)
            return 0;
        items = 2;
    }
    else if (header == 0xf9)
    {
        if (len < 259 || (uint32_t)(node[1] << 8 | node[2]) != len - 3)
            return 0;
        items = 3;
    }
    else
    {
        return 0;
    }
    const uint8_t *at = node + items, *end = node + len;
    uint32_t children = 0;
    /* Written out sixteen times: a loop that returns from its middle is not unrolled. */
#define MPT_PLAIN_SLOT                                                                                               \
    {                                                                                                                \
        if (at >= end)                                                                                               \
            return 0;                                                                                                \
        uint32_t prefix = at[0];                                                                                     \
        if (prefix == 0xa0)                                                                                          \
        {                                                                                                            \
            at += 33;                                                                                                \
            children++;                                                                                              \
        }                                                                                                            \
        else if (prefix == 0x80)                                                                                     \
        {                                                                                                            \
            at++;                                                                                                    \
        }                                                                                                            \
        else                                                                                                         \
        {                                                                                                            \
            return 0;                                                                                                \
        }                                                                                                            \
    }
    MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT
    MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT MPT_PLAIN_SLOT
#undef MPT_PLAIN_SLOT
    if (at >= end)
        return 0;
    return at + 1 == end && at[0] == 0x80 && children >= 2;
}

static mnode *mpt_fail(int e)
{
    MptErr = e;
    return 0;
}

/* The node r, a hash or an embedded encoding, refers to, `depth` nibbles below the root: built from its encoding,
   which is checked the first time; 0 with MptErr set if it cannot be. */
static mnode *mpt_build(mref r, uint32_t depth)
{
    const uint8_t *rlp;
    uint32_t len, widx = ~0u;
    const uint32_t *hash = 0;
    int items;
    if ((r.t & 3) == THash)
    {
        uint32_t h[8];
        mpt_load32(h, (const uint8_t *)r.p);
        int k = mpt_lookup_words(h, 1);
        if (k < 0)
            return mpt_fail(MptMissing);
        widx = (uint32_t)k;
        hash = mpt_w_hash(widx);
        rlp = mpt_w_rlp(widx, &len);
        int info = MptInfo[k];
        if (info == 0)
        {
            if ((len == MptFullBranch && mpt_is_full_branch(rlp)) || mpt_is_plain_branch(rlp, len))
                items = depth >= MptKeyNibbles ? -1 : 0;
            else
                items = mpt_validate(rlp, (int)len, (int)(MptKeyNibbles - depth), depth > 0);
            if (items < 0)
                return mpt_fail(MptInvalid);
            MptInfo[k] = (uint8_t)(items + 1);
        }
        else
        {
            /* Checked already: only what depends on the depth is left. */
            if ((depth > 0 && len < 32) || (info == 1 && depth >= MptKeyNibbles))
                return mpt_fail(MptInvalid);
            items = info - 1;
        }
    }
    else
    {
        rlp = (const uint8_t *)r.p;
        len = r.t >> 8;
        items = mpt_validate(rlp, (int)len, (int)(MptKeyNibbles - depth), 0);
        if (items < 0)
            return mpt_fail(MptInvalid);
    }
    mnode *n;
    if (items == 0)
    {
        n = mpt_alloc(sizeof(mnode));
        n->state = KBranch | (hash ? 1u << 24 : 0);
        n->shape = len == MptFullBranch ? 1u << 16 : 0;
        n->mask = 0;
        n->dmask = 0;
        n->pk = 0;
        n->val = 0;
        n->vlen = 0;
    }
    else
    {
        const uint8_t *at = rlp + items, *pk, *second;
        uint32_t pp = at[0], bytes;
        if (pp < 0x80)
        {
            bytes = 1;
            pk = at;
            second = at + 1;
        }
        else
        {
            bytes = pp - 0x80;
            pk = at + 1;
            second = at + 1 + bytes;
        }
        uint32_t flags = pk[0], odd = (flags >> 4) & 1, nibbles = 2 * bytes - 2 + odd;
        int leaf = (flags & 0x20) != 0;
        if (leaf ? depth + nibbles != MptKeyNibbles : nibbles == 0 || depth + nibbles >= MptKeyNibbles)
            return mpt_fail(MptInvalid);
        n = mpt_alloc((uint32_t)(offsetof(mnode, child) + (leaf ? 0 : sizeof(mref))));
        n->state = (leaf ? KLeaf : KExt) | (hash ? 1u << 24 : 0);
        n->shape = nibbles | (odd ? 1u : 2u) << 8;
        n->mask = 1;
        n->dmask = 1;
        n->pk = pk;
        if (leaf)
        {
            n->val = mpt_string(second, &n->vlen);
        }
        else
        {
            n->val = 0;
            n->vlen = 0;
            n->child[0] = mpt_ref_item(second);
        }
    }
    n->elen = len;
    n->slen = len;
    n->src = rlp;
    n->eraw = rlp;
    n->hash = hash;
    n->widx = widx;
    return n;
}

static inline __attribute__((always_inline)) mnode *mpt_node(mref r, uint32_t depth)
{
    return (r.t & 3) == TNode ? (mnode *)r.p : mpt_build(r, depth);
}

/* The key of the write under way, kept once a new node takes its path from it. */
static const uint8_t *MptKey, *MptKeyKept;

static const uint8_t *mpt_kept_key(void)
{
    if (!MptKeyKept)
        MptKeyKept = mpt_keep(MptKey, 32);
    return MptKeyKept;
}

static mnode *mpt_new_leaf(uint32_t depth, const uint8_t *val, uint32_t vlen)
{
    mnode *n = mpt_new(KLeaf, 0);
    n->pk = mpt_kept_key();
    n->poff = (uint8_t)depth;
    n->plen = (uint8_t)(MptKeyNibbles - depth);
    n->val = val;
    n->vlen = vlen;
    return n;
}

/* count nibbles of the key from depth, then an extension to child. */
static mnode *mpt_new_ext(uint32_t depth, uint32_t count, mref child)
{
    mnode *n = mpt_new(KExt, 1);
    n->pk = mpt_kept_key();
    n->poff = (uint8_t)depth;
    n->plen = (uint8_t)count;
    n->child[0] = child;
    return n;
}

/* The path of n behind the `count` nibbles from nibble `pos` of `front`. */
static void mpt_prepend(mnode *n, const uint8_t *front, uint32_t pos, uint32_t count)
{
    uint32_t total = count + n->plen, off = total & 1, size = (total + off + 1) / 2 + 1;
    uint8_t *b = mpt_alloc(size);
    for (uint32_t i = 0; i < size; i++)
        b[i] = 0;
    for (uint32_t i = 0; i < total; i++)
    {
        uint32_t v = i < count ? mpt_nib(front, pos + i) : mpt_nib(n->pk, n->poff + i - count);
        uint32_t at = off + i;
        b[at >> 1] |= (uint8_t)((at & 1) ? v : v << 4);
    }
    n->pk = b;
    n->poff = (uint8_t)off;
    n->plen = (uint8_t)total;
    n->dirty = 1;
    n->changed = 1;
}

/* ---- reads --------------------------------------------------------------------------------------- */

static int mpt_get(mtrie *t, const uint8_t *key, const uint8_t **val, uint32_t *vlen)
{
    mref r = t->root;
    if ((r.t & 3) == TEmpty)
        return MptAbsent;
    mnode *n = (mnode *)r.p;
    if ((r.t & 3) != TNode)
    {
        n = mpt_build(r, 0);
        if (!n)
            return MptErr;
        t->root = mpt_ref_node(n);
    }
    uint32_t depth = 0;
    for (;;)
    {
        if (n->kind == KBranch)
        {
            if (depth >= MptKeyNibbles)
                return MptInvalid;
            uint32_t i = mpt_nib(key, depth);
            depth++;
            mref c = mpt_child(n, i);
            uint32_t tag = c.t & 3;
            if (tag == TNode)
            {
                n = (mnode *)c.p;
                continue;
            }
            if (tag == TEmpty)
                return MptAbsent;
            mnode *b = mpt_build(c, depth);
            if (!b)
                return MptErr;
            n->child[i] = mpt_ref_node(b);
            n->mask |= 1u << i;
            n = b;
            continue;
        }
        if (mpt_common(n->pk, n->poff, n->plen, key, depth) != n->plen)
            return MptAbsent;
        if (n->kind == KLeaf)
        {
            if (depth + n->plen != MptKeyNibbles)
                return MptInvalid;
            *val = n->val;
            *vlen = n->vlen;
            return MptOk;
        }
        depth += n->plen;
        if (depth >= MptKeyNibbles)
            return MptInvalid;
        mref c = n->child[0];
        mnode *b = (mnode *)c.p;
        if ((c.t & 3) != TNode)
        {
            b = mpt_build(c, depth);
            if (!b)
                return MptErr;
            n->child[0] = mpt_ref_node(b);
        }
        /* The child of an extension is a branch: a leaf or another extension there would be one path. */
        if (b->kind != KBranch)
            return MptInvalid;
        n = b;
    }
}

/* ---- writes -------------------------------------------------------------------------------------- */

/* A step of a write's way down: a branch and the child taken, or an extension, where child is 16. */
typedef struct
{
    mnode *n;
    uint8_t child, depth;
} mstep;

/* The trie of root with the key set to the vlen bytes at val, or without the key when vlen is 0: its new root.
   MptChanged says whether anything changed; MptErr that a node it needed is missing or malformed. */
static mref mpt_set(mref root, const uint8_t *val, uint32_t vlen)
{
    const uint8_t *key = MptKey;
    mstep path[MptKeyNibbles + 2];
    int top = 0;
    uint32_t depth = 0;
    mref r = root, c;
    /* Down to where the key's path ends or leaves the trie: c is what takes that place. */
    for (;;)
    {
        if ((r.t & 3) == TEmpty)
        {
            if (!vlen)
                goto unchanged;
            c = mpt_ref_node(mpt_new_leaf(depth, val, vlen));
            break;
        }
        mnode *n = (mnode *)r.p;
        if ((r.t & 3) != TNode)
        {
            n = mpt_build(r, depth);
            if (!n)
                return root;
            if (top)
            {
                /* The branch above keeps what was built. An extension's child is built before the step down. */
                mstep *up = &path[top - 1];
                up->n->child[up->child] = mpt_ref_node(n);
                up->n->mask |= 1u << up->child;
            }
            else
            {
                root = mpt_ref_node(n);
            }
        }
        if (n->kind == KBranch)
        {
            if (depth >= MptKeyNibbles)
            {
                MptErr = MptInvalid;
                return root;
            }
            uint32_t i = mpt_nib(key, depth);
            path[top].n = n;
            path[top].child = (uint8_t)i;
            path[top].depth = (uint8_t)depth;
            top++;
            r = mpt_child(n, i);
            depth++;
            continue;
        }
        uint32_t cp = mpt_common(n->pk, n->poff, n->plen, key, depth);
        if (n->kind == KLeaf)
        {
            if (depth + n->plen != MptKeyNibbles)
            {
                MptErr = MptInvalid;
                return root;
            }
            if (cp == n->plen)
            {
                if (!vlen)
                {
                    c.p = 0;
                    c.t = TEmpty;
                    break;
                }
                n->val = val;
                n->vlen = vlen;
                n->dirty = 1;
                n->changed = 1;
                c = mpt_ref_node(n);
                break;
            }
            if (!vlen)
                goto unchanged;
            mnode *b = mpt_new(KBranch, 16);
            uint32_t at = mpt_nib(n->pk, n->poff + cp);
            n->poff = (uint8_t)(n->poff + cp + 1);
            n->plen = (uint8_t)(n->plen - cp - 1);
            n->dirty = 1;
            n->changed = 1;
            b->child[at] = mpt_ref_node(n);
            b->child[mpt_nib(key, depth + cp)] = mpt_ref_node(mpt_new_leaf(depth + cp + 1, val, vlen));
            c = cp ? mpt_ref_node(mpt_new_ext(depth, cp, mpt_ref_node(b))) : mpt_ref_node(b);
            break;
        }
        /* An extension. */
        if (cp == n->plen)
        {
            if (depth + cp >= MptKeyNibbles)
            {
                MptErr = MptInvalid;
                return root;
            }
            mnode *below = mpt_node(n->child[0], depth + cp);
            if (!below)
                return root;
            /* The child of an extension is a branch: a leaf or another extension there would be one path. */
            if (below->kind != KBranch)
            {
                MptErr = MptInvalid;
                return root;
            }
            n->child[0] = mpt_ref_node(below);
            path[top].n = n;
            path[top].child = 16;
            path[top].depth = (uint8_t)depth;
            top++;
            depth += cp;
            r = mpt_ref_node(below);
            continue;
        }
        if (!vlen)
            goto unchanged;
        mnode *b = mpt_new(KBranch, 16);
        uint32_t at = mpt_nib(n->pk, n->poff + cp), rest = n->plen - cp - 1;
        if (rest)
        {
            n->poff = (uint8_t)(n->poff + cp + 1);
            n->plen = (uint8_t)rest;
            n->dirty = 1;
            n->changed = 1;
            b->child[at] = mpt_ref_node(n);
        }
        else
        {
            b->child[at] = n->child[0];
        }
        b->child[mpt_nib(key, depth + cp)] = mpt_ref_node(mpt_new_leaf(depth + cp + 1, val, vlen));
        c = cp ? mpt_ref_node(mpt_new_ext(depth, cp, mpt_ref_node(b))) : mpt_ref_node(b);
        break;
    }
    MptChanged = 1;
    /* Back up: each node on the way takes its changed child, and folds where a child went. */
    while (top > 0)
    {
        mstep *s = &path[--top];
        mnode *n = s->n;
        if (s->child == 16)
        {
            if ((c.t & 3) == TEmpty)
                continue;
            mnode *cn = (mnode *)c.p;
            if (cn->kind != KBranch)
            {
                /* The branch folded into a leaf or an extension, which takes this path in front of its own. */
                mpt_prepend(cn, n->pk, n->poff, n->plen);
                continue;
            }
            n->child[0] = c;
            n->dirty = 1;
            n->changed = 1;
            c = mpt_ref_node(n);
            continue;
        }
        uint32_t i = s->child;
        n->child[i] = c;
        n->mask |= 1u << i;
        n->dmask |= 1u << i;
        n->dirty = 1;
        n->changed = 1;
        if ((c.t & 3) != TEmpty)
        {
            c = mpt_ref_node(n);
            continue;
        }
        /* A child went: with one left, the branch folds into it. */
        mref all[16];
        mpt_children(n, all);
        int only = -1, count = 0;
        for (int j = 0; j < 16; j++)
        {
            if ((all[j].t & 3) != TEmpty)
            {
                only = j;
                count++;
            }
        }
        if (count != 1)
        {
            c = mpt_ref_node(n);
            continue;
        }
        mnode *o = mpt_node(all[only], (uint32_t)s->depth + 1);
        if (!o)
            return root;
        uint8_t front = (uint8_t)(only << 4);
        if (o->kind == KBranch)
        {
            mnode *x = mpt_new(KExt, 1);
            uint8_t *b = mpt_alloc(1);
            b[0] = front;
            x->pk = b;
            x->poff = 0;
            x->plen = 1;
            x->child[0] = mpt_ref_node(o);
            c = mpt_ref_node(x);
            continue;
        }
        mpt_prepend(o, &front, 0, 1);
        c = mpt_ref_node(o);
    }
    return c;
unchanged:
    return root;
}

/* ---- the root ------------------------------------------------------------------------------------ */

/* Room for an encoding at any alignment modulo 32, with 64 bytes to spare on either side. */
static uint32_t MptRoom[(64 + 32 + MptBufferBytes + 64) / 4] __attribute__((aligned(32)));

static inline uint32_t mpt_ref_size(mref r)
{
    switch (r.t & 3)
    {
    case TEmpty:
        return 1;
    case THash:
        return 33;
    case TEmbed:
        return r.t >> 8;
    default:
    {
        const mnode *n = (const mnode *)r.p;
        return n->hashed ? 33 : n->elen;
    }
    }
}

static inline uint8_t *mpt_put_ref(uint8_t *at, mref r)
{
    switch (r.t & 3)
    {
    case TEmpty:
        *at = 0x80;
        return at + 1;
    case THash:
    {
        uint32_t h[8];
        mpt_load32(h, (const uint8_t *)r.p);
        *at = 0xa0;
        mpt_store32(at + 1, h);
        return at + 33;
    }
    case TEmbed:
    {
        const uint8_t *from = (const uint8_t *)r.p;
        uint32_t len = r.t >> 8;
        for (uint32_t i = 0; i < len; i++)
            at[i] = from[i];
        return at + len;
    }
    default:
    {
        const mnode *n = (const mnode *)r.p;
        if (n->hashed)
        {
            *at = 0xa0;
            mpt_store32(at + 1, n->hash);
            return at + 33;
        }
        const uint8_t *from = n->eraw;
        uint32_t len = n->elen;
        for (uint32_t i = 0; i < len; i++)
            at[i] = from[i];
        return at + len;
    }
    }
}

static inline uint8_t *mpt_put_list(uint8_t *at, uint32_t payload)
{
    if (payload < 56)
    {
        *at = (uint8_t)(0xc0 + payload);
        return at + 1;
    }
    if (payload < 256)
    {
        at[0] = 0xf8;
        at[1] = (uint8_t)payload;
        return at + 2;
    }
    at[0] = 0xf9;
    at[1] = (uint8_t)(payload >> 8);
    at[2] = (uint8_t)payload;
    return at + 3;
}

/* The path of n as a hex-prefix string item. */
static uint8_t *mpt_put_path(uint8_t *at, const mnode *n, int leaf)
{
    uint32_t odd = n->plen & 1, bytes = n->plen / 2 + 1;
    uint32_t first = (leaf ? 0x20 : 0) | (odd ? 0x10 | mpt_nib(n->pk, n->poff) : 0);
    if (bytes == 1 && first < 0x80)
    {
        *at = (uint8_t)first;
        return at + 1;
    }
    *at++ = (uint8_t)(0x80 + bytes);
    *at++ = (uint8_t)first;
    uint32_t pos = n->poff + odd;
    if ((pos & 1) == 0)
    {
        memcpy(at, n->pk + (pos >> 1), bytes - 1);
    }
    else
    {
        for (uint32_t i = 0; i + 1 < bytes; i++)
            at[i] = (uint8_t)(mpt_nib(n->pk, pos + 2 * i) << 4 | mpt_nib(n->pk, pos + 2 * i + 1));
    }
    return at + bytes - 1;
}

static inline uint32_t mpt_path_size(const mnode *n, int leaf)
{
    uint32_t odd = n->plen & 1, bytes = n->plen / 2 + 1;
    uint32_t first = (leaf ? 0x20 : 0) | (odd ? 0x10 | mpt_nib(n->pk, n->poff) : 0);
    return bytes == 1 && first < 0x80 ? 1 : 1 + bytes;
}

/* The encoding of n, whose children are all computed, at `start`; its length, or 0 with MptErr set. */
static uint32_t mpt_encode(const mnode *n, uint8_t *start)
{
    uint8_t *at;
    if (n->kind == KBranch)
    {
        mref all[16];
        mpt_children(n, all);
        uint32_t payload = 1;
        for (int i = 0; i < 16; i++)
            payload += mpt_ref_size(all[i]);
        at = mpt_put_list(start, payload);
        for (int i = 0; i < 16; i++)
            at = mpt_put_ref(at, all[i]);
        *at++ = 0x80;
        return (uint32_t)(at - start);
    }
    if (n->kind == KExt)
    {
        at = mpt_put_list(start, mpt_path_size(n, 0) + mpt_ref_size(n->child[0]));
        at = mpt_put_path(at, n, 0);
        at = mpt_put_ref(at, n->child[0]);
        return (uint32_t)(at - start);
    }
    const uint8_t *val = n->val;
    uint32_t vlen = n->vlen, single = vlen == 1 && val[0] < 0x80;
    uint32_t vsize = single ? 1 : vlen < 56 ? 1 + vlen : vlen < 256 ? 2 + vlen : 3 + vlen;
    uint32_t payload = mpt_path_size(n, 1) + vsize;
    if (payload + 3 > MptBufferBytes)
    {
        MptErr = MptTooLong;
        return 0;
    }
    at = mpt_put_list(start, payload);
    at = mpt_put_path(at, n, 1);
    if (!single)
    {
        if (vlen < 56)
        {
            *at++ = (uint8_t)(0x80 + vlen);
        }
        else if (vlen < 256)
        {
            *at++ = 0xb8;
            *at++ = (uint8_t)vlen;
        }
        else
        {
            *at++ = 0xb9;
            *at++ = (uint8_t)(vlen >> 8);
            *at++ = (uint8_t)vlen;
        }
    }
    memcpy(at, val, vlen);
    return (uint32_t)(at + vlen - start);
}

/* The digest of a branch built from an encoding, whose `changed` children were hashes there and are hashes still
   and whose `emptied` ones were empty there: the encoding with the new hashes over the old. 0 if n is not so. */
static int mpt_finish_spliced(mnode *n, uint32_t changed, uint32_t emptied, uint32_t *out)
{
    const uint8_t *src = n->src;
    uint8_t *to = (uint8_t *)MptRoom + 64 + ((uintptr_t)src & 31);
    if (n->full)
    {
        if (emptied)
            return 0;
        uint32_t first = (uint32_t)__builtin_ctz(changed);
        int block = (int)((4 + 33 * first) / 136);
        mpt_copy(to + 136 * block, src + 136 * block, MptFullBranch - 136 * (uint32_t)block);
        for (uint32_t m = changed; m; m &= m - 1)
        {
            uint32_t i = (uint32_t)__builtin_ctz(m);
            mpt_store32(to + 4 + 33 * i, ((const mnode *)n->child[i].p)->hash);
        }
        mpt_hash_edited(n->widx, block, src, to, out);
        n->elen = MptFullBranch;
        return 1;
    }
    /* Where the changed children lie, and that an emptied one was empty already. */
    const uint8_t *at = mpt_items(src);
    uint32_t offset[16];
    for (uint32_t i = 0; (changed | emptied) >> i; i++)
    {
        uint32_t prefix = at[0];
        if ((changed >> i) & 1)
        {
            if (prefix != 0xa0)
                return 0;
            offset[i] = (uint32_t)(at + 1 - src);
        }
        else if (((emptied >> i) & 1) && prefix != 0x80)
        {
            return 0;
        }
        at += mpt_item_size(prefix);
    }
    uint32_t len = n->slen;
    mpt_copy(to, src, len);
    for (uint32_t i = 0, m = changed; m; i++, m >>= 1)
    {
        if (m & 1)
            mpt_store32(to + offset[i], ((const mnode *)n->child[i].p)->hash);
    }
    mpt_keccak(to, len, out);
    n->elen = len;
    return 1;
}

/* Computes the encoding of n and of every changed node below it. */
static void mpt_finish(mnode *n, int is_root)
{
    uint32_t changed = 0, emptied = 0, unhashed = 0;
    if (n->kind == KBranch)
    {
        for (uint32_t m = n->dmask; m; m &= m - 1)
        {
            uint32_t i = (uint32_t)__builtin_ctz(m);
            mref c = n->child[i];
            uint32_t tag = c.t & 3;
            if (tag == TNode)
            {
                mnode *cn = (mnode *)c.p;
                if (cn->dirty)
                {
                    mpt_finish(cn, 0);
                    if (MptErr)
                        return;
                }
                if (cn->changed)
                {
                    if (cn->hashed)
                        changed |= 1u << i;
                    else
                        unhashed = 1;
                }
            }
            else if (tag == TEmpty)
            {
                emptied |= 1u << i;
            }
        }
    }
    else if (n->kind == KExt)
    {
        mref c = n->child[0];
        if ((c.t & 3) == TNode && ((mnode *)c.p)->dirty)
        {
            mpt_finish((mnode *)c.p, 0);
            if (MptErr)
                return;
        }
    }
    uint32_t *h = mpt_alloc(32);
    if (changed && !unhashed && n->src && n->hashed && mpt_finish_spliced(n, changed, emptied, h))
    {
        n->hash = h;
    }
    else
    {
        uint8_t *start = (uint8_t *)MptRoom + 64;
        uint32_t len = mpt_encode(n, start);
        if (MptErr)
            return;
        if (len >= 32 || is_root)
        {
            mpt_keccak(start, len, h);
            n->hash = h;
            n->hashed = 1;
        }
        else
        {
            n->eraw = mpt_keep(start, len);
            n->hashed = 0;
        }
        n->elen = len;
    }
    n->dirty = 0;
}

static int mpt_root(mtrie *t, uint32_t out[8])
{
    const uint32_t *h;
    uint32_t tag = t->root.t & 3;
    if (tag == TEmpty)
    {
        h = MptEmptyRoot;
    }
    else if (tag == TNode)
    {
        mnode *n = (mnode *)t->root.p;
        if (n->dirty || !n->hashed)
        {
            MptErr = 0;
            mpt_finish(n, 1);
            if (MptErr)
                return MptErr;
        }
        h = n->hash;
    }
    else
    {
        h = t->rooth;
    }
    for (int i = 0; i < 8; i++)
        out[i] = h[i];
    return MptOk;
}

static mtrie *mpt_open(const uint32_t root[8])
{
    mtrie *t = mpt_alloc(sizeof(mtrie));
    uint32_t differs = 0;
    for (int i = 0; i < 8; i++)
    {
        t->rooth[i] = root[i];
        differs |= root[i] ^ MptEmptyRoot[i];
    }
    if (differs)
    {
        t->root.p = (uintptr_t)t->rooth;
        t->root.t = THash;
    }
    else
    {
        t->root.p = 0;
        t->root.t = TEmpty;
    }
    return t;
}

static int mpt_write(mtrie *t, const uint8_t *key, const uint8_t *val, uint32_t vlen)
{
    MptErr = 0;
    MptChanged = 0;
    MptKey = key;
    MptKeyKept = 0;
    mref r = mpt_set(t->root, vlen ? mpt_keep(val, vlen) : 0, vlen);
    t->root = r;
    return MptErr;
}

/* ---- keys ---------------------------------------------------------------------------------------- */

enum
{
    MptMemoSize = 8192
};

static struct
{
    uint32_t pre[8];
    uint32_t len;
    uint32_t h[8];
} MptMemo[MptMemoSize];

/* The keccak of the len bytes, 20 or 32, in the words at pre, which are zero past them. */
static const uint8_t *mpt_key(const uint32_t pre[8], uint32_t len)
{
    uint32_t x = pre[0] ^ pre[2] ^ pre[4] ^ pre[7];
    x ^= x >> 15;
    x ^= x >> 7;
    __typeof__(MptMemo[0]) *e = &MptMemo[x & (MptMemoSize - 1)];
    if (e->len == len
        && ((e->pre[0] ^ pre[0]) | (e->pre[1] ^ pre[1]) | (e->pre[2] ^ pre[2]) | (e->pre[3] ^ pre[3])
            | (e->pre[4] ^ pre[4]) | (e->pre[5] ^ pre[5]) | (e->pre[6] ^ pre[6]) | (e->pre[7] ^ pre[7])) == 0)
        return (const uint8_t *)e->h;
    for (int i = 0; i < 8; i++)
        e->pre[i] = pre[i];
    e->len = len;
    mpt_keccak_key(e->pre, len, e->h);
    return (const uint8_t *)e->h;
}

static inline uint32_t mpt_bswap(uint32_t x)
{
#ifdef MPT_REV8
    return MPT_REV8(x);
#else
    return __builtin_bswap32(x);
#endif
}

static const uint8_t *mpt_account_key(const uint8_t *address)
{
    uint32_t pre[8];
    if (((uintptr_t)address & 3) == 0)
    {
        const uint32_t *w = (const uint32_t *)address;
#pragma clang loop unroll(full)
        for (int i = 0; i < 5; i++)
            pre[i] = w[i];
    }
    else
    {
        memcpy(pre, address, 20);
    }
    pre[5] = pre[6] = pre[7] = 0;
    return mpt_key(pre, 20);
}

/* The key of a storage slot: the keccak of its index, big-endian. index holds eight little-endian words. */
static const uint8_t *mpt_slot_key(const uint32_t *index)
{
    uint32_t pre[8];
#pragma clang loop unroll(full)
    for (int i = 0; i < 8; i++)
        pre[i] = mpt_bswap(index[7 - i]);
    return mpt_key(pre, 32);
}

/* ---- the calls ----------------------------------------------------------------------------------- */

void *airbender_mpt_open(const uint8_t *root)
{
    uint32_t h[8];
    memcpy(h, root, 32);
    return mpt_open(h);
}

void airbender_mpt_clear(void *trie)
{
    mtrie *t = trie;
    t->root.p = 0;
    t->root.t = TEmpty;
}

/* Whether the trie's root is the empty tree or a witness node. */
int airbender_mpt_has_root(const uint8_t *root)
{
    uint32_t h[8], differs = 0;
    memcpy(h, root, 32);
    for (int i = 0; i < 8; i++)
        differs |= h[i] ^ MptEmptyRoot[i];
    return !differs || mpt_lookup_words(h, 0) >= 0;
}

int airbender_mpt_root(void *trie, uint8_t *out)
{
    uint32_t h[8];
    int e = mpt_root(trie, h);
    if (e)
        return e;
    memcpy(out, h, 32);
    return MptOk;
}

/* The account's encoding where it lies: valid until the trie is written. */
int airbender_mpt_get_account(void *trie, const uint8_t *address, const uint8_t **value, uint32_t *len)
{
    return mpt_get(trie, mpt_account_key(address), value, len);
}

/* Sets the account to the len bytes at value, or removes it when len is 0. */
int airbender_mpt_set_account(void *trie, const uint8_t *address, const uint8_t *value, uint32_t len)
{
    return mpt_write(trie, mpt_account_key(address), value, len);
}

/* The big-endian bytes of a string item of at most `max` bytes at *at, as little-endian words in out, which
   holds max / 4 of them; *at moves past the item. 0 if it is not such an item. */
static inline int mpt_take_number(const uint8_t **at, const uint8_t *end, uint32_t *out, uint32_t max)
{
    const uint8_t *p = *at;
    for (uint32_t i = 0; i < max / 4; i++)
        out[i] = 0;
    if (p >= end)
        return 0;
    uint32_t v = p[0], count;
    if (v < 0x80)
    {
        out[0] = v;
        *at = p + 1;
        return 1;
    }
    count = v - 0x80;
    if (count > max || p + 1 + count > end)
        return 0;
    p++;
    for (uint32_t i = 0; i < count; i++)
    {
        uint32_t to = count - 1 - i;
        out[to >> 2] |= (uint32_t)p[i] << (8 * (to & 3));
    }
    *at = p + count;
    return 1;
}

static const uint32_t MptEmptyCode[8] = {0x0146d2c5, 0x3c23f786, 0xb27d7e92, 0xc003c7dc,
                                         0x53b600e5, 0x3b2782ca, 0x04d8fa7b, 0x70a4855d};

/* The account as fields: two words of nonce, eight of balance, the storage root, the code hash, then a word whose
   bit 0 says the storage is empty and bit 1 that the code is. 2, with its encoding at *value, when the encoding is
   not the four items an account has in the state trie. */
int airbender_mpt_get_account_fields(void *trie, const uint8_t *address, uint32_t *out, const uint8_t **value,
                                     uint32_t *len)
{
    int e = mpt_get(trie, mpt_account_key(address), value, len);
    if (e)
        return e;
    const uint8_t *at = *value, *end = at + *len;
    if (*len < 70 || at[0] != 0xf8 || (uint32_t)at[1] + 2 != *len)
        return 2;
    at += 2;
    if (!mpt_take_number(&at, end, out, 8) || !mpt_take_number(&at, end, out + 2, 32))
        return 2;
    if (end - at != 66 || at[0] != 0xa0 || at[33] != 0xa0)
        return 2;
    mpt_load32(out + 10, at + 1);
    mpt_load32(out + 18, at + 34);
    uint32_t storage = 0, code = 0;
#pragma clang loop unroll(full)
    for (int i = 0; i < 8; i++)
    {
        storage |= out[10 + i] ^ MptEmptyRoot[i];
        code |= out[18 + i] ^ MptEmptyCode[i];
    }
    out[26] = (storage == 0) | (code == 0) << 1;
    return MptOk;
}

static inline uint8_t *mpt_put_number(uint8_t *at, const uint32_t *words, int count)
{
    int top = count - 1;
    while (top >= 0 && words[top] == 0)
        top--;
    if (top < 0)
    {
        *at = 0x80;
        return at + 1;
    }
    uint32_t w = words[top];
    uint32_t bytes = (uint32_t)top * 4 + (w >> 24 ? 4 : w >> 16 ? 3 : w >> 8 ? 2 : 1);
    if (!(bytes == 1 && w < 0x80))
        *at++ = (uint8_t)(0x80 + bytes);
    for (uint32_t i = 0; i < bytes; i++)
    {
        uint32_t from = bytes - 1 - i;
        at[i] = (uint8_t)(words[from >> 2] >> (8 * (from & 3)));
    }
    return at + bytes;
}

/* Sets the account from its fields: two words of nonce, eight of balance, and the 32 bytes of the storage root
   and of the code hash. */
int airbender_mpt_set_account_fields(void *trie, const uint8_t *address, const uint32_t *nonce,
                                     const uint32_t *balance, const uint8_t *storage_root, const uint8_t *code_hash)
{
    uint32_t words[32], h[8];
    uint8_t *enc = (uint8_t *)words, *at = enc + 2;
    at = mpt_put_number(at, nonce, 2);
    at = mpt_put_number(at, balance, 8);
    *at++ = 0xa0;
    mpt_load32(h, storage_root);
    mpt_store32(at, h);
    at += 32;
    *at++ = 0xa0;
    mpt_load32(h, code_hash);
    mpt_store32(at, h);
    at += 32;
    enc[0] = 0xf8;
    enc[1] = (uint8_t)(at - enc - 2);
    return mpt_write(trie, mpt_account_key(address), enc, (uint32_t)(at - enc));
}

/* The trie of the account's storage: of its storage root, or empty when there is no such account. 0 with *status
   negative when a witness node is missing or malformed, or 2 when the account's encoding is not the usual one. */
void *airbender_mpt_open_storage(void *trie, const uint8_t *address, int *status)
{
    const uint8_t *value;
    uint32_t len, h[8];
    int e = mpt_get(trie, mpt_account_key(address), &value, &len);
    *status = e < 0 ? e : MptOk;
    if (e < 0)
        return 0;
    if (e == MptAbsent)
        return mpt_open(MptEmptyRoot);
    if (len < 70 || value[len - 66] != 0xa0 || value[len - 33] != 0xa0)
    {
        /* Not the encoding of an account in the state trie: for the caller to read. */
        *status = 2;
        return 0;
    }
    mpt_load32(h, value + len - 65);
    return mpt_open(h);
}

/* The slot's value as eight little-endian words, zero when absent. */
int airbender_mpt_get_slot(void *trie, const uint32_t *index, uint32_t *value)
{
    const uint8_t *v;
    uint32_t len;
#pragma clang loop unroll(full)
    for (int i = 0; i < 8; i++)
        value[i] = 0;
    int e = mpt_get(trie, mpt_slot_key(index), &v, &len);
    if (e)
        return e == MptAbsent ? MptOk : e;
    /* The leaf holds the RLP of the value's big-endian bytes. */
    if (len == 0)
        return MptInvalid;
    uint32_t count;
    if (v[0] < 0x80)
    {
        if (len != 1)
            return MptInvalid;
        count = 1;
    }
    else
    {
        count = v[0] - 0x80u;
        if (count > 32 || len != 1 + count)
            return MptInvalid;
        v++;
    }
    for (uint32_t i = 0; i < count; i++)
    {
        uint32_t at = count - 1 - i;
        value[at >> 2] |= (uint32_t)v[i] << (8 * (at & 3));
    }
    return MptOk;
}

/* Sets the slot to the eight little-endian words at value, or removes it when they are zero. */
int airbender_mpt_set_slot(void *trie, const uint32_t *index, const uint32_t *value)
{
    uint32_t words[12];
    uint8_t *enc = (uint8_t *)words;
    int top = 7;
    while (top >= 0 && value[top] == 0)
        top--;
    uint32_t len = 0;
    if (top >= 0)
    {
        uint32_t w = value[top];
        uint32_t count = (uint32_t)top * 4 + (w >> 24 ? 4 : w >> 16 ? 3 : w >> 8 ? 2 : 1);
        uint8_t *at = enc;
        if (!(count == 1 && w < 0x80))
            *at++ = (uint8_t)(0x80 + count);
        for (uint32_t i = 0; i < count; i++)
        {
            uint32_t from = count - 1 - i;
            at[i] = (uint8_t)(value[from >> 2] >> (8 * (from & 3)));
        }
        len = (uint32_t)(at + count - enc);
    }
    return mpt_write(trie, mpt_slot_key(index), enc, len);
}

/* For tests: the calls on a key given as its 32 bytes. */
int airbender_mpt_get_raw(void *trie, const uint8_t *key, const uint8_t **value, uint32_t *len)
{
    return mpt_get(trie, key, value, len);
}

int airbender_mpt_set_raw(void *trie, const uint8_t *key, const uint8_t *value, uint32_t len)
{
    return mpt_write(trie, key, value, len);
}
