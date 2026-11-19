/*
 * The inversion hint of src/modexp.c on the host: modexp.c and the u256 delegation's C model are compiled in, and this
 * file owns the field hint tape and records what modexp gives out. Every case compares with the plain exponentiation
 * (no tape), so a hint changes how the result is found, never the result.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../src/modexp.c"
#include "../../src/u256_soft.c"

field_hints FieldHints;
static u256d Emitted[8];
static int EmittedCount;

void field_hint_emit(const u256d *v)
{
    if (EmittedCount < 8)
        Emitted[EmittedCount] = *v;
    EmittedCount++;
}

static int Failures;

static void check(int ok, const char *name)
{
    if (!ok)
    {
        printf("FAIL %s\n", name);
        Failures++;
    }
}

/* 32 big-endian bytes of the eight little-endian limbs at v. */
static void to_bytes(uint8_t out[32], const uint32_t v[8])
{
    for (int i = 0; i < 32; i++)
        out[i] = (uint8_t)(v[(31 - i) / 4] >> (8 * ((31 - i) % 4)));
}

static void from_bytes(uint32_t v[8], const uint8_t in[32])
{
    memset(v, 0, 32);
    for (int i = 0; i < 32; i++)
        v[(31 - i) / 4] |= (uint32_t)in[i] << (8 * ((31 - i) % 4));
}

static void add_small(uint32_t r[8], const uint32_t a[8], int64_t delta)
{
    int64_t carry = delta;
    for (int i = 0; i < 8; i++)
    {
        int64_t t = (int64_t)a[i] + carry;
        r[i] = (uint32_t)t;
        carry = t >> 32;
    }
}

/* The plain result, with no tape. */
static void plain(uint8_t out[32], const uint8_t *base, size_t base_len, const uint8_t *exp, size_t exp_len, const uint8_t m[32])
{
    FieldHints.active = 0;
    zkvm_modexp(base, base_len, exp, exp_len, m, 32, out);
}

/* Runs with a tape of count values from tape; returns how many it took. */
static long hinted(uint8_t out[32], const uint8_t *base, size_t base_len, const uint8_t *exp, size_t exp_len, const uint8_t m[32],
                   const u256d *tape, int count)
{
    static u256d area[4] __attribute__((aligned(32)));
    memcpy(area, tape, (size_t)count * sizeof(u256d));
    FieldHints.active = 1;
    FieldHints.next = area;
    FieldHints.end = area + count;
    EmittedCount = 0;
    zkvm_modexp(base, base_len, exp, exp_len, m, 32, out);
    return FieldHints.next - area;
}

static void case_prime(const char *name, const uint32_t prime[8])
{
    uint8_t m[32], e[33], x[32], expected[32], out[32];
    uint32_t pm2[8], value[8], inverse[8];
    to_bytes(m, prime);
    add_small(pm2, prime, -2);
    e[0] = 0;
    to_bytes(e + 1, pm2);
    /* x = (m - 1) / 3 + 12345, below m */
    for (int i = 0; i < 8; i++)
        value[i] = prime[i] / 3 + (i == 0 ? 12345u : 0);
    to_bytes(x, value);
    plain(expected, x, 32, e + 1, 32, m);
    from_bytes(inverse, expected);

    char label[128];
    u256d tape[1];
    memcpy(tape[0].w, inverse, 32);
    snprintf(label, sizeof label, "%s: right hint is taken", name);
    check(hinted(out, x, 32, e + 1, 32, m, tape, 1) == 1 && !memcmp(out, expected, 32) && EmittedCount == 0, label);

    snprintf(label, sizeof label, "%s: exponent with a leading zero byte", name);
    check(hinted(out, x, 32, e, 33, m, tape, 1) == 1 && !memcmp(out, expected, 32) && EmittedCount == 0, label);

    add_small(tape[0].w, inverse, 1);
    snprintf(label, sizeof label, "%s: wrong hint is computed and given out", name);
    check(hinted(out, x, 32, e + 1, 32, m, tape, 1) == 1 && !memcmp(out, expected, 32) && EmittedCount == 1
              && !memcmp(Emitted[0].w, inverse, 32), label);

    memcpy(tape[0].w, prime, 32);
    snprintf(label, sizeof label, "%s: hint not below m is computed", name);
    check(hinted(out, x, 32, e + 1, 32, m, tape, 1) == 1 && !memcmp(out, expected, 32) && EmittedCount == 1, label);

    snprintf(label, sizeof label, "%s: empty tape is computed and given out", name);
    check(hinted(out, x, 32, e + 1, 32, m, tape, 0) == 0 && !memcmp(out, expected, 32) && EmittedCount == 1
              && !memcmp(Emitted[0].w, inverse, 32), label);

    uint8_t zero[1] = {0};
    snprintf(label, sizeof label, "%s: zero base gives zero and takes nothing", name);
    uint8_t zeros[32] = {0};
    check(hinted(out, zero, 1, e + 1, 32, m, tape, 1) == 0 && !memcmp(out, zeros, 32) && EmittedCount == 0, label);

    /* x + m: reduced before the check */
    uint32_t over[8];
    uint64_t carry = 0;
    for (int i = 0; i < 8; i++)
    {
        uint64_t t = (uint64_t)value[i] + prime[i] + carry;
        over[i] = (uint32_t)t;
        carry = t >> 32;
    }
    uint8_t xo[33];
    xo[0] = (uint8_t)carry;
    to_bytes(xo + 1, over);
    memcpy(tape[0].w, inverse, 32);
    snprintf(label, sizeof label, "%s: base above m", name);
    check(hinted(out, xo, 33, e + 1, 32, m, tape, 1) == 1 && !memcmp(out, expected, 32), label);

    uint32_t pm3[8];
    add_small(pm3, prime, -3);
    uint8_t e3[32], other[32];
    to_bytes(e3, pm3);
    plain(other, x, 32, e3, 32, m);
    snprintf(label, sizeof label, "%s: another exponent takes no hint", name);
    check(hinted(out, x, 32, e3, 32, m, tape, 1) == 0 && !memcmp(out, other, 32) && EmittedCount == 0, label);
}

int main(void)
{
    const char *names[] = {"bn254 r", "bn254 p", "secp256k1 p", "secp256k1 n", "p256 p", "p256 n", "bls12-381 r"};
    for (size_t p = 0; p < sizeof(InversePrimes) / sizeof(InversePrimes[0]); p++)
        case_prime(names[p], InversePrimes[p]);

    /* An odd modulus that is not one of the primes: no hint, even with exponent m - 2. */
    uint32_t composite[8], cm2[8], value[8] = {7};
    add_small(composite, InversePrimes[0], 2);
    add_small(cm2, composite, -2);
    uint8_t m[32], e[32], x[32], expected[32], out[32];
    to_bytes(m, composite);
    to_bytes(e, cm2);
    to_bytes(x, value);
    plain(expected, x, 32, e, 32, m);
    u256d tape[1] = {{{1}}};
    check(hinted(out, x, 32, e, 32, m, tape, 1) == 0 && !memcmp(out, expected, 32) && EmittedCount == 0,
          "composite modulus takes no hint");

    printf("%s: %d failures\n", Failures ? "FAILED" : "passed", Failures);
    return Failures != 0;
}
