/* Exercise the witness state transitions without consuming the tape during point decoding. */
#include <assert.h>
#include <stdio.h>
#include "../src/ec.c"

zkvm_status zkvm_bn254_pairing(const uint8_t *, size_t, bool *);
zkvm_status zkvm_bls12_pairing(const uint8_t *, size_t, bool *);

static void miller(fp12 *f, const point p[4], const point q[4])
{
    fp12_one(f);
    pairing_add_batch(f, p, q, 4);
    pairing_add_batch(f, p, q, 2);
}

static void check_witness(int bls, const uint8_t *pair)
{
    bool empty;
    zkvm_status status = bls ? zkvm_bls12_pairing(NULL, 0, &empty) : zkvm_bn254_pairing(NULL, 0, &empty);
    assert(status == ZKVM_EOK && empty);
    point p[4], q[4];
    assert(pt_load(p, pair, 1, false) && pt_load(q, pair + 8 * C->n, 2, !bls));
    assert(!pt_is_inf(p, 1) && !pt_is_inf(q, 2));
    pt_neg(&p[1], &p[0], 1);
    p[2] = p[0];
    p[3] = p[1];
    q[1] = q[2] = q[3] = q[0];

    /* Three cancelling pairs (P,Q),(-P,Q). Their combined Miller value f belongs to Fp6. */
    fp12 f, baseline, c, cinv, w, power;
    pairing_begin();
    miller(&f, p, q);
    baseline = f;
    assert(pairing_is_one(&baseline) && !pairing_retry());

    /* c = 2 lies in Fp, so all its Frobenius powers equal c and the required w is c^(loop+1)/f. */
    fp12_one(&c);
    fp_add(&c.c[0].c[0], &C->one, &C->one);
    fp12_inv(&cinv, &c);
    fp12_one(&power);
    for (size_t i = 0; i < sizeof C->loop; i++)
        for (int bit = 7; bit >= 0; bit--)
        {
            fp12_mul(&power, &power, &power);
            if (C->loop[i] & (1 << bit))
                fp12_mul(&power, &power, &c);
        }
    fp12_mul(&power, &power, &c);
    fp12_inv(&w, &f);
    fp12_mul(&w, &w, &power);
    for (int i = 3; i < 6; i++)
        assert(fq_is_zero(&w.c[i], 2));
    u256d tape[60];
    const fp12 *parts[] = {&c, &cinv, &w};
    size_t used = 0, witness_start = 0;
    for (int part = 0; part < 3; part++)
    {
        if (part == 2)
            witness_start = used;
        for (int i = 0; i < (part == 2 ? 3 : 6); i++)
            for (int j = 0; j < 2; j++)
                for (uint32_t k = 0; k < (C->n + 7) / 8; k++)
                    memcpy(&tape[used++], &parts[part]->c[i].c[j].l[8 * k], 32);
    }

    /* Force the next miss to start at the cache slot held by the first pair. */
    unsigned slot = 0;
    while (slot < LineCachePoints && !(LineCache[slot].of == C && fq_eq(&LineCache[slot].x, &q[0].x, 2)
                                       && fq_eq(&LineCache[slot].y, &q[0].y, 2)))
        slot++;
    assert(slot < LineCachePoints);
    LineCacheNext = slot;
    point opposite[2] = {q[0], q[0]}, same[2] = {p[0], p[0]};
    pt_neg(&opposite[1], &q[0], 2);
    pairing_begin();
    fp12_one(&f);
    pairing_add_batch(&f, same, opposite, 2);
    assert(pairing_is_one(&f) && !pairing_retry());
    assert(fq_eq(&LineCache[slot].y, &q[0].y, 2));

    for (int retry = 0; retry <= 1; retry++)
    {
        if (retry)
            memset(tape + witness_start, 0, (used - witness_start) * sizeof(*tape)); /* w = 0 */
        FieldHints = (field_hints){tape, tape + used, 1};
        pairing_begin();
        miller(&f, p, q);
        assert(Pw.state == PwMerged && FieldHints.next == FieldHints.end);
        assert(pairing_is_one(&f) == !retry);
        assert(pairing_retry() == retry);
        if (retry)
        {
            miller(&f, p, q);
            assert(pairing_is_one(&f) && !pairing_retry());
        }
        FieldHints = (field_hints){0};
    }
}

int main(void)
{
    char op[8], hex[577];
    int seen = 0;
    while (scanf("%7s %576s", op, hex) == 2)
    {
        int bls = op[1] == 'l', bytes = bls ? 288 : 192;
        uint8_t pair[288];
        for (int i = 0; i < bytes; i++)
        {
            unsigned byte;
            assert(sscanf(hex + 2 * i, "%2x", &byte) == 1);
            pair[i] = (uint8_t)byte;
        }
        check_witness(bls, pair);
        seen |= 1 << bls;
    }
    assert(seen == 3);
    puts("Pairing witnesses on BN254 and BLS12-381: passed");
}
