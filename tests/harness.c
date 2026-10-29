/* Reads one operation per line, "<op> <hex> <hex> ...", and prints its result as hex. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/zkvm.h"

zkvm_status zkvm_keccak256(const uint8_t *data, size_t len, uint8_t output[32]);
zkvm_status zkvm_sha256(const uint8_t *data, size_t len, uint8_t output[32]);
zkvm_status zkvm_ripemd160(const uint8_t *data, size_t len, uint8_t output[32]);
zkvm_status zkvm_blake2f(uint32_t rounds, uint8_t h[64], const uint8_t m[128], const uint8_t t[16], uint8_t f);
zkvm_status zkvm_secp256k1_ecrecover(const uint8_t msg[32], const uint8_t sig[64], uint8_t recid, uint8_t out[64]);
zkvm_status zkvm_secp256r1_verify(const uint8_t msg[32], const uint8_t sig[64], const uint8_t pk[64], bool *ok);
zkvm_status zkvm_modexp(const uint8_t *b, size_t bl, const uint8_t *e, size_t el, const uint8_t *m, size_t ml,
                        uint8_t *out);
zkvm_status zkvm_bn254_g1_add(const uint8_t *p1, const uint8_t *p2, uint8_t *out);
zkvm_status zkvm_bn254_g1_mul(const uint8_t *p, const uint8_t *k, uint8_t *out);
zkvm_status zkvm_bn254_pairing(const uint8_t *pairs, size_t n, bool *ok);
zkvm_status zkvm_kzg_point_eval(const uint8_t *c, const uint8_t *z, const uint8_t *y, const uint8_t *proof, bool *ok);
zkvm_status zkvm_bls12_g1_add(const uint8_t *p1, const uint8_t *p2, uint8_t *out);
zkvm_status zkvm_bls12_g1_msm(const uint8_t *pairs, size_t n, uint8_t *out);
zkvm_status zkvm_bls12_g2_add(const uint8_t *p1, const uint8_t *p2, uint8_t *out);
zkvm_status zkvm_bls12_g2_msm(const uint8_t *pairs, size_t n, uint8_t *out);
zkvm_status zkvm_bls12_pairing(const uint8_t *pairs, size_t n, bool *ok);
zkvm_status zkvm_bls12_map_fp_to_g1(const uint8_t *fe, uint8_t *out);
zkvm_status zkvm_bls12_map_fp2_to_g2(const uint8_t *fe, uint8_t *out);

static uint8_t Args[6][65536];
static size_t ArgLen[6];

static size_t unhex(uint8_t *out, const char *s)
{
    size_t n = 0;
    if (strcmp(s, "-") == 0)
        return 0;
    for (; s[0] && s[1]; s += 2)
    {
        unsigned v;
        sscanf(s, "%2x", &v);
        out[n++] = (uint8_t)v;
    }
    return n;
}

static void print_hex(const uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; i++)
        printf("%02x", p[i]);
}

int main(void)
{
    static char line[300000];
    while (fgets(line, sizeof(line), stdin))
    {
        char *op = strtok(line, " \n");
        if (!op)
            continue;
        int argc = 0;
        for (char *tok; argc < 6 && (tok = strtok(NULL, " \n"));)
        {
            ArgLen[argc] = unhex(Args[argc], tok);
            argc++;
        }
        uint8_t out[4096];
        zkvm_status st;
        if (!strcmp(op, "keccak"))
        {
            st = zkvm_keccak256(Args[0], ArgLen[0], out);
            print_hex(out, 32);
        }
        else if (!strcmp(op, "sha256"))
        {
            st = zkvm_sha256(Args[0], ArgLen[0], out);
            print_hex(out, 32);
        }
        else if (!strcmp(op, "ripemd"))
        {
            st = zkvm_ripemd160(Args[0], ArgLen[0], out);
            print_hex(out, 32);
        }
        else if (!strcmp(op, "blake2f"))
        {
            /* EIP-152 input: rounds(4 BE) h(64) m(128) t(16) f(1) */
            uint8_t *in = Args[0];
            uint32_t rounds = ((uint32_t)in[0] << 24) | (in[1] << 16) | (in[2] << 8) | in[3];
            memcpy(out, in + 4, 64);
            st = zkvm_blake2f(rounds, out, in + 68, in + 196, in[212]);
            print_hex(out, 64);
        }
        else if (!strcmp(op, "ecrecover"))
        {
            st = zkvm_secp256k1_ecrecover(Args[0], Args[1], Args[2][0], out);
            if (st == ZKVM_EOK)
                print_hex(out, 64);
            else
                printf("fail");
        }
        else if (!strcmp(op, "ecrecoverbatch"))
        {
            /* Hashes, signatures and recovery ids, concatenated; each key or "fail", comma-separated. */
            static uint8_t keys[64 * 1024], ok[1024];
            uint32_t n = (uint32_t)ArgLen[2];
            airbender_secp256k1_ecrecover_batch(Args[0], Args[1], Args[2], keys, ok, n);
            for (uint32_t i = 0; i < n; i++)
            {
                if (i)
                    printf(",");
                if (ok[i])
                    print_hex(keys + 64 * i, 64);
                else
                    printf("fail");
            }
        }
        else if (!strcmp(op, "verifybatch"))
        {
            /* Hashes, signatures, recovery ids, keys, R's ys and weights, each concatenated. */
            printf("%d", airbender_secp256k1_verify_batch(Args[0], Args[1], Args[2], Args[3], Args[4], Args[5],
                                                          (uint32_t)ArgLen[2]));
        }
        else if (!strcmp(op, "p256"))
        {
            bool ok;
            st = zkvm_secp256r1_verify(Args[0], Args[1], Args[2], &ok);
            printf("%d", ok);
        }
        else if (!strcmp(op, "modexp"))
        {
            st = zkvm_modexp(Args[0], ArgLen[0], Args[1], ArgLen[1], Args[2], ArgLen[2], out);
            if (ArgLen[2] == 0)
                printf("-");
            print_hex(out, ArgLen[2]);
        }
        else if (!strcmp(op, "bnadd") || !strcmp(op, "bnmul") || !strcmp(op, "blsg1add") || !strcmp(op, "blsg2add") ||
                 !strcmp(op, "blsg1msm") || !strcmp(op, "blsg2msm") || !strcmp(op, "blsmapg1") ||
                 !strcmp(op, "blsmapg2"))
        {
            /* One operand holding the whole input, split here as the precompile does. */
            size_t n = 64;
            if (!strcmp(op, "bnadd"))
                st = zkvm_bn254_g1_add(Args[0], Args[0] + 64, out);
            else if (!strcmp(op, "bnmul"))
                st = zkvm_bn254_g1_mul(Args[0], Args[0] + 64, out);
            else if ((n = 96, !strcmp(op, "blsg1add")))
                st = zkvm_bls12_g1_add(Args[0], Args[0] + 96, out);
            else if (!strcmp(op, "blsg1msm"))
                st = zkvm_bls12_g1_msm(Args[0], ArgLen[0] / 128, out);
            else if (!strcmp(op, "blsmapg1"))
                st = zkvm_bls12_map_fp_to_g1(Args[0], out);
            else if ((n = 192, !strcmp(op, "blsg2add")))
                st = zkvm_bls12_g2_add(Args[0], Args[0] + 192, out);
            else if (!strcmp(op, "blsg2msm"))
                st = zkvm_bls12_g2_msm(Args[0], ArgLen[0] / 224, out);
            else
                st = zkvm_bls12_map_fp2_to_g2(Args[0], out);
            if (st == ZKVM_EOK)
                print_hex(out, n);
            else
                printf("fail");
        }
        else if (!strcmp(op, "bnpair") || !strcmp(op, "blspair") || !strcmp(op, "kzg"))
        {
            bool ok = false;
            if (!strcmp(op, "bnpair"))
                st = zkvm_bn254_pairing(Args[0], ArgLen[0] / 192, &ok);
            else if (!strcmp(op, "blspair"))
                st = zkvm_bls12_pairing(Args[0], ArgLen[0] / 288, &ok);
            else
                st = zkvm_kzg_point_eval(Args[0], Args[1], Args[2], Args[3], &ok);
            if (st == ZKVM_EOK)
                printf("%d", ok);
            else
                printf("fail");
        }
        else
        {
            printf("unknown");
            st = ZKVM_EFAIL;
        }
        (void)st;
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
