/* Host input stream for the Airbender reader. Delegations use their portable implementations. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int allocations_until_failure = -1;
static void *input_malloc(size_t size)
{
    if (allocations_until_failure == 0)
        return NULL;
    if (allocations_until_failure > 0)
        allocations_until_failure--;
    return malloc(size);
}

#define malloc input_malloc
#include "../../src/platform_airbender.c"
#undef malloc

field_hints FieldHints;
static uint32_t stream[96], cursor, length;

uint32_t csr_read(void)
{
    assert(cursor <= length);
    return cursor++ ? stream[cursor - 2] : length * 4;
}

void airbender_finish_error(uint32_t code)
{
    fprintf(stderr, "unexpected guest exit %u\n", code);
    abort();
}

static void check_read(int prehash)
{
    const uint8_t *data;
    size_t size;
    cursor = 0;
    if (prehash)
        airbender_read_prehashed_input(&data, &size);
    else
        read_input(&data, &size);
    assert(size == length * 4 && cursor == length + 1);
    assert(memcmp(data, stream, size) == 0);
}

int main(void)
{
    stream[0] = TrustedMagic;
    stream[1] = 1;
    uint32_t i = TrustedHeaderWords;
    stream[TrustedNodesField] = i * 4;
    stream[i] = 1;
    i += 2;
    stream[i] = 1;
    i += TrustedItemPrefixWords;
    stream[i] = 0x80;
    i += 2;
    stream[TrustedCodesField] = i * 4;
    i += 2; /* A valid empty code list. */
    stream[TrustedTotalField] = i * 4;
    stream[i] = 1;
    i += 2;
    uint32_t div = i;
    for (int j = 0; j < 16; j++)
        stream[i++] = j + 1;
    stream[i] = 1;
    i += 2;
    uint32_t field = i;
    for (int j = 0; j < 8; j++)
        stream[i++] = j + 17;
    length = i;

    check_read(0);
    assert(!Prehashed.node_count && !FieldHints.active);
    allocations_until_failure = 1; /* The raw buffer fits but the retained states do not. */
    check_read(1);
    assert(!Prehashed.node_count && !FieldHints.active);
    allocations_until_failure = -1;
    check_read(1);
    assert(Prehashed.node_count == 1 && Prehashed.code_count == 0);
    assert(DivisionHints.count == 1 && ((uintptr_t)DivisionHints.start & 31) == 0);
    assert(memcmp(DivisionHints.start, stream + div, 64) == 0);
    assert(FieldHints.active && FieldHints.end - FieldHints.next == 1);
    assert(((uintptr_t)FieldHints.next & 31) == 0 && memcmp(FieldHints.next, stream + field, 32) == 0);

    uint8_t address[20] = {1}, bloom[256] = {0}, expected[256] = {0}, digest[32];
    keccak256_long(address, sizeof(address), digest);
    for (int j = 0; j < 6; j += 2)
    {
        unsigned bit = ((digest[j] & 7u) << 8) | digest[j + 1];
        expected[255 - bit / 8] |= 1u << (bit % 8);
    }
    airbender_bloom_set(bloom, NULL, address, sizeof(address));
    assert(memcmp(bloom, expected, sizeof(bloom)) == 0);
    puts("Airbender input, hints and bloom: passed");
}
