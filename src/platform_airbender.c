/*
 * The I/O half for the Airbender zkVM. Everything goes through the non-determinism CSR 0x7c0:
 *   - a read pops the next input word; the input is one word of byte length, then the bytes
 *     packed little-endian into words;
 *   - a write feeds the simulator's quasi-UART: 0xffffffff, the number of words that follow,
 *     the byte length, then the bytes packed little-endian into words. Debug output only; the
 *     prover ignores it.
 * The public output is the eight words left in x10-x17 at the end of the run, which
 * airbender_finish_success (platform_airbender.S) loads in the exact sequence the Airbender
 * tooling looks for in the binary.
 */
#include "zkvm.h"

void *malloc(size_t size);
__attribute__((noreturn)) void airbender_finish_success(const uint32_t words[16]);
__attribute__((noreturn)) void airbender_finish_error(uint32_t code);

/* x10-x17 the output, x18-x25 zero, as the recursion convention expects. */
static uint32_t Output[16];

static inline uint32_t csr_read(void)
{
    uint32_t value;
    __asm__ volatile("csrrw %0, 0x7c0, x0" : "=r"(value));
    return value;
}

static inline void csr_write(uint32_t value)
{
    __asm__ volatile("csrrw x0, 0x7c0, %0" : : "r"(value));
}

void read_input(const uint8_t **buf_ptr, size_t *buf_size)
{
    uint32_t len = csr_read();
    uint32_t *words = malloc(((size_t)len + 3) & ~(size_t)3);
    for (uint32_t i = 0; i < (len + 3) / 4; i++)
        words[i] = csr_read();
    *buf_ptr = (const uint8_t *)words;
    *buf_size = len;
}

void write_output(const uint8_t *output, size_t size)
{
    if (size > 32)
        airbender_fail(0x201);
    for (size_t i = 0; i < 8; i++)
    {
        uint32_t w = 0;
        for (size_t b = 0; b < 4; b++)
            if (4 * i + b < size)
                w |= (uint32_t)output[4 * i + b] << (8 * b);
        Output[i] = w;
    }
}

void sys_write(uint32_t fd, const uint8_t *ptr, size_t len)
{
    (void)fd;
    csr_write(0xffffffff);
    csr_write((uint32_t)(1 + (len + 3) / 4));
    csr_write((uint32_t)len);
    for (size_t i = 0; i < len; i += 4)
    {
        uint32_t w = 0;
        for (size_t b = 0; b < 4 && i + b < len; b++)
            w |= (uint32_t)ptr[i + b] << (8 * b);
        csr_write(w);
    }
}

/* The end of a run: success publishes Output, anything else proves nothing. */
void airbender_exit(int status)
{
    if (status == 0)
        airbender_finish_success(Output);
    airbender_finish_error((uint32_t)status);
}

void airbender_fail(uint32_t code)
{
    airbender_finish_error(code);
}
