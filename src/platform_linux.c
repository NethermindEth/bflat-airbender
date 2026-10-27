/*
 * The I/O half for running a guest as a Linux process (qemu-riscv32, or natively for tests):
 * the input is stdin, the public output is printed as hex, and a failure is an exit status.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "zkvm.h"

void read_input(const uint8_t **buf_ptr, size_t *buf_size)
{
    size_t cap = 1 << 20, len = 0;
    uint8_t *buf = malloc(cap);
    for (;;)
    {
        if (len == cap)
            buf = realloc(buf, cap *= 2);
        ssize_t got = read(0, buf + len, cap - len);
        if (got <= 0)
            break;
        len += (size_t)got;
    }
    *buf_ptr = buf;
    *buf_size = len;
}

void write_output(const uint8_t *output, size_t size)
{
    printf("output: ");
    for (size_t i = 0; i < size; i++)
        printf("%02x", output[i]);
    printf("\n");
    fflush(stdout);
}

void sys_write(uint32_t fd, const uint8_t *ptr, size_t len)
{
    (void)write((int)fd, ptr, len);
}

void airbender_fail(uint32_t code)
{
    fprintf(stderr, "airbender_fail 0x%x\n", code);
    _exit(100);
}
