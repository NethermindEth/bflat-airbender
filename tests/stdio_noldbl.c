/*
 * vfprintf and vfscanf for the harness, in place of musl's: on riscv32, long double is 128-bit and musl's
 * floating-point formatting and scanning need the quad-precision helpers of compiler-rt (__multf3, ...), which it
 * builds only where the compiler has __int128, so not for riscv32. The harness and the linux platform format and
 * scan integers and strings only: %s, %c, %d, %u and %x with a width and a zero flag; %x with a width in scanf.
 * musl's printf, fprintf and sscanf call these.
 */
#include <stdarg.h>
#include <stdio.h>

int vfprintf(FILE *restrict f, const char *restrict fmt, va_list ap)
{
    char buf[256];
    size_t n = 0;
    int total = 0;
#define PUT(c)                              \
    do                                      \
    {                                       \
        if (n == sizeof(buf))               \
        {                                   \
            total += (int)fwrite(buf, 1, n, f); \
            n = 0;                          \
        }                                   \
        buf[n++] = (char)(c);               \
    } while (0)
    for (; *fmt; fmt++)
    {
        if (*fmt != '%')
        {
            PUT(*fmt);
            continue;
        }
        fmt++;
        char pad = ' ';
        if (*fmt == '0')
            pad = '0', fmt++;
        int width = 0;
        while (*fmt >= '0' && *fmt <= '9')
            width = width * 10 + (*fmt++ - '0');
        char digits[16];
        int len = 0;
        switch (*fmt)
        {
        case '%':
            PUT('%');
            continue;
        case 'c':
            PUT(va_arg(ap, int));
            continue;
        case 's':
            for (const char *s = va_arg(ap, const char *); *s; s++)
                PUT(*s);
            continue;
        case 'd':
        case 'u':
        case 'x':
        {
            unsigned v;
            int neg = 0;
            if (*fmt == 'd')
            {
                int d = va_arg(ap, int);
                neg = d < 0;
                v = neg ? 0u - (unsigned)d : (unsigned)d;
            }
            else
                v = va_arg(ap, unsigned);
            unsigned base = *fmt == 'x' ? 16 : 10;
            do
                digits[len++] = "0123456789abcdef"[v % base];
            while (v /= base);
            if (neg)
                width--;
            if (neg && pad == '0')
                PUT('-');
            for (int i = len; i < width; i++)
                PUT(pad);
            if (neg && pad != '0')
                PUT('-');
            while (len)
                PUT(digits[--len]);
            continue;
        }
        default:
            return -1;
        }
    }
    total += (int)fwrite(buf, 1, n, f);
#undef PUT
    return total;
}

int vfscanf(FILE *restrict f, const char *restrict fmt, va_list ap)
{
    int assigned = 0;
    for (; *fmt; fmt++)
    {
        if (*fmt != '%')
        {
            int c = getc(f);
            if (c != (unsigned char)*fmt)
            {
                if (c != EOF)
                    ungetc(c, f);
                return assigned;
            }
            continue;
        }
        fmt++;
        int width = 0;
        while (*fmt >= '0' && *fmt <= '9')
            width = width * 10 + (*fmt++ - '0');
        if (*fmt != 'x')
            return -1;
        unsigned v = 0;
        int got = 0, c;
        while ((width == 0 || got < width) && (c = getc(f)) != EOF)
        {
            int d = c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
            if (d < 0)
            {
                ungetc(c, f);
                break;
            }
            v = v * 16 + (unsigned)d;
            got++;
        }
        if (!got)
            return assigned ? assigned : EOF;
        *va_arg(ap, unsigned *) = v;
        assigned++;
    }
    return assigned;
}
