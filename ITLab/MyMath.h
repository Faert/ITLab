#pragma once

//by Faert

template<typename T>
T MyPow(T a, size_t pow)
{
    T result = T(1);
    while (pow != 0)
    {
        if (pow & 1)
        {
            result *= a;
        }
        a *= a;
        pow >>= 1;
    }
    return result;
}

unsigned long long int ModMul(unsigned long long int a, unsigned long long int b, unsigned long long int M)
{
    if (b == 0)
    {
        return 0;
    }

    unsigned long long int z = ModMul(a, b / 2, M);

    if (b % 2 == 0)
    {
        return ((2 * z) % M);
    }
    else
    {
        return ((a + ((2 * z) % M)) % M);
    }
}

unsigned long long int ModExp(unsigned long long int x, unsigned long long int n, unsigned long long int M)
{
    if (n == 0)
    {
        return 1;
    }

    unsigned long long int z = ModExp(x, n / 2, M);

    if (n % 2 == 0)
    {
        return ModMul(z, z, M);
    }
    else
    {
        return ModMul(ModMul(z, z, M), x, M);
    }
}

unsigned long long int gcd(unsigned long long int a, unsigned long long int b)
{
    return b ? gcd(b, a % b) : a;
}


unsigned long long int extended_gcd(unsigned long long int a, unsigned long long int b, long long int&x, long long int&y)
// a >= b >= 0
{
    if (b == 0)
    {
        x = 1; y = 0;
        return a;
    }
    long long int x1 = 0, y1 = 0;
    unsigned long long int d = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a/b)*y1;
    return d;
}

unsigned long long int inverse_element_by_mod(unsigned long long int a, unsigned long long int N)
//it may be that the inverse element does not exist(gcd(a, N) != 1)
{
    long long int x = 0, y = 0;
    extended_gcd(a, N, x, y);
    while (x < 0)
    {
        x += N;
    }
    return x;
}