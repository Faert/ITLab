#pragma once

#include <iostream>
#include <chrono>
#include <random>
#include <array>
#include "MyMath.h"

using namespace std;

void MillerRabin(unsigned long long int N, vector<unsigned long long int>& vec)
{
    mt19937 gen(N);
    uniform_int_distribution<unsigned long long int> uid(2, N - 1);

    while (N != 1)
    {
        if (N % 2 == 0)
        {
            vec.push_back(2);
            N /= 2;
            continue;
        }

        unsigned long long int X = 0;
        uid = uniform_int_distribution<unsigned long long int>(2, N - 1);
        X = uid(gen);
        if (gcd(X, N) != 1)
        {
            vec.push_back(gcd(X, N));
            N /= vec.back();
            continue;
        }

        unsigned long long int u = N - 1;
        unsigned long long int t = 0;

        while (u % 2 == 0) //N - 1 = (2^t)*u
        {
            u /= 2;
            t++;
        }

        unsigned long long int temp = u;

        for (unsigned long long int i = 0; i <= t; i++)
        {
            if (ModExp(X, u, N) == 1)
            {
                if (i == t)
                {
                    vec.push_back(N);
                    N /= N;
                    break;
                }

                if (u != N - 1 && i != 0 && ModExp(temp, u, N) != N - 1 && N % temp == 0)
                {
                    vec.push_back(temp);
                    N /= temp;
                }

                if (i != 0)
                {
                    break;
                }
            }

            temp = u;
            u *= 2;
        }
    }
}

array<unsigned long long int, 3> RSAkey(size_t bit = 16) //m < 2^(2*bit - 2) && bit < 32
{
    mt19937 gen(bit);
    uniform_int_distribution<unsigned long long int> uid((1i64 << (bit - 1)), (1i64 << (bit)));
    unsigned long long int p = uid(gen);
    unsigned long long int q = uid(gen);

    if (p % 2 == 0)
    {
        p++;
    }
    if (q % 2 == 0)
    {
        q++;
    }

    vector<unsigned long long int> v{};
    while (v.size() != 1)
    {
        v.clear();
        MillerRabin(p, v);
        p += 2;
    }
    v.clear();
    while (v.size() != 1)
    {
        v.clear();
        MillerRabin(q, v);
        q += 2;
    }
    p -= 2;
    q -= 2;

    unsigned long long int fn = (p - 1) * (q - 1);
    unsigned long long int n = p * q;

    uid = uniform_int_distribution<unsigned long long int>(2, fn - 1);
    unsigned long long int e = uid(gen);
    while (gcd(e, fn) != 1)
    {
        e++;
    }

    unsigned long long int d = inverse_element_by_mod(e, fn);
    return array<unsigned long long, 3>{ e, d, n }; 
}

array<unsigned long long int, 3> RSA_encryption(unsigned long long int m, unsigned long long int bit = 16) // out{c, d, n}
{
    array<unsigned long long int, 3> a = RSAkey(bit);
    unsigned long long int c = ModExp(m, a[0], a[2]);
    a[0] = c;
    return a;
}

unsigned long long int RSA_decryption(array<unsigned long long int, 3> a) // a{c, d, n}
{
    return ModExp(a[0], a[1], a[2]);
}