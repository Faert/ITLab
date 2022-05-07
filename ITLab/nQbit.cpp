#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <chrono>
#include "Qbit.h"
#include "RM_RSA.h"

using namespace std;

int main(int argc, char* argv[])
{
    auto start = std::chrono::steady_clock::now();

    //Shor test
    size_t n = 4;

    size_t a = 2;
    size_t N = 15;

    double error = 0;//in %, if in bar chart then int

    cout << "P, CP, CPP error = " << error << "%\n";

    Qbit<double> qb(4*n + 3);

    qb.Shor(a, N, 0, 4 * n + 3, error);

    qb.condition_exp_cout(0, 2 * n, (1i64 << (n * 2 + 3)));
    
    size_t res = 1;

    //qb.cleaning_up_small_errors();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n\n";

    ofstream out;
    out.open("Shor_result.txt");

    if (out.is_open())
    {
        //qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);
        vector<size_t> temp = qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1] * 2)) && (temp[i] > (temp[i - 1] * 2)))
            {
                res++;
            }
        }

        cout << res << '\n';

        out << a << ' ' << N << ' ' << error << ' ' << res << '\n';
    }

    out.close();


    //sample RSA
    cout << RSA_decryption(RSA_encryption(15, 10)) << '\n' << '\n';

    //RSA hacking
    
    //RSA open key{e, pq}
    size_t pq = 15;
    size_t e = 7;

    //RSA encryption
    size_t m = 3;//message
    size_t c = ModExp(m, e, pq);

    Qbit<double> SHOR_RSA(19);

    size_t res1 = 1;

    while (true)
    {
        SHOR_RSA.clear();
        SHOR_RSA.Shor(e, pq, 0, 19, 1);

        vector<size_t> temp1 = SHOR_RSA.condition_exp(0, 8, 1i64 << 11);

        for (size_t i = 1; i < temp1.size(); i++)
        {
            if ((temp1[i] > (temp1[0] / 2)) && (temp1[i] > (temp1[i + 1] * 2)) && (temp1[i] > (temp1[i - 1] * 2)))
            {
                res1++;
            }
        }

        if (res1 & 1)
        {
            while (gcd(e, pq) != 1)
            {
                e++;
                if (e >= pq)
                {
                    e = 2;
                }
            }
        }
        else
        {
            break;
        }
    }

    //if res1%2 = 0, else start shor(a, pq) with random a: gcd(a, pq) = 1 
    size_t p = gcd(MyPow(e, res1 / 2) + 1, pq);
    size_t q = pq / p;
    size_t fn = (p - 1) * (q - 1);

    size_t d = inverse_element_by_mod(e, fn);

    cout << m << ' ' <<  ModExp(c, d, pq) << '\n';

    //Sample dense coding
    Qbit<double> Dense(4); //0, 1 Alice; 2, 3 Bob in 0
    size_t st = 1;//st < 2^2
    Dense[st] = 1;
    
    Dense.H(2);
    Dense.CNOT(2, 3);

    Dense.CNOT(1, 2);
    Dense.H(2);
    Dense.CNOT(0, 2);
    Dense.H(2);

    Dense.CNOT(2, 3);
    Dense.H(2);

    Dense.condition_exp_cout(0, 4, 100);//Alice = Bob -> res = st + (st << 2)

    //system("pause");
}
