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
    size_t n = 5;

    size_t a = 30;
    size_t N = 31;

    double error = 5;//in %, if in bar chart then int

    cout << "P, CP, CPP error = " << error << "%\n";

    Qbit<double> qb(4*n + 2);

    qb.Shor(a, N, 0, 4 * n + 2, error);

    qb.condition_exp_cout(0, 2 * n, (1i64 << (n * 2 + 3)));//
    
    size_t res = 1;

    //qb.cleaning_up_small_errors();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n\n";

    ofstream out;
    out.open("Shor_result.txt");

    if (out.is_open())
    {
        vector<size_t> temp = qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1] * 2)) && (temp[i] > (temp[i - 1] * 2)))
            {
                res++;
            }
        }

        cout << res << '\n' << '\n';

        out << a << ' ' << N << ' ' << error << ' ' << res << ' ' << round(elapsed_seconds.count()) << '\n';
    }

    out.close();


    //sample RSA
    cout << RSA_decryption(RSA_encryption(15, 10)) << '\n' << '\n';


    auto start1 = std::chrono::steady_clock::now();
    //RSA hacking
    
    //RSA open key{e, pq}
    size_t pq = 10;
    size_t e = 3;

    size_t n1 = 4;

    //RSA encryption
    size_t m = 3;//message
    size_t c = ModExp(m, e, pq);

    Qbit<double> SHOR_RSA(4*n1 + 2);

    size_t res1 = 1;

    while (true)
    {
        SHOR_RSA.clear();
        SHOR_RSA.Shor(e, pq, 0, 4 * n1 + 2, 1);

        vector<size_t> temp1 = SHOR_RSA.condition_exp(0, 2 * n1, 1i64 << 2 * n1 + 3);

        for (size_t i = 1; i < temp1.size(); i++)
        {
            if ((temp1[i] > (temp1[0] / 2)) && (temp1[i] > (temp1[i + 1] * 2)) && (temp1[i] > (temp1[i - 1] * 2)))
            {
                res1++;
            }
        }

        size_t e_ = e;

        size_t tempp = (gcd(MyPow(e, res1 / 2) + 1, pq));
        size_t tempm = (gcd(MyPow(e, res1 / 2) - 1, pq));
        if ((res1 & 1) || ((tempp == pq || tempp == 1) && (tempm == pq || tempm == 1)))
        {
            e++;
            while (gcd(e, pq) != 1)
            {
                e++;
                if (e >= pq)
                {
                    e = 2;
                }
                
            }

            if (e = e_)
            {
                cout << "Wrong pq\n";
                break;
            }

            res1 = 1;
        }
        else
        {
            break;
        }
    }

    //if res1%2 = 0, else start shor(a, pq) with random a: gcd(a, pq) = 1 
    size_t p = gcd(MyPow(e, res1 / 2) + 1, pq);
    if (p == pq)
    {
        p = gcd(MyPow(e, res1 / 2) - 1, pq);
    }
    size_t q = pq / p;
    size_t fn = (p - 1) * (q - 1);

    size_t d = inverse_element_by_mod(e, fn);

    cout << m << ' ' <<  ModExp(c, d, pq) << '\n';

    auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
    std::cout << "elapsed time: " << elapsed_seconds1.count() << "s\n\n";

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
