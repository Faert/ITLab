#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <omp.h>//in the properties
#include "Qbit.h"
#include "RM_RSA.h"
//#include "OLD.h"

using namespace std;

int main(int argc, char* argv[])
{
    //argc and argv when final release

    auto start = std::chrono::steady_clock::now();

    //Shor test
    size_t n = 4;

    size_t a = 7;
    size_t N = 15;

    double error = 0;//in %, if in hist then int

    cout << "P, CP, CPP error = " << error << "%\n";

    Qbit<double> qb(4*n + 3);

    qb.Shor(a, N, 0, 4 * n + 3, error);

    //Instead of qb.condition_exp(0, 2*n, (1i64 << (n * 2 + 3))) because we count the res
    
    size_t res = 1;

    vector<size_t> temp(1i64 << n * 2);

    for (size_t i = 0; i < (1i64 << (n * 2 + 3)); i++)
    {
        double r = double(rand()) / RAND_MAX;
        double sum = 0;
        size_t k = 0;
        while (sum < r)
        {
            sum += norm(qb[k]);
            k++;
        }
        temp[(k-1) % (1i64 << n * 2)]++;
    }

    for (size_t i = 0; i < temp.size(); i++)
    {
        if (temp[i] != 0)
        {
            cout << '(' << i << "; " << temp[i] << ")\n";
        }
    }

    for (size_t i = 1; i < temp.size(); i++)
    {
        if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1] * 2)) && (temp[i] > (temp[i - 1] * 2)))
        {
            res++;
        }
    }

    cout << res << '\n';

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n\n";

    ofstream out;
    out.open("Shor_result.txt");

    if (out.is_open())
    {
        qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);
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
    SHOR_RSA.Shor(e, pq, 0, 19, 2);

    res = 1;

    vector<size_t> temp1(1i64 << 8);

    for (size_t i = 0; i < (1i64 << 11); i++)
    {
        double r = double(rand()) / RAND_MAX;
        double sum = 0;
        size_t k = 0;
        while (sum < r)
        {
            sum += norm(qb[k]);
            k++;
        }
        temp1[(k - 1) % (1i64 << 8)]++;
    }

    for (size_t i = 1; i < temp1.size(); i++)
    {
        if ((temp1[i] > (temp1[0] / 2)) && (temp1[i] > (temp1[i + 1] * 2)) && (temp1[i] > (temp1[i - 1] * 2)))
        {
            res++;
        }
    }

    //in file for bar chart
    ofstream outRSA;
    out.open("Shor_RSA_result.txt");

    if (outRSA.is_open())
    {
        qb.condition_exp_in_file(0, 8, (1i64 << 11), outRSA);
        outRSA << a << ' ' << N << ' ' << error << ' ' << res << '\n';
    }

    outRSA.close();

    //if res%2 = 0, else start shor(a, pq) with random a: gcd(a, pq) = 1 
    size_t p = gcd(MyPow(e, res / 2) + 1, pq);
    size_t q = pq / p;
    size_t fn = (p - 1) * (q - 1);

    size_t d = inverse_element_by_mod(e, fn);

    cout << m << ' ' <<  ModExp(c, d, pq) << '\n';

    //Sample dense coding
    Qbit<double> Dense(4); //0, 1 Alice, 2, 3 Bob
    size_t st = 1;
    Dense[st] = 1;
    
    Dense.H(2);
    Dense.CNOT(2, 3);

    Dense.CNOT(1, 2);
    Dense.H(2);
    Dense.CNOT(0, 2);
    Dense.H(2);

    Dense.CNOT(2, 3);
    Dense.H(2);

    Dense.condition_exp(0, 4, 100);
    //
    //system("pause");
}
