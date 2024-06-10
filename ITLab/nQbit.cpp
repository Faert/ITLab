#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <chrono>
#include "Qbit.h"
#include "RM_RSA.h"

using namespace std;

size_t Shor(size_t a, size_t N, size_t n, double error = 0, bool print = false, string out_file = "Shor_result.txt")
{
    if (print)
    {
        cout << "P, CP, CPP error = " << error << "%\n";
    }

    auto start = std::chrono::steady_clock::now();

    Qbit<double> qb(4 * n + 2);

    qb.Shor(a, N, 0, 4 * n + 2, error);

    if (print)
    {
        qb.condition_exp_cout(0, 2 * n, (1i64 << (n * 2 + 3)));
    }

    size_t res = 1;

    //qb.cleaning_up_small_errors();

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    if (print)
    {
        std::cout << "Shor " << a << ' ' << N << " time: " << elapsed_seconds.count() << "s\n";
    }

    ofstream out;
    out.open(out_file);

    if (out.is_open())
    {
        vector<size_t> temp = qb.condition_exp_in_file(0, 2 * n, (1i64 << (n * 2 + 3)), out);

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1])) && (temp[i] > (temp[i - 1])))
            {
                res++;
            }
        }

        if (print)
        {
            cout << res << "\n\n";
        }

        out << a << ' ' << N << ' ' << error << ' ' << res << ' ' << round(elapsed_seconds.count()) << '\n';
    }
    else
    {
        vector<size_t> temp = qb.condition_exp(0, 2 * n, (1i64 << (n * 2 + 3)));

        for (size_t i = 1; i < temp.size(); i++)
        {
            if ((temp[i] > (temp[0] / 2)) && (temp[i] > (temp[i + 1])) && (temp[i] > (temp[i - 1])))
            {
                res++;
            }
        }

        if (print)
        {
            cout << res << "\n\n";
        }
    }

    out.close();

    return res;
}

size_t RSA_hacking(size_t e, size_t pq, size_t n1, size_t c)
{
    auto start = std::chrono::steady_clock::now();

    size_t res1 = 1;
    size_t e_ = e;

    while (true)
    {
        res1 = Shor(e, pq, n1, 0, false, "RSA_Shor_result.txt");

        size_t tempp = (gcd(MyPow(e, res1 / 2) + 1, pq));
        size_t tempm = (gcd(MyPow(e, res1 / 2) - 1, pq));
        if ((res1 & 1) || ((tempp == pq && tempp == 1) && (tempm == pq && tempm == 1)))
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

            if (e == e_)
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

    size_t p = gcd(MyPow(e, res1 / 2) + 1, pq);
    if (p == pq || p == 1)
    {
        p = gcd(MyPow(e, res1 / 2) - 1, pq);
    }
    size_t q = pq / p;
    size_t fn = (p - 1) * (q - 1);

    size_t d = inverse_element_by_mod(e_, fn);

    size_t m_ = ModExp(c, d, pq);

    cout << m_ << '\n';

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "RSA hacking time: " << elapsed_seconds.count() << "s\n\n";

    return (m_);
}

void Dense_coding(Qbit<double> Dense, size_t count_exp = 100)//Dense 2,3 qbits = 0
{
    Dense.H(2);
    Dense.CNOT(2, 3);

    Dense.CNOT(1, 2);
    //Z
    Dense.H(2);
    Dense.CNOT(0, 2);
    Dense.H(2);

    Dense.CNOT(2, 3);
    Dense.H(2);

    Dense.condition_exp_cout(2, 4, count_exp);//Alice = Bob -> res = st
}

void Bells_inequality(Qbit<double> Bell, double a = 0, double b = 0, size_t count_exp = 100)//Bell(2 qubit) in 11
{
    Bell.H(1);
    Bell.CNOT(1, 0);

    Bell.RZ(0, a);
    Bell.RZ(1, b);

    Bell.H(0);
    Bell.H(1);

    for(int i =0; i < 4; i++)
    cout << Bell[i] << '\n';

    cout << '\n';

    Bell.condition_exp_cout(0, 2, count_exp);
}

int main(int argc, char* argv[])
{
    //Shor
    
    size_t n = 4;

    size_t a = 11;//gcd(a, N) == 1, a < N
    size_t N = 15;//N < 2^n

    double error = 0;//in %, if in bar chart then int

    Shor(a, N, n, error, true);
    

    /*
    //sample RSA
    size_t mes = 15;//message
    size_t bit = 10;

    cout << RSA_decryption(RSA_encryption(mes, bit)) << '\n' << '\n';

    //RSA hacking
    
    //RSA open key{e, pq}
    size_t pq = 14;
    size_t e = 5;//gcd(e, fn) = 1, fn = (p-1)(q-1) 

    size_t n1 = 4;

    //RSA encryption
    size_t m = 3;//message
    size_t c = ModExp(m, e, pq);

    RSA_hacking(e, pq, n1, c);


    //Sample dense coding
    Dense_coding(1);

    //system("pause");
    */

    /*
    //Dense coding
    Qbit<double> Dense(4); //0, 1 Alice; 2, 3 Bob in 0

    Dense[1] = im;
    Dense[2] = 3;
    Dense.normalization();

    cout << "Dense coding:\n";
    Dense_coding(Dense, 1000);
    cout << "\n\n";

    //Quantum Fourier transform
    cout << "Quantum Fourier transform\n";

    Qbit<double> QFT(3);

    QFT[3] = 1;

    QFT.condition_exp_cout(0, 3, 1000);

    cout << '\n';

    QFT.QFT(0, 3);

    for (size_t i = 0; i < 8; i++)
        cout << i << ' ' << QFT[i] << '\n';
    cout << '\n';

    cout << "Second QFT error 38.5 \n";
    QFT.QFT(0, 3, 38.5);

    for (size_t i = 0; i < 8; i++)
        cout << i << ' ' << QFT[i] << '\n';
    cout << '\n';


    QFT.condition_exp_cout(0, 3, 1000);
    cout << '\n';

    //size_t a = 2;
    //size_t start = 0;
    //size_t end = 3;
    //double error = 0;

    //сумматор с обратным к a
    /*
    for (size_t j = end-1; (j >= start) && (j < end); j--)
    {
        if ((a >> (j - start)) & 1)
        {
            for (size_t i = j; (i >= start) && (i <= j); i--)
            {
                QFT.P(i, PI / (1 << (j - i)), error);
            }
        }
    }
    */

    /*
    //сумматор
    for (size_t j = end - 1; (j >= start) && (j < end); j--)
    {
        if ((a >> ((end - 1) - j)) & 1)
        {
            for (size_t i = j; (i >= start) && (i <= j); i--)
            {
                QFT.P(i, PI / (1 << (j - i)), error);
            }
        }
    }
    */

    //QFT.RQFT(0, 3);
    
    //QFT.condition_exp_cout(0, 3, 1000);



    //Bells inequality
    
    /*
    cout << "\nBells inequality\n";

    Qbit<double> Bell(2);
    Bell[3] = 1;

    double a = PI / 3;
    double b = 0;

    Bells_inequality(Bell, a, b, 8192);
    */
}
