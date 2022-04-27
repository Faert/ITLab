#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "MyMath.h"

using namespace std;

const double PI = 3.141592653589793;// 238463;
complex<double> im(0, 1);

template<typename T>
class Qbit
{
private:
	vector<complex<T>> data;
	size_t size;

	const size_t get_size() const
	{
		return data.size();
	}

public:
	Qbit(size_t s = 4): data(vector<complex<T>>(1i64 << s)), size(s) {} //max s = 31
	Qbit(const Qbit& other): data(vector<complex<T>>(other.data)), size(other.size) {}
	Qbit(const vector<complex<T>>& other): data(vector<complex<T>>(other)) 
	{
		size_t DataSize = data.size();
		while (DataSize != 0)
		{
			DataSize /= 2;
			size++;
		}
		size--;
	}

	const size_t qsize() const
	{
		return size;
	}

	complex<T>& operator[](size_t i)
	{
		return data.at(i);
	}

	const complex<T>& operator[](size_t i) const
	{
		return data.at(i);
	}

	void normalization()
	{
		complex<T> len = 0;

		for (size_t j = 0; j < data.size(); j++)
		{
			len += norm((*this)[j]);
		}

		len = sqrt(len);

		if (len.real() != 1)
		{
			for (size_t j = 0; j < data.size(); j++)
			{
				(*this)[j] /= len;
			}
		}
	}

	void H(size_t n)
	{
		size_t step = 1 << n;

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (!((i % (1 << (n + 1))) >> n))
			{
				complex<T> temp = (*this)[i];
				(*this)[i] = (temp + (*this)[i + step]) / sqrt(2);
				(*this)[i + step] = (temp - (*this)[i + step]) / sqrt(2);
			}
		}
	}

	void X(size_t n) // NOT
	{
		size_t step = 1 << n;

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (!((i % (1 << (n + 1))) >> n))
			{
				swap((*this)[i], (*this)[i + step]);
			}
		}
	}

	void CNOT(size_t h, size_t l)//h - control
	{
		size_t lstep = 1 << l;

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (((i % (1 << (h + 1))) >> h) && !((i % (1 << (l + 1))) >> l))
			{
				swap((*this)[i], (*this)[i + lstep]);
			}
		}
	}

	void CCNOT(size_t h1, size_t h2, size_t l)//h1, h2 - control
	{
		size_t lstep = 1 << l;

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (((i % (1 << (h1 + 1))) >> h1) && ((i % (1 << (h2 + 1))) >> h2) && !((i % (1 << (l + 1))) >> l))
			{
				swap((*this)[i], (*this)[i + lstep]);
			}
		}
	}

	void P(size_t n, double fi = (PI/2), double error = 0)
	{
		size_t step = 1 << n;
		fi *= (1 + (error / 100));

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if(!((i % (1 << (n+1))) >> n))
			{
				(*this)[i + step] = (*this)[i + step] * (cos(fi)+im*sin(fi));
			}
		}
	}

	void CP(size_t h, size_t l, double fi = (PI / 2.0), double error = 0)
	{
		size_t lstep = 1 << l;
		fi *= (1 + (error / 100));

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (((i % (1 << (h + 1))) >> h) && !((i % (1 << (l + 1))) >> l))
			{
				(*this)[i + lstep] = (*this)[i + lstep] * (cos(fi) + im * sin(fi));
			}
		}
	}

	void CCP(size_t h1, size_t h2, size_t l, double fi = (PI / 2.0), double error = 0)
	{
		size_t lstep = 1 << l;
		fi *= (1 + (error / 100));

#pragma omp parallel for
		for (int i = 0; i < data.size(); i++)
		{
			if (((i % (1 << (h1 + 1))) >> h1) && ((i % (1 << (h2 + 1))) >> h2) && !((i % (1 << (l + 1))) >> l))
			{
				(*this)[i + lstep] = (*this)[i + lstep] * (cos(fi) + im * sin(fi));
			}
		}
	}

	//logic op (l - res, at the beginning l = 0 if strict calculations)
	// X(n) = NOT(n) in n
	//CCNOT - AND where l = h1 AND h2
	void NOT(size_t h, size_t l)//NOT(h) in l
	{
		X(l);
		CNOT(h, l);
	}

	void XOR(size_t h1, size_t h2, size_t l)
	{
		CNOT(h1, l);
		CNOT(h2, l);
	}

	void OR(size_t h1, size_t h2, size_t l)//de Morgan
	{
		X(h1);
		X(h2);
		CCNOT(h1, h2, l);
		X(h1);
		X(h2);
		X(l);
	}

	void SWAP(size_t start, size_t end)//(end-start)%2 == 0
	{
		for (size_t i = start; i < start + ((end - start) / 2); i++)
		{
			CNOT(i + ((end - start) / 2), i);
			CNOT(i, i + ((end - start) / 2));
			CNOT(i + ((end - start) / 2), i);
		}
	}

	void CSWAP(size_t start, size_t end, size_t u)//(end-start)%2 == 0
	{
		for (size_t i = start; i < start + ((end - start) / 2); i++)
		{
			CNOT(i + ((end - start) / 2), i);
			CCNOT(u, i, i + ((end - start) / 2));
			CNOT(i + ((end - start) / 2), i);
		}
	}

	//end = last+1
	void QFT(size_t start, size_t end, double error = 0)
	{
		for (size_t i = end - 1; i >= start && i + 1 != 0; i--)
		{
			this->H(i);
			for (size_t j = i - 1; j >= start && j + 1 != 0; j--)
			{
				this->CP(j, i, PI / (1 << (i - j)), error);
			}
		}
	}

	void RQFT(size_t start, size_t end, double error = 0)
	{
		for (size_t i = start; i < end; i++)
		{
			for (size_t j = start; j < i; j++)
			{
				this->CP(j, i, -PI / (1 << (i - j)), error);
			}
			this->H(i);
		}
	}

	//for all ADD and SUB, qb in Fourier
	void ADD(size_t a, size_t start, size_t end, double error = 0)//a < 2^n, a + b mod 2^n
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->P(i, PI / (1 << (i - j)), error);
				}
			}
		}
	}

	void CADD(size_t a, size_t start, size_t end, size_t u, double error = 0)
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->CP(u, i, PI / (1 << (i - j)), error);
				}
			}
		}
	}

	void CCADD(size_t a, size_t start, size_t end, size_t u1, size_t u2, double error = 0)
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->CCP(u1, u2, i, PI / (1 << (i - j)), error);
				}
			}
		}
	}

	void SUB(size_t a, size_t start, size_t end, double error = 0)//a < 2^n, a-b mod 2^n if b < a -> 2^(n-1) - (a-b)
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->P(i, -PI / (1 << (i - j)), error);
				}
			}
		}
	}

	void CSUB(size_t a, size_t start, size_t end, size_t u, double error = 0)
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->CP(u, i, -PI / (1 << (i - j)), error);
				}
			}
		}
	}

	void CCSUB(size_t a, size_t start, size_t end, size_t u1, size_t u2, double error = 0)
	{
		for (size_t j = start; j < end; j++)
		{
			for (size_t i = j; i < end; i++)
			{
				if ((a >> (j - start)) & 1)
				{
					this->CCP(u1, u2, i, -PI / (1 << (i - j)), error);
				}
			}
		}
	}

	//for all ModADD and ModSUB QFT(start, end-1)

	void ModADD(size_t a, size_t N, size_t start, size_t end, double error = 0)
	//a, b < N < 2^(end-start - 1), end-start >= max(bitsize(a, b))+2
	{
		ADD(a, start, end - 1, error);
		SUB(N, start, end - 1, error);

		RQFT(start, end - 1, error);
		CNOT(end - 2, end - 1);
		QFT(start, end - 1, error);

		CADD(N, start, end - 1, end - 1, error);
		SUB(a, start, end - 1, error);

		RQFT(start, end - 1, error);
		X(end - 2);
		CNOT(end - 2, end - 1);
		X(end - 2);
		QFT(start, end - 1, error);

		ADD(a, start, end - 1, error);
	}

	void CModADD(size_t a, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		CADD(a, start, end - 1, u, error);
		SUB(N, start, end - 1, error);

		RQFT(start, end - 1, error);
		CNOT(end - 2, end - 1);
		QFT(start, end - 1, error);

		CADD(N, start, end - 1, end - 1, error);
		CSUB(a, start, end - 1, u, error);

		RQFT(start, end - 1, error);
		X(end - 2);
		CNOT(end - 2, end - 1);
		X(end - 2);
		QFT(start, end - 1, error);

		CADD(a, start, end - 1, u, error);
	}

	void CCModADD(size_t a, size_t N, size_t start, size_t end, size_t u1, size_t u2, double error = 0)
	{
		CCADD(a, start, end - 1, u1, u2, error);
		SUB(N, start, end - 1, error);

		RQFT(start, end - 1, error);
		CNOT(end - 2, end - 1);
		QFT(start, end - 1, error);

		CADD(N, start, end - 1, end - 1, error);
		CCSUB(a, start, end - 1, u1, u2, error);

		RQFT(start, end - 1, error);
		X(end - 2);
		CNOT(end - 2, end - 1);
		X(end - 2);
		QFT(start, end - 1, error);

		CCADD(a, start, end - 1, u1, u2, error);
	}

	void ModSUB(size_t a, size_t N, size_t start, size_t end, double error = 0)
	//|a|, b < N < 2^(end-start - 1), end-start >= max(bitsize(a, b))+2
	{
		a = (N - a);
		ModADD(a, N, start, end, error);
	}

	void CModSUB(size_t a, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		a = (N - a);
		CModADD(a, N, start, end, u, error);
	}

	void CCModSUB(size_t a, size_t N, size_t start, size_t end, size_t u1, size_t u2, double error = 0)
	{
		a = (N - a);
		CCModADD(a, N, start, end, u1, u2, error);
	}

	void ModMULX(size_t a,size_t x, size_t N, size_t start, size_t end, double error = 0)
	//a, b < N < 2^(end-start - 1), end-start >= max(bitsize(a, b))+2
	{
		for (size_t j = start; j < end; j++)
		{
			if ((x >> (j - start)) & 1)
			{
				ModADD(a, N, start, end, error);
			}
			a = (a << 1) % N;
		}
	}

	void CModMULX(size_t a, size_t x, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		for (size_t j = start; j < end; j++)
		{
			if ((x >> (j - start)) & 1)
			{
				CModADD(a, N, start, end, u, error);
			}
			a = (a << 1) % N;
		}
	}

	//for all ModMUL and RModMUL and unitMUL QFT(start + ((end-start)/2), end-1)

	void ModMUL(size_t a, size_t N, size_t start, size_t end, double error = 0)
	//a, b < N < 2^(end-start - 1), end-start >= 2*max(bitsize(a, b))+3, end - start % 2 == 1
	{
		for (size_t j = start; j < start + ((end - start) / 2); j++)
		{
			CModADD(a, N, start + ((end - start) / 2), end, j, error);
			a = (a << 1) % N;
		}
	}

	void CModMUL(size_t a, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		for (size_t j = start; j < start + ((end - start) / 2); j++)
		{
			CCModADD(a, N, start + ((end - start) / 2), end, j, u, error);
			a = (a << 1) % N;
		}
	}

	void RModMUL(size_t a, size_t N, size_t start, size_t end, double error = 0)
	//a, b < N < 2^(end-start - 1), end-start >= 2*max(bitsize(a, b))+3, end - start % 2 == 1
	{
		for (size_t j = start; j < start + ((end - start) / 2); j++)
		{
			CModSUB(a, N, start + ((end - start) / 2), end, j, error);
			a = (a << 1) % N;
		}
	}

	void CRModMUL(size_t a, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		for (size_t j = start; j < start + ((end - start) / 2); j++)
		{
			CCModSUB(a, N, start + ((end - start) / 2), end, j, u, error);
			a = (a << 1) % N;
		}
	}

	void unitMUL(size_t a, size_t N, size_t start, size_t end, double error = 0)//(gcd(a, N) == 1, else output -> a*x mod N, x)
	//a, b < N < 2^(end-start - 1), end-start >= 2*max(bitsize(a, b))+3, end - start % 2 == 1
	{
		ModMUL(a, N, start, end, error);
		RQFT(start + ((end - start) / 2), end-1, error);
		SWAP(start, end-1);
		QFT(start + ((end - start) / 2), end-1, error);
		RModMUL(inverse_element_by_mod(a, N), N, start, end, error);
	}

	void CunitMUL(size_t a, size_t N, size_t start, size_t end, size_t u, double error = 0)
	{
		CModMUL(a, N, start, end, u, error);
		RQFT(start + ((end - start) / 2), end-1, error);
		CSWAP(start, end-1, u);
		QFT(start + ((end - start) / 2), end-1, error);
		CRModMUL(inverse_element_by_mod(a, N), N, start, end, u, error);
	}

	void Shor(size_t a, size_t N, size_t start, size_t end, double error = 0)//4n+3 zero qbit without QFT -> 0-127 and gcd(a, N) == 1
	{
		size_t n = (end - start - 3) / 4;
		(*this)[1 << (n*2)] = 1;
		size_t a_ = a;

		for (size_t i = 0; i < n*2; i++)
		{
			H(i);
		}

		QFT(n*3 + 1, end - 1, error);
		for (size_t i = 0; i < n*2; i++)
		{
			CunitMUL(a, N, n*2, end, (n*2 - i - 1), error);
			a = ModMul(a, a, N);
		}
		RQFT(n * 3 + 1, end - 1, error);
		
		RQFT(start, n*2, error);
	}

	void condition_exp(size_t start, size_t end, unsigned long int count = 10)
	{
		vector<size_t> temp(1i64 << (end - start));
		for(size_t i = 0; i < count; i++)
		{
			T r = T(rand()) / RAND_MAX;
			T sum = 0;
			size_t j = 0;
			while (sum < r)
			{
				sum += norm((*this)[j]);
				j++;
			}
			temp[((j-1) % (1i64 << end))>>start]++;
		}

		for (size_t i = 0; i < temp.size(); i++)
		{
			if (temp[i] != 0)
			{
				cout << '(' << i << "; " << temp[i] << ")\n";
			}
		}
	}

	void condition_exp_in_file(size_t start, size_t end, unsigned long int count, ofstream& out)
	{
		if (out.is_open())
		{
			vector<size_t> temp(1i64 << (end - start));
			for (size_t i = 0; i < count; i++)
			{
				T r = T(rand()) / RAND_MAX;
				T sum = 0;
				size_t j = 0;
				while (sum < r)
				{
					sum += norm((*this)[j]);
					j++;
				}
				temp[((j - 1) % (1i64 << end)) >> start]++;
			}

			out << size << ' ' << start << ' ' << end << ' ' << count << '\n';

			for (size_t i = 0; i < temp.size(); i++)
			{
				out << temp[i] << ' ';
			}
			out << '\n';
		}
		else
		{
			cout << "The files are not open!\n";
		}
	}

	friend istream& operator>>(istream& in, Qbit& v)
	{
		for (size_t i = 0; i < v.get_size(); i++)
		{
			in >> v[i];
		}
		return in;
	}
};
