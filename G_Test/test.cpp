#include "pch.h"

#include "..//ITLab/Qbit.h"
#include "..//ITLab/RM_RSA.h"

TEST(Qbit, normalization)
{
	Qbit<double> test(2);
	test[1] = 2;
	test[2] = 2;
	test.normalization();
	EXPECT_EQ(test[1], test[2]);
	EXPECT_FLOAT_EQ(norm(test[1]) + norm(test[2]), 1);
}

TEST(Qbit, attach)
{
	Qbit<double> test1(2);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[3] = 1;
	test2[0] = -1;
	test2[1] = 1;
	test1.normalization();
	test2.normalization();
	test1.attach(test2);
	Qbit<double> test3(4);
	test3[0] = -1;
	test3[3] = -1;
	test3[4] = 1;
	test3[7] = 1;
	test3.normalization();
	for (size_t i = 0; i < 8; i++)
	{
		EXPECT_TRUE(norm(test1[i] - test3[i]) < 0.00000001);
	}
}

TEST(Qbit, gate_H)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[0] = 1;
	test2[1] = 1;
	test1.H(0);
	test2.normalization();
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit, gate_X)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[1] = 1;
	test1.X(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])),  0.0);
	}
}

TEST(Qbit, gate_Y)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test2[1] = im;
	test1.Y(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit, gate_Z)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[1] = 1;
	test1.normalization();
	test2[0] = 1;
	test2[1] = -1;
	test2.normalization();
	test1.Z(0);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit, gate_P)
{
	Qbit<double> test1(1);
	Qbit<double> test2(1);
	test1[0] = 1;
	test1[1] = 1;
	test1.normalization();
	test2[0] = 1;
	test2[1] = im;
	test2.normalization();
	test1.P(PI/4);
	for (size_t i = 0; i < 2; i++)
	{
		EXPECT_TRUE((norm(test1[i] - test2[i])) < 0.00000001);
	}
}

TEST(Qbit, gate_CU)
{
	Qbit<double> test1(2);
	Qbit<double> test2(2);
	test1[2] = 1;
	test1[3] = 1;
	test1.normalization();
	test2[1] = 1;
	test2[2] = 1;
	test2.normalization();
	test1.CNOT(0, 1);
	for (size_t i = 0; i < 4; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])),  0.0);
	}
}

TEST(Qbit, SWAP)
{
	Qbit<double> test1(3);
	Qbit<double> test2(3);
	test1[3] = 1;
	test1[4] = 1;
	test1.normalization();
	test2[1] = 1;
	test2[6] = 1;
	test2.normalization();
	test1.SWAP(0, 2);
	for (size_t i = 0; i < 8; i++)
	{
		EXPECT_DOUBLE_EQ((norm(test1[i] - test2[i])), 0.0);
	}
}

TEST(Qbit, QFT_RQFT)
{
	Qbit<double> test1(3);
	test1[0] = 1;
	test1.QFT(0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[0]) < 0.00000001);
}

TEST(Qbit, ADD_SUB)
{
	Qbit<double> test1(3);
	test1[3] = 1;
	test1.QFT(0, 3);
	test1.ADD(2, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[5]) < 0.00000001);
	test1.QFT(0, 3);
	test1.SUB(3, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);
	test1.QFT(0, 3);
	test1.ADD(7, 0, 3);
	test1.RQFT(0, 3);
	EXPECT_TRUE(1.0 - norm(test1[1]) < 0.00000001);
	test1.QFT(1, 3);
	test1.SUB(2, 1, 3);
	test1.RQFT(1, 3);
	EXPECT_TRUE(1.0 - norm(test1[5]) < 0.00000001);
}

TEST(Qbit, ADD_SUB_BY_MODULE)
{
	Qbit<double> test1(5);
	test1[1] = 1;
	test1.QFT(0, 4);
	test1.ModADD(1, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);

	test1.QFT(0, 4);
	test1.ModADD(3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[1]) < 0.00000001);

	test1.QFT(0, 4);
	test1.ModSUB(3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[2]) < 0.00000001);
}

TEST(Qbit, MUL_RMUL_MULX_BY_MODULE)
{
	Qbit<double> test1(5);
	test1[1] = 1;
	test1.QFT(0, 4);
	test1.ModMULX(2, 3, 4, 0, 5);
	test1.RQFT(0, 4);
	EXPECT_TRUE(1.0 - norm(test1[3]) < 0.00000001);

	Qbit<double> test2(6);
	test2[2] = 1;
	test2.QFT(2, 5);
	test2.ModMUL(2, 3, 0, 6);
	test2.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test2[2 + (1<<2)]) < 0.00000001);

	Qbit<double> test3(6);
	test3[3] = 1;
	test3.QFT(2, 5);
	test3.RModMUL(3, 4, 0, 6);
	test3.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test3[3 + (3 << 2)]) < 0.00000001);
}

TEST(Qbit, unitMUL)
{
	Qbit<double> test2(6);
	test2[3] = 1;
	test2.QFT(2, 5);
	test2.unitMUL(3, 4, 0, 6);
	test2.RQFT(2, 5);
	EXPECT_TRUE(1.0 - norm(test2[1]) < 0.00000001);
}

TEST(Qbit, SHOR)
{
	Qbit<double> test2(18);
	test2.Shor(2, 15, 0, 18, 0);
	vector<size_t> temp = test2.condition_exp(0, 8, (1i64 << 11));
	for (size_t i = 0; i < 256; i++)
	{
		if (temp[i] != 0)
		{
			EXPECT_TRUE((i % 64) == 0);
		}
	}
}

TEST(RMR_SA, MillerRabin)
{
	vector<unsigned long long int> RM;
	vector<unsigned long long int> res;
	res.push_back(2);
	res.push_back(2);
	res.push_back(3);
	res.push_back(7);
	res.push_back(5);
	res.push_back(7);
	MillerRabin(2940, RM);
	for (size_t i = 0; i < RM.size(); i++)
	{
		EXPECT_EQ(RM[i], res[i]);
	}
	RM.erase(RM.begin(), RM.end());
	res.erase(res.begin(), res.end());
	res.push_back(163);
	res.push_back(878167);
	MillerRabin(143141221, RM);
	for (size_t i = 0; i < RM.size(); i++)
	{
		EXPECT_EQ(RM[i], res[i]);
	}
}

TEST(RMR_SA, RSA)
{
	unsigned int m = 4;
	EXPECT_EQ(RSA_decryption(RSA_encryption(m)), m);
}