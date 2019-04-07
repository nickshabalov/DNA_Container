#include "pch.h"
//#include "../DNK/dnk.h"
#include "../DNK/dnk.cpp"

using namespace laboratory;

TEST(RNKTest, RNKisComplimentary)
{
	RNK rnk1(10, Nucleotide::A);
	RNK rnk2(10, Nucleotide::T);
	EXPECT_TRUE(rnk1.isComplimantary(rnk2));
	rnk2[15] = nucleotide("11");
	EXPECT_FALSE(rnk1.isComplimantary(rnk2));
	rnk1[15] = nucleotide("00");
	EXPECT_FALSE(rnk1.isComplimantary(rnk2));
	for (size_t i = 10; i <= 15; i++)
	{
		rnk2[i] = nucleotide("11");
	}
	EXPECT_TRUE(rnk1.isComplimantary(rnk2));
}

TEST(RNKTest, RNKComplimentary)
{
	RNK rnk1(10, Nucleotide::A);
	RNK rnk2(10, Nucleotide::T);
	!rnk1;
	EXPECT_TRUE(rnk1 == rnk2);
	!rnk1;
	!rnk2;
	EXPECT_TRUE(rnk2 == rnk1);
	EXPECT_EQ(!(!rnk1), rnk1);
}

TEST(RNKTest, RNKisEqual)
{
	RNK rnk1(10, Nucleotide::A);
	RNK rnk2(10, Nucleotide::A);
	RNK rnk3(10, Nucleotide::T);
	EXPECT_TRUE(rnk1 == rnk2);
	EXPECT_TRUE(rnk1 != rnk3);
	EXPECT_TRUE(rnk1 == (!rnk3));
	EXPECT_TRUE(rnk1 == rnk3);
	rnk2[15] = nucleotide("11");
	EXPECT_TRUE(rnk1 != rnk2);
	for (size_t i = 10; i <= 15; i++)
	{
		rnk1[i] = nucleotide("11");
	}
	EXPECT_TRUE(rnk1 != rnk2);
	for (size_t i = 10; i <= 15; i++)
	{
		rnk2[i] = nucleotide("11");
	}
	EXPECT_TRUE(rnk1 == rnk2);
}

TEST(RNKTest, RNKLenght)
{
	RNK rnk(10, Nucleotide::A);
	EXPECT_EQ(rnk.Length(), 10);
	EXPECT_EQ(rnk.Capacity(), 16);
}

TEST(RNKTest, RNKCardinality)
{
	RNK rnk(10, Nucleotide::A);
	EXPECT_EQ(rnk.Cardinality(Nucleotide::A), 10);
	EXPECT_EQ(rnk.Cardinality(Nucleotide::T), 0);
	EXPECT_EQ(rnk.Cardinality(Nucleotide::C), 0);
	for (size_t i = 1; i < rnk.Length(); i = i + 4)
	{
		rnk[i] = Nucleotide::G;
	}
	EXPECT_EQ(rnk.Cardinality(Nucleotide::G), 3);
}

TEST(RNKTest, RNKSum)
{

	RNK rnk1(10, Nucleotide::A);
	RNK rnk2(5, Nucleotide::T);
	RNK rnk3(10, Nucleotide::A);
	for (size_t i = 10; i <= 15; i++)
	{
		rnk3[i] = Nucleotide::T;
	}
	rnk1 + rnk2;
	EXPECT_EQ(rnk1.Length(), 15);
	EXPECT_EQ(rnk1.Capacity(), 16);
	
	EXPECT_TRUE(rnk1 == rnk3);
}
TEST(RNKTest, RNKSplit)
{
	RNK rnk1(10, Nucleotide::A);
	RNK rnk2(5, Nucleotide::A);
	rnk1.split(5);
	EXPECT_EQ(rnk1.Length(), 5);
	EXPECT_TRUE(rnk1 == rnk2);
}
TEST(RNKTest, RNKControlTests)
{
	RNK rnk(10, Nucleotide::A);
	rnk[2] = Nucleotide::T;
	rnk[3] = rnk[1];
	nucleotide n = rnk[1];

	EXPECT_TRUE(rnk[2] == nucleotide("11"));
	EXPECT_TRUE(rnk[3] == nucleotide("00"));
	EXPECT_EQ(n, nucleotide("11"));

	const RNK rnk2(10, Nucleotide::A);
	nucleotide n2 = rnk2[1];
	EXPECT_EQ(n2, nucleotide("00"));

	
}

TEST(RNKTest, RNKLeak)
{
	RNK rnk(10, Nucleotide::A);
	for (size_t i = 0; i < 1000000; i++)
	{
		rnk.push(Nucleotide::T);
	}
	rnk.split(500000);
}