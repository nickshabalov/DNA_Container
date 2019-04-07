#include<iostream>
#include <stdexcept>
#include <exception>
#include<bitset>
#include <string>
#include<iomanip>
#include"dnk.h"

using namespace std;
using namespace laboratory;


nucleotide  Set(Nucleotide value) {
	switch (value) {
	case Nucleotide::A: return nucleotide("00"); break;
	case Nucleotide::G: return nucleotide("01"); break;
	case Nucleotide::C: return nucleotide("10"); break;
	case Nucleotide::T: return nucleotide("11"); break;
	}
}

Nucleotide get(nucleotide value) {
	switch (value.to_ulong()) {
	case 0: return Nucleotide::A; break;
	case 1: return Nucleotide::G; break;
	case 2: return Nucleotide::C; break;
	case 3: return Nucleotide::T; break;
	default: return Nucleotide::A; break;
	}
}

RNK::RNK(size_t lenght, Nucleotide value)
{

	lenghtRNK = lenght;
	while (lenghtRNK > capacityRNK)
	{
		capacityRNK = capacityRNK * 2;
	}
	this->rnk = new nucleotide[capacityRNK];
	for (size_t i = 0; i < lenghtRNK; i++)
	{
		this->rnk[i] = Set(value);
	};
}

RNK::RNK(const RNK &a)
{
	lenghtRNK = a.Length();
	while (lenghtRNK > capacityRNK)
	{
		capacityRNK = capacityRNK * 2;
	}
	this->rnk = new nucleotide[capacityRNK];
	for (size_t i = 0; i < a.Length(); i++)
	{
		this->rnk[i] = a.rnk[i];
	}
}

RNK::~RNK()
{
	delete[] rnk;
}

RNK &RNK::Reallocate()
{
	nucleotide* temp = new nucleotide[capacityRNK];
	for (size_t i = 0; i < lenghtRNK; i++)
	{
		temp[i] = rnk[i];
	}
	delete[] rnk;
	rnk = temp;
	return *this;
}

size_t RNK::Capacity()const
{
	return capacityRNK;
}

size_t RNK::Length() const
{
	return lenghtRNK;
}

size_t RNK::Cardinality(Nucleotide value)
{
	size_t count = 0;
	for (size_t i = 0; i < lenghtRNK; i++)
	{
		if (get(rnk[i]) == value)
		{
			count++;
		}
	}
	return count;
}

RNK &RNK::operator+(RNK &a)
{
	lenghtRNK = lenghtRNK + a.Length();
	if (capacityRNK < lenghtRNK)
	{
		while (lenghtRNK > capacityRNK)
		{
			capacityRNK = capacityRNK * 2;
		}
	}
	nucleotide* temp = new nucleotide[capacityRNK];
	for (size_t i = 0; i < (lenghtRNK - a.Length()); i++)
	{
		temp[i] = rnk[i];
	}
	for (size_t i = lenghtRNK - a.Length(); i < lenghtRNK; i++)
	{
		temp[i] = a.rnk[i - (lenghtRNK - a.Length())];
	}
	delete[] rnk;
	rnk = temp;
	return *this;
}

bool RNK::operator==(const RNK &a) const
{
	if (lenghtRNK != a.lenghtRNK)
	{
		return false;
	}
	else
	{
		for (size_t i = 0; i < lenghtRNK; i++)
		{
			if ((this->rnk[i] != a.rnk[i]) && (this->rnk[i] != a.rnk[lenghtRNK - i - 1]))
			{
				return false;
			}
		}
	}
	return true;
}

bool RNK::operator!=(const RNK &a) const
{
	return !(operator==(a));
}

RNK &RNK::operator!()
{
	for (size_t i = 0; i < lenghtRNK; i++)
	{
		switch (get(rnk[i]))
		{
		case Nucleotide::A:
		{
			this->rnk[i] = nucleotide("11");
			break;
		}
		case Nucleotide::G:
		{
			this->rnk[i] = nucleotide("10");
			break;
		}
		case Nucleotide::C:
		{
			this->rnk[i] = nucleotide("01");
			break;
		}
		case Nucleotide::T:
		{
			this->rnk[i] = nucleotide("00");
			break;
		}

		}
	}
	return *this;
}

RNK &RNK::operator=(RNK const &a)
{
	if (this != &a)
	{
		nucleotide* temp = new nucleotide[capacityRNK];
		delete[] rnk;
		this->capacityRNK = a.Capacity();
		this->lenghtRNK = a.Length();
		rnk = temp;
		for (size_t i = 0; i < lenghtRNK; i++)
		{
			rnk[i] = a.rnk[i];
		}
	}
	return *this;
}

bool RNK::isComplimantary(RNK &a) const
{
	!a;
	if (a == *this)
	{
		!a;
		return true;
	}

}

RNK &RNK::split(size_t indx)
{
	if (indx > lenghtRNK)
	{
		throw length_error("Can't split RNA that haven't enough Nucleotides");
	}
	else
	{
		lenghtRNK = indx;
		capacityRNK = 1;
		if (lenghtRNK >= capacityRNK)
		{
			while (lenghtRNK >= capacityRNK)
			{
				capacityRNK = capacityRNK * 2;
			}
			this->Reallocate();
		}
		nucleotide* temp = new nucleotide[lenghtRNK];
		for (size_t i = 0; i < lenghtRNK; i++)
		{
			temp[i] = rnk[i];
		}
		delete[] rnk;
		rnk = temp;
		return *this;
	}
}

void RNK::print() const
{
	for (size_t i = 0; i < lenghtRNK; i++)
	{
		cout << "|" << rnk[i];
	}
}

RNK::Reference RNK::operator[](const size_t indx)
{
	if (indx > capacityRNK)
	{
		while (indx > capacityRNK)
		{
			capacityRNK = capacityRNK * 2;
		}
		this->Reallocate();
	}
	if (indx + 1 > lenghtRNK)
	{
		for (size_t j = lenghtRNK; j < indx; j++)
		{
			rnk[j] = nucleotide("00");
		}
		lenghtRNK = indx;
	}
	return Reference(*this, indx);
}

nucleotide RNK::operator[](size_t const indx) const
{
		return rnk[indx];
}
void RNK::Reference::operator=(nucleotide assgn)
{
	parent.rnk[indx - 1] = assgn;
}

void RNK::Reference::operator=(const RNK::Reference &assgn)
{
	parent.rnk[indx-1] = assgn.Get();
}

void RNK::Reference::operator=(const Nucleotide assgn)
{
	parent.rnk[indx - 1] = Set(assgn);
}

void RNK::push(const Nucleotide value)
{
	if (lenghtRNK >= capacityRNK)
	{
		{
			capacityRNK = capacityRNK * 2;
		}
		this->Reallocate();
	}
	this->rnk[lenghtRNK] = Set(value);
	lenghtRNK++;
}

DNK::DNK(RNK &left, RNK &right)
{
	if (left.isComplimantary(right) == false)
	{
		throw length_error("Impossible to make dnk; RNK aren`t complantary");

	}
	else
	{
		left_rnk = new RNK(left);
	}
}

DNK::DNK(const DNK &right)
{
	delete left_rnk;
	this->left_rnk = new RNK(*right.left_rnk);
}

DNK &DNK::operator=(const DNK &right)
{
	if (&right != this) {
		delete left_rnk;
		left_rnk = new RNK(*right.left_rnk);
	}
	return *this;
}

RNK &DNK::getRight()
{

	return (*left_rnk);
}

RNK &DNK::getLeft() {
	return *left_rnk;
}

RNK::Reference::operator nucleotide const&()
{
	return parent.rnk[indx];
}

nucleotide laboratory::RNK::Reference::Get() const
{
	return parent.rnk[indx];
}

RNK::Reference::operator Nucleotide const&()
{
	return get(parent.rnk[indx]);
}

bool RNK::Reference::operator==(const Reference &assgn) const
{
	return (parent.rnk == assgn.parent.rnk);
}

bool laboratory::RNK::Reference::operator==(const nucleotide value) const
{
	return (parent.rnk[indx] == value);
}

bool RNK::Reference::operator!=(const Reference &assgn) const
{
	return (parent.rnk != assgn.parent.rnk);
}

void pout(nucleotide value)
{
	switch (value.to_ulong()) {
	case 0: cout << "A"; break;
	case 1: cout << "G"; break;
	case 2: cout << "C"; break;
	case 3: cout << "T"; break;
	default: break;
	}
}