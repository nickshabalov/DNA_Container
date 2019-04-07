#include<iostream>
#include<bitset>
#include <string>
#include <stdexcept>
#include <exception>
#include<iomanip>

using namespace std;
using nucleotide = bitset<2>;

enum class Nucleotide { A, G, C, T }; //����� �����������

nucleotide Set(Nucleotide value);

Nucleotide get(nucleotide value);

void pout(nucleotide value);
namespace laboratory {




	class RNK
	{
	private:

		size_t lenghtRNK;//���������� ����������� � ��� 

		nucleotide *rnk; //��������� ��� / ������ ����������

		size_t capacityRNK = 1;//������� ���������� 

	public:

		RNK(size_t lenght = 1000, Nucleotide value = Nucleotide::A); // ����������� (������ ��������� ����� � ��������) ����� ������� ������

		RNK(RNK const &a);

		RNK &Reallocate(); // ������������ � ����� ���������

		size_t Capacity() const; //���������� �������

		size_t Length() const; //���������� ����� 

		size_t Cardinality(Nucleotide value); //���������� ����������� ��������� ��������

		RNK &operator=(RNK const &a);

		RNK &operator+(RNK &a); // ��������� ����� ������ ��� � ������� ������

		bool operator==(const RNK &a) const; // ��������� �� ��������� ������ ������� ����

		bool operator!=(const RNK &a) const;// ������� ==

		bool isComplimantary(RNK &a) const; //��������� �� �������������� ������ ������� ���� � ������ � ��������, � ����� �����

		RNK &operator!();//��������� ������ ������� ���� � ������������

		nucleotide operator[](size_t const indx) const;

		RNK &split(size_t indx); //�������� ����� ��� �� ���������� ������

		void print() const;//����� ��� �� �����

		explicit RNK(size_t capacityRNK);

		~RNK();// ����������

		class Reference
		{
		private:
			RNK &parent;

			size_t indx;
		public:
			Reference(RNK &p, size_t x) : parent(p), indx(x)
			{
			}

			void operator = (const nucleotide assgn);

			void operator = (const Reference &assgn);

			void operator = (const Nucleotide assgn);

			operator nucleotide const &();

			operator Nucleotide const &();

			bool operator == (const Reference & assgn) const;

			bool operator == (const nucleotide value) const;

			bool operator != (const Reference & assgn) const;

			nucleotide Get() const;
		};
		Reference operator[](const size_t indx);

		void push(Nucleotide value);
	};

	class DNK
	{
	private:
		RNK * left_rnk;
	public:

		DNK(RNK &left, RNK &right);
		DNK(const DNK &right);
		DNK &operator=(const DNK &right);
		RNK &getRight();
		RNK &getLeft();
	};

}