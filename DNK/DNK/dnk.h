#include<iostream>
#include<bitset>
#include <string>
#include <stdexcept>
#include <exception>
#include<iomanip>

using namespace std;
using nucleotide = bitset<2>;

enum class Nucleotide { A, G, C, T }; //Набор нуклеотидов

nucleotide Set(Nucleotide value);

Nucleotide get(nucleotide value);

void pout(nucleotide value);
namespace laboratory {




	class RNK
	{
	private:

		size_t lenghtRNK;//Количество нуклеотидов в РНК 

		nucleotide *rnk; //Контейнер РНК / Хранит нуклеотиды

		size_t capacityRNK = 1;//Ёмкость контейнера 

	public:

		RNK(size_t lenght = 1000, Nucleotide value = Nucleotide::A); // Конструктор (задает начальную длину и значение) здесь создаем массив

		RNK(RNK const &a);

		RNK &Reallocate(); // переписывает в новый контейнер

		size_t Capacity() const; //Возвращает ёмкость

		size_t Length() const; //Возвращает длину 

		size_t Cardinality(Nucleotide value); //Количество нуклеотидов заданного значения

		RNK &operator=(RNK const &a);

		RNK &operator+(RNK &a); // соединить конец первой рнк с началом второй

		bool operator==(const RNK &a) const; // проверить на равенство каждый элемент цепи

		bool operator!=(const RNK &a) const;// обратно ==

		bool isComplimantary(RNK &a) const; //проверить на компланарность каждый элемент цепи в прямой и обратной, а также длину

		RNK &operator!();//перевести каждый элемент цепи в компланарный

		nucleotide operator[](size_t const indx) const;

		RNK &split(size_t indx); //Отрезает часть РНК до указанного номера

		void print() const;//Вывод РНК на экран

		explicit RNK(size_t capacityRNK);

		~RNK();// Деструктор

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