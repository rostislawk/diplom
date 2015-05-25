#include "utils.h"

uint64 RDTSC()
{
    __asm _emit 0x0F __asm _emit 0x31
}

void printBinaryRepresentation2(word *a, size_t arr_size, size_t m) 
{
	size_t i, length = m;
	for (i = 0; i < arr_size; ++i) {
		if (length > B_PER_W) {
			length = length - B_PER_W;
			printBinaryRepresentation1(a[i], B_PER_W);
			printf(" ");
		}
		else {
			printBinaryRepresentation1(a[i], length);
			break;
		}
	}
	printf("end\n");
}

void printBinaryRepresentation1(word a, size_t m)
{
	size_t index = 0;
	if (m > B_PER_W) {
		m = B_PER_W;
	}
	for (index = m; index > 0; --index) {
		printf("%u", (a>>(index-1))&1);
	}
}

void shiftRight(word *a, size_t n, size_t m, size_t shift)
{
	word carry;
	word rest;
	word abc;
	//normalize(a, n, m);
	if (shift > B_PER_W) {
		carry = wordShHiCarry(a, n, B_PER_W, 0);
		wordSetBits(a, m - B_PER_W ,B_PER_W, carry);
		shiftRight(a, n, m, shift - B_PER_W);
	}
	else {
		abc = (m >= B_PER_W) ? m - B_PER_W : 0; 
		carry = wordShHiCarry(a, n, shift, 0);
		if (m < B_PER_W) {
			wordShHi(&carry, 1, B_PER_W - m);
		}
		rest = wordGetBits(a, abc, B_PER_W);
		carry = carry | rest;
		wordSetBits(a, abc, B_PER_W, carry);
	}
}

void normalize(word *a, size_t n, size_t m) 
{
	wordTrimHi(a, n, m);
}

size_t bits_in_number(size_t number) 
{
	size_t index = 0;
	while (number != 0) {
		index++;
		number = number >> 1;
	}
	return index;
}

size_t next_power_of_two(size_t number)
{
	size_t index = 1;
	while (number != 0) {
		index << 1;
		number = number >> 1;
	}
	return number;	
}

size_t size_in_words(size_t number)
{
	word rest = number % B_PER_W;
	size_t size_in_words = number / B_PER_W;
	return rest == 0 ? size_in_words : size_in_words + 1;
};

octet reverseOctet(octet a)
{
	a = (a & 0xF0) >> 4 | (a & 0x0F) << 4;
	a = (a & 0xCC) >> 2 | (a & 0x33) << 2;
	a = (a & 0xAA) >> 1 | (a & 0x55) << 1;
	return a;
}

word reverseWord(word a)
{
	a = (a & 0xFFFF0000) >> 16 | (a & 0x0000FFFF) << 16;
	a = (a & 0xFF00FF00) >> 8 | (a & 0x00FF00FF) << 8;
	a = (a & 0xF0F0F0F0) >> 4 | (a & 0x0F0F0F0F) << 4;
	a = (a & 0xCCCCCCCC) >> 2 | (a & 0x33333333) << 2;
	a = (a & 0xAAAAAAAA) >> 1 | (a & 0x55555555) << 1;
	return a;
}

void reverse(word *a, word *b, size_t n)
{
	size_t index = 0;
	for (index = 0; index < n / 2; ++index) {
		b[index] = reverseWord(a[n-index-1]);
		b[n-index-1] = reverseWord(a[index]);
	}
}
