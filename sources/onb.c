#include "onb.h"
#include "utils.h"

static octet f[256] = {
	0x00, 0x01, 0x02, 0x03, 0x05, 0x04, 0x07, 0x06, 0x08, 0x09, 0x0A, 0x0B, 0x0D, 0x0C, 0x0F, 0x0E,
	0x15, 0x14, 0x17, 0x16, 0x10, 0x11, 0x12, 0x13, 0x1D, 0x1C, 0x1F, 0x1E, 0x18, 0x19, 0x1A, 0x1B,
	0x22, 0x23, 0x20, 0x21, 0x27, 0x26, 0x25, 0x24, 0x2A, 0x2B, 0x28, 0x29, 0x2F, 0x2E, 0x2D, 0x2C,
	0x37, 0x36, 0x35, 0x34, 0x32, 0x33, 0x30, 0x31, 0x3F, 0x3E, 0x3D, 0x3C, 0x3A, 0x3B, 0x38, 0x39,
	0x51, 0x50, 0x53, 0x52, 0x54, 0x55, 0x56, 0x57, 0x59, 0x58, 0x5B, 0x5A, 0x5C, 0x5D, 0x5E, 0x5F,
	0x44, 0x45, 0x46, 0x47, 0x41, 0x40, 0x43, 0x42, 0x4C, 0x4D, 0x4E, 0x4F, 0x49, 0x48, 0x4B, 0x4A,
	0x73, 0x72, 0x71, 0x70, 0x76, 0x77, 0x74, 0x75, 0x7B, 0x7A, 0x79, 0x78, 0x7E, 0x7F, 0x7C, 0x7D,
	0x66, 0x67, 0x64, 0x65, 0x63, 0x62, 0x61, 0x60, 0x6E, 0x6F, 0x6C, 0x6D, 0x6B, 0x6A, 0x69, 0x68,
	0x80, 0x81, 0x82, 0x83, 0x85, 0x84, 0x87, 0x86, 0x88, 0x89, 0x8A, 0x8B, 0x8D, 0x8C, 0x8F, 0x8E,
	0x95, 0x94, 0x97, 0x96, 0x90, 0x91, 0x92, 0x93, 0x9D, 0x9C, 0x9F, 0x9E, 0x98, 0x99, 0x9A, 0x9B,
	0xA2, 0xA3, 0xA0, 0xA1, 0xA7, 0xA6, 0xA5, 0xA4, 0xAA, 0xAB, 0xA8, 0xA9, 0xAF, 0xAE, 0xAD, 0xAC,
	0xB7, 0xB6, 0xB5, 0xB4, 0xB2, 0xB3, 0xB0, 0xB1, 0xBF, 0xBE, 0xBD, 0xBC, 0xBA, 0xBB, 0xB8, 0xB9,
	0xD1, 0xD0, 0xD3, 0xD2, 0xD4, 0xD5, 0xD6, 0xD7, 0xD9, 0xD8, 0xDB, 0xDA, 0xDC, 0xDD, 0xDE, 0xDF,
	0xC4, 0xC5, 0xC6, 0xC7, 0xC1, 0xC0, 0xC3, 0xC2, 0xCC, 0xCD, 0xCE, 0xCF, 0xC9, 0xC8, 0xCB, 0xCA,
	0xF3, 0xF2, 0xF1, 0xF0, 0xF6, 0xF7, 0xF4, 0xF5, 0xFB, 0xFA, 0xF9, 0xF8, 0xFE, 0xFF, 0xFC, 0xFD,
	0xE6, 0xE7, 0xE4, 0xE5, 0xE3, 0xE2, 0xE1, 0xE0, 0xEE, 0xEF, 0xEC, 0xED, 0xEB, 0xEA, 0xE9, 0xE8};

// Возведение в квадрат в оптимальном нормальном базисе
void sqr(word *res, word *a, size_t n, size_t m)
{
	wordCopy(res, a, n);
	shiftRight(res, n, m, 1);
}

void generate_b(word *b, size_t m)
{
	size_t index = 0, k, bit_size = bits_in_number(m);
	word M, N;
	BOOL res;
	for (index = 0; index < m + 1; ++index) {
		k = (index + m) / 2;
		M = index;
		N = k;
		M = ~M;
		M = ~(M | N);
		normalize(&M, 1, bit_size);
		res = (M == 0) ? TRUE : FALSE;
		wordSetBit(b, index, res);
	}
	normalize(b, bit_size, m + 1);
}

void generate_pi(word *pi, size_t m)
{
	size_t index = 0;
	size_t p = 2 * m + 1;
	pi[0] = 1;
	for (index = 1; index < m; ++index) {
		pi[index] = (pi[index-1] * 2) % p;
		if (pi[index] > m) {
			pi[index] = p - pi[index];
		}
		pi[index - 1] -= 1;
	}
	pi[m - 1] -= 1;
}

void applyPi(word *a, word *b, word *pi, size_t m)
{
	size_t index = 0;
	BOOL bit = FALSE;
	for (index = 0; index < m; ++index) {
		bit = wordGetBits(a, index, 1) == 1 ? TRUE : FALSE;
		wordSetBit(b, pi[index], bit); 
	}
}

void applyFToWord(word *a)
{
	octet *b, *c, carry = 0;
	b = (octet *)a;
	c = (octet *)malloc(2 * sizeof(octet));

	//apply f to 4 octets
	c[0] = reverseOctet(b[3]);
	c[1] = reverseOctet(b[2]);
	c[0] <<= 1;
	carry = c[1] & 0x80;
	c[1] <<= 1;
	c[0] += carry;
	b[0] ^= reverseOctet(c[0]);
	b[1] ^= reverseOctet(c[1]);
	
	//apply f to 2,2 octets
	c[0] = reverseOctet(b[1]);
	c[0] <<= 1;
	b[0] ^= reverseOctet(c[0]);

	c[1] = reverseOctet(b[3]);
	c[1] <<= 1;
	b[2] ^= reverseOctet(c[1]);

	//apply f to 1,1,1,1 octets
	b[0] = f[b[0]];
	b[1] = f[b[1]];
	b[2] = f[b[2]];
	b[3] = f[b[3]];

	free(c);
}

void applyF(word *a, word *mem, size_t n)
{
	size_t index = 0, m = n / 2;
	if (m > 1) {
		wordCopy(mem, a + m, m);
		reverse(mem, mem, m);
		wordShHi(mem, m, 1);
		for (index = 0; index < m; ++index) {
			a[index] ^= mem[index];
		}
		applyF(a, mem, m);
		applyF(a + m, mem, m);
	}
	else {
		applyFToWord(a);	
	}
}

void generateONB2_A(word *a, size_t m)
{
	size_t p = 2 * m + 1;
	size_t ksigma, kmu;
	word *pi = (word *)malloc(p * sizeof(word));
	word *pi_inv = (word *)malloc(p * sizeof(word));
	word *sigma = (word *)malloc((p - 1) * sizeof(word));
	word *mu = (word *)malloc((p - 1) * sizeof(word));
	size_t index;
	pi[0] = 1;
	pi[p-1] = 1;
	pi_inv[0] = WORD_MAX;
	pi_inv[1] = 0;
	for (index = 1; index < p - 1; ++index) {
		pi[index] = (2 * pi[index - 1]) % p;
		pi_inv[pi[index]] = index;
	}
	sigma[0] = 1;
	mu[0] = WORD_MAX;
	for (index = 1; index < p; ++index) {
		if (index > 0) {
			ksigma = pi_inv[index - 1];
			if (ksigma != WORD_MAX && ksigma < m) {
				sigma[ksigma] = pi_inv[1+pi[ksigma]] % m;
			}
			else if (ksigma < m) {
				sigma[ksigma] = WORD_MAX;
			}
		}
		if (index < p - 1) {
			kmu = pi_inv[index + 1];
			if (kmu != WORD_MAX && kmu < m) {
				mu[kmu] = pi_inv[-1+pi[kmu]] % m;
			}
			else if (kmu < m) {
				mu[kmu] = WORD_MAX;
			}
		}
	}
	a[0] = sigma[0];
	for (index = 1; index < m; ++index) {
		a[index * 2 - 1] = sigma[index];
		a[index * 2] = mu[index];
	}
	free(sigma);
	free(mu);
	free(pi);
	free(pi_inv);
}

void generateONB3_A(word *a, size_t m)
{
	size_t p = 2 * m + 1;
	size_t ksigma, kmu;
	word *pi = (word *)malloc(p * sizeof(word));
	word *pi_inv = (word *)malloc(p * sizeof(word));
	word *sigma = (word *)malloc((p - 1) * sizeof(word));
	word *mu = (word *)malloc((p - 1) * sizeof(word));
	size_t index;
	pi[0] = 1;
	pi[p-1] = 1;
	for (index = 0; index < p; ++index) {
		pi_inv[index] = WORD_MAX;
	}
	pi_inv[1] = 0;
	for (index = 1; index < m; ++index) {
		pi[index] = (2 * pi[index - 1]) % p;
		pi[index + m];
		pi_inv[pi[index]] = index;
	}
	sigma[0] = 1;
	mu[0] = WORD_MAX;
	
	for (index = 1; index < m; ++index) {
		ksigma = pi[index] + 1;
		if (pi_inv[ksigma] == WORD_MAX) {
			ksigma = p - pi[index] - 1;
		}
		sigma[index] = pi_inv[ksigma];
		kmu = pi[index] - 1;
		if (pi_inv[kmu] == WORD_MAX) {
			kmu = p - pi[index] + 1;
		}
		mu[index] = pi_inv[kmu];
	}
	a[0] = sigma[0];
	for (index = 1; index < m; ++index) {
		a[index * 2 - 1] = sigma[index];
		a[index * 2] = mu[index];
	}
	free(sigma);
	free(mu);
	free(pi);
	free(pi_inv);
}

void mul(word *res, word *a, word *b, word *A, size_t n, size_t m)
{
	size_t step = 0;
	size_t p = 2 * m - 1;
	size_t index = 0;
	size_t apos;
	size_t bpos;
	size_t aposm, bposm;
	BOOL xi, yj, sum;
	BOOL bit = FALSE;
	wordSetZero(res, n);
	for (step = 0; step < m; ++step) {
		bit = FALSE;
		for (index = 0; index < p; ++index) {
			aposm = ((index + 1) >> 1) + step;
			bposm = A[index] + step;
			apos = aposm > m ? aposm - m: aposm;
			bpos = bposm > m ? bposm - m : bposm;
			xi = wordTestBit(a, apos);
			yj = wordTestBit(b, bpos); 
			bit = bit ^ (xi && yj);
		}
		wordSetBit(res, step, bit);
	}
}

void inv(word *res, word *a, word *A, size_t n, size_t m) 
{
	word *stack1 = (word *)malloc(n * sizeof(word));
	word degree = m - 1;
	size_t bitSize = wordBitSize(&degree, 1);
	size_t index;
	sqr(stack1, a, n, m);
	for (index = bitSize - 1; index > 0; --index) { 
		if (wordGetBits(&degree, index - 1, 1) & 1) {
			mul(res, a, stack1, A, n, m);
		}
		else {
			wordCopy(res, stack1, n);
		}
		sqr(stack1, res, n, m);
	}
	free(stack1);
}

void div_onb(word *res, word *a, word *b, word *A, size_t n, size_t m)
{
	word *b_inv = (word *)malloc(n * sizeof(word));
	inv(b_inv, b, A, n, m);
	mul(res, a, b_inv, A, n, m);
	free(b_inv);
}

void fromONB2ToStandard(word *onb2, word *st, word *b, word *pi, size_t m)
{
	size_t ext = (next_power_of_two(m) > B_PER_W) ? next_power_of_two(m) : B_PER_W;
	word *mem = (word *)malloc(ext);
	wordSetZero(st, ext / B_PER_W);
	applyPi(onb2, st, pi, m);
	applyF(st, mem, ext / B_PER_W);
	if (wordGetBits(st, m - 1, 1) == 1) { 
		wordShHi(st, ext / B_PER_W, 1);
		wordXor(st, st, b, ext / B_PER_W);
	}
	free(mem);
}