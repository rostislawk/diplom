#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <time.h>

using namespace std;
typedef unsigned long uint64;
uint64 RDTSC()
{
    __asm _emit 0x0F __asm _emit 0x31
}

static const int _primes[] = { 3...8167};
//остальное вырезано, так как много лишнего

int gcdex (int a, int b, int & x, int & y) {
    if (a == 0) {
        x = 0; y = 1;
        return b;
    }
    int x1, y1;
    int d = gcdex (b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

int inverse(int a, int m) {
    int x, y;
    int g = gcdex (a, m, x, y);
    if (g != 1)
        return -1;
    else {
        x = (x % m + m) % m;
        return x;
    }
}

int inv2p(int p) {
    return (p + 1) / 2;

}

int powmod(unsigned long long a, unsigned long long pow, unsigned int mod) {
    unsigned long long ans = 1;
    int b = a % mod;
    while (pow != 0) {
        if (pow % 2 == 1) {
            ans *= b;
            ans %= mod;
        }
        b *= b;
        b %= mod;
        pow /= 2;
    }
    return ans;
}

bool miller_rabin (int n)
{
    // сначала проверяем тривиальные случаи
    if (n == 2)
        return true;
    if (n < 2 || n % 2 == 0)
        return false;

    // разлагаем n - 1 = 2^s * t
    int t = n - 1, s = 0;
    while (t % 2 == 0) {
        t /= 2;
        ++s;
    }
    int b[3] = {2, 7, 61};
    int rounds = 0;
    while (rounds < 3) {
        // вычисляем b^q mod n, если оно равно 1 или n-1, то n простое (или псевдопростое)
        int rem = powmod (b[rounds], t, n);
        if ((rem == 1 || rem == n - 1) && rounds == 1) {
            return true;
        }

        // теперь вычисляем b^2q, b^4q, ... , b^((n-1)/2)
        // если какое-либо из них равно n-1, то n простое (или псевдопростое)
        for (int i = 1; i < s; ++i)
        {
            rem = powmod (rem, 2, n);
            if (rem == n - 1 && rounds == 1) {
                return true;
            }
        }
        ++rounds;
    }
    return false;
}

vector<int> fast_check(int n, int step) {
    set<int> bolter;
    vector<int> answer;
    for (int i = 0; i < step; ++i) {
        int p = _primes[i];
        int r = p - n % p;
        int x = 0;
        if (r % 2 == 0)
            x += r;
        else 
            x += r + p;
        if (r == p) 
            x = 0;
        for (int j = x; j < _primes[step-1]; j = j + 2 * p) {
            bolter.insert(j);
        }
    }

    for (set<int>::iterator it = bolter.begin(); it != bolter.end(); ) {
        set<int>::iterator i1 = it, i2 = ++it;
        if (i2 != bolter.end()) {
            int rest = (*i2 - *i1) / 2 - 1;
            while (rest > 0) {
                answer.push_back(n + *i2 - 2 * rest);
                --rest;
            }
        }
    }

    return answer;
}

void func1(int count, int step) {
    srand(time(0));
    vector<int>* a = new vector<int>();
    do {
        do {
            int r = rand() % (INT_MAX - 65536) + 65536; 
            if (r % 2 == 0) 
                r += 1;
            *a = fast_check(r, step);
        } while (a->size() == 0);
        for (int i = 0; i < a->size(); ++i) {
            if (!miller_rabin(a->at(i)))	{
            //	cout << "prime number = " << a->at(i) << endl;
                count--;
            }
        }
    } while (count > 0);
    
}

void func2(int count) {
    int r = 0;
    do {
        do {
            r = rand () % (INT_MAX - 65536) + 65536;
            if (r % 2 == 0) 
                r += 1;
        } while (!miller_rabin(r));
        //cout << "prime number = " << r << endl;
    } while (--count > 0);
}

int main() {
    vector<uint64> f1, f2;
    fstream out;
    out.open("freq1.txt", fstream::out);
    uint64 start, overhead, freq1, freq2;
    for (int i = 2; i < 1024; i += 1) {
        start = RDTSC();
        overhead = RDTSC() - start;
        func1(1000, i);
        freq1 = RDTSC() - start - overhead;
        f1.push_back(freq1 / 1000);
        out << freq1 / 1000 << endl;
        cout << i << "-> OK " << endl;
    }
    start = RDTSC();
    func2(1000);
    freq2 = RDTSC() - start - overhead;
    cout << "Processor ticks, elapsed on  miller-rabin = " << freq2 / 1000 << endl;
    out << endl << freq2 / 1000;
    out.close();
    return 0;
}