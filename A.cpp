#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cassert>
#include <functional>
#include <complex>

using namespace std;

const long double PI = std::acos(-1.0L);

int char_to_number(char symbol) {
    if (symbol >= 48 && symbol <= 57)
        return symbol - 48;
    if (symbol >= 65 && symbol <= 90)
        return symbol - 55;
    if (symbol >= 97 && symbol <= 122)
        return symbol - 61;
    if (symbol == 32)
        return 62;
    if (symbol == 46)
        return 63;
    return 64;
}

struct Big_Int {
    static const int BASE = (int) 1e9; // Основание системы счисления
    static const int WIDTH = 9;       // Количество десятичных цифр, которые хранятся в одной цифре

    // Вектор под цифры числа:
    std::vector<int> digits;

    // Конструкторы
    Big_Int(int64_t number = 0);

    Big_Int(const std::string &s);

    Big_Int(const std::vector<int> &digits);

    // Методы нормализации и сравнения:
    Big_Int &normalize(); // удаление лидирующих нулей и проверка на принадлежность цифр диапазону [0, BASE)
    int compare(const Big_Int &other) const; // Сравнение (меньше = -1, равно = 0, больше = 1)

    // Методы умножения:
    Big_Int
    slow_mult(
            const Big_Int &other) const; // Медленное произведение (работает довольно быстро на числах небольшой длины)
    Big_Int fast_mult(
            const Big_Int &other) const; // Быстрое произведение (на основе Быстрого Преобразования Фурье комплексные числа)
    Big_Int mult(const Big_Int &other) const; // Комбинированный метод умножения на основе экспериментальных данных

    // Метод деления:
    std::pair<Big_Int, Big_Int> div_mod(const Big_Int &other) const; // Целая часть и остаток от деления

    // Операторы:
    Big_Int &operator+=(const int num);     // Прибавление короткого
    Big_Int &operator+=(const Big_Int &other); // Прибавление длинного
    Big_Int &operator-=(const int num);     // Вычитание короткого
    Big_Int &operator-=(const Big_Int &other); // Вычитание длинного
    Big_Int &operator*=(const int num);     // Умножение на короткое
    Big_Int &operator*=(const Big_Int &other); // Умножение на длинное
    Big_Int &operator/=(const int num);     // Деление на короткое
    Big_Int &operator/=(const Big_Int &other); // Деление на длинное
    Big_Int &operator%=(const Big_Int &other); // Остаток от деления на длинное
};

std::istream &operator>>(std::istream &, Big_Int &); // Ввод из потока
std::ostream &operator<<(std::ostream &, const Big_Int &); // Вывод в поток

Big_Int pow(Big_Int, int); // Возведение в степень
Big_Int gcd(Big_Int, Big_Int); // Наибольший общий делитель

Big_Int operator+(const Big_Int &, const Big_Int &);

Big_Int operator-(const Big_Int &, const Big_Int &);

Big_Int operator*(const Big_Int &, const Big_Int &);

Big_Int operator/(const Big_Int &, const Big_Int &);

Big_Int operator%(const Big_Int &, const Big_Int &);

Big_Int operator+(const Big_Int &, const int);

Big_Int operator+(const int, const Big_Int &);

Big_Int operator-(const Big_Int &, const int);

Big_Int operator*(const Big_Int &, const int);

Big_Int operator*(const int, const Big_Int &);

Big_Int operator/(const Big_Int &, const int);

Big_Int operator^(const Big_Int &, const int); // возведение в степень

bool operator<(const Big_Int &, const Big_Int &);

bool operator>(const Big_Int &, const Big_Int &);

bool operator<=(const Big_Int &, const Big_Int &);

bool operator>=(const Big_Int &, const Big_Int &);

bool operator==(const Big_Int &, const Big_Int &);

bool operator!=(const Big_Int &, const Big_Int &);

Big_Int &Big_Int::normalize() {
    while (digits.back() == 0 && (int) digits.size() > 1) digits.pop_back();
    for (auto d: digits) assert(0 <= d && d < BASE);
    return *this;
}

// Конструктор от короткого целого
Big_Int::Big_Int(int64_t number) {
    assert(number >= 0);
    do {
        digits.push_back(number % BASE);
        number /= BASE;
    } while (number > 0);
    normalize();
}

// Конструктор от вектора из цифр:
Big_Int::Big_Int(const std::vector<int> &digits) : digits(digits) {
    normalize();
}

// Конструктор от строчки:
Big_Int::Big_Int(const std::string &s) {
    const int size = (int) s.size();
    for (int idGroup = 1, nGroups = size / WIDTH; idGroup <= nGroups; ++idGroup) {
        digits.push_back(std::stoi(s.substr(size - idGroup * WIDTH, WIDTH)));
    }
    if (size % WIDTH != 0) {
        digits.push_back(std::stoi(s.substr(0, size % WIDTH)));
    }
    normalize();
}

// Прибавление короткого:
Big_Int &Big_Int::operator+=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this += Big_Int(num);
    }
    int rem = num;
    for (int i = 0; rem > 0; ++i) {
        if (i >= (int) digits.size()) digits.push_back(0);
        rem += digits[i];
        if (rem >= BASE) {
            digits[i] = rem - BASE;
            rem = 1;
        } else {
            digits[i] = rem;
            rem = 0;
        }
    }
    return this->normalize();
}

// Прибавление длинного:
Big_Int &Big_Int::operator+=(const Big_Int &other) {
    if (other.digits.size() == 1u) {
        return *this += other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    int rem = 0;
    for (int i = 0; i < s1 || i < s2 || rem > 0; ++i) {
        int d1 = i < s1 ? this->digits[i] : (digits.push_back(0), 0);
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += d1 + d2;
        auto div = rem / BASE;
        digits[i] = rem - div * BASE;
        rem = div;
    }
    return this->normalize();
}

// Вычитание короткого:
Big_Int &Big_Int::operator-=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this -= Big_Int(num);
    }
    int rem = -num;
    for (int i = 0; i < (int) digits.size() && rem < 0; ++i) {
        rem += digits[i];
        if (rem < 0) { // Занимаем разряд
            digits[i] = rem + BASE;
            rem = -1;
        } else {
            digits[i] = rem;
            rem = 0;
        }
    }
    assert(rem == 0);
    return this->normalize();
}

// Вычитание длинного:
Big_Int &Big_Int::operator-=(const Big_Int &other) {
    if (other.digits.size() == 1u) {
        return *this -= other.digits[0];
    }
    const int s1 = this->digits.size();
    const int s2 = other.digits.size();
    assert(s1 >= s2);
    int rem = 0;
    for (int i = 0; i < s1; ++i) {
        int d2 = i < s2 ? other.digits[i] : 0;
        rem += this->digits[i] - d2;
        if (rem < 0) {
            digits[i] = rem + BASE;
            rem = -1;
        } else {
            digits[i] = rem;
            rem = 0;
            if (i >= s2) break;
        }
    }
    assert(rem == 0); // Иначе *this < other
    return this->normalize();
}

// Умножение на короткое:
Big_Int &Big_Int::operator*=(const int num) {
    assert(num >= 0);
    if (num >= BASE) {
        return *this *= Big_Int(num);
    }
    int64_t rem = 0;
    for (auto &d: digits) {
        rem += 1LL * d * num;
        auto div = rem / BASE;
        d = rem - div * BASE;
        rem = div;
    }
    if (rem > 0) digits.push_back(rem);
    return this->normalize();
}

// Медленное произведение:
Big_Int Big_Int::slow_mult(const Big_Int &other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }
    const int s1 = (int) this->digits.size();
    const int s2 = (int) other.digits.size();
    std::vector<int> temp(s1 + s2);
    for (int i = 0; i < s1; ++i) {
        int64_t rem = 0;
        for (int j = 0; j < s2; ++j) {
            rem += temp[i + j] + 1LL * this->digits[i] * other.digits[j];
            auto div = rem / BASE;
            temp[i + j] = rem - div * BASE;
            rem = div;
        }
        if (rem > 0) temp[i + s2] += rem;
    }
    return Big_Int(temp);
}

// Быстрое умножение на основе быстрого преобразования Фурье:
Big_Int Big_Int::fast_mult(const Big_Int &other) const {
    if (other.digits.size() == 1u) {
        return *this * other.digits[0];
    }

    // Разворот битов в числе num:
    std::function<int(int, int)> reverse = [](int number, int nBits) {
        int res = 0;
        for (int i = 0; i < nBits; ++i) {
            if (number & (1 << i)) {
                res |= 1 << (nBits - 1 - i);
            }
        }
        return res;
    };

    typedef std::complex<long double> complex;
    // Быстрое преобразование Фурье:
    std::function<void(std::vector<complex> &, bool)> fft = [&reverse](std::vector<complex> &a, bool invert) {
        const int n = (int) a.size();
        int nBits = 0;
        while ((1 << nBits) < n) ++nBits;

        for (int i = 0; i < n; ++i) {
            if (i < reverse(i, nBits)) {
                std::swap(a[i], a[reverse(i, nBits)]);
            }
        }

        for (int len = 2; len <= n; len <<= 1) {
            auto ang = 2 * PI / len * (invert ? -1 : 1);
            complex wlen(std::cos(ang), std::sin(ang));
            for (int i = 0; i < n; i += len) {
                complex w(1);
                for (int j = 0; j < len / 2; ++j) {
                    complex u = a[i + j];
                    complex v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (int i = 0; i < n; ++i) {
                a[i] /= n;
            }
        }
    };

    // Подготавливаем вектора из комплексных коэффициентов fa и fb:
    // Так как происходит потеря точности из-за арифметики с плавающей точкой, основание системы необходимо понизить:
    assert(BASE == 1000 * 1000 * 1000);
    std::function<std::vector<complex>(const Big_Int &)> prepare = [](const Big_Int &number) {
        std::vector<complex> result;
        result.reserve(3 * number.digits.size());
        for (auto d: number.digits) {
            result.push_back(d % 1000);
            result.push_back(d / 1000 % 1000);
            result.push_back(d / 1000000);
        }
        return result;
    };

    auto fa = prepare(*this);
    auto fb = prepare(other);

    // Округляем размер векторов до ближайшей степени двойки:
    int n = 1;
    while (n < (int) std::max(fa.size(), fb.size())) n *= 2;
    n *= 2;
    fa.resize(n);
    fb.resize(n);

    // Вызываем прямое преобразование Фурье:
    fft(fa, false);
    fft(fb, false);
    // Перемножаем результаты:
    for (int i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }
    // Вызываем обратное преобразование Фурье:
    fft(fa, true);
    // Копируем ответ с округлениями:
    std::vector<int64_t> temp(n);
    for (int i = 0; i < (int) fa.size(); ++i) {
        temp[i] = int64_t(fa[i].real() + 0.5);
    }
    // Не забываем про переносы в старшие разряды:
    int64_t carry = 0;
    for (int i = 0; i < n || carry > 0; ++i) {
        if (i >= n) temp.push_back(0);
        temp[i] += carry;
        carry = temp[i] / 1000;
        temp[i] -= carry * 1000;
        assert(temp[i] >= 0);
    }
    // Формируем ответ:
    std::vector<int> res;
    res.reserve(this->digits.size() + other.digits.size());

    for (int i = 0; i < n; i += 3) {
        int c = temp[i];
        int b = i + 1 < n ? temp[i + 1] : 0;
        int a = i + 2 < n ? temp[i + 2] : 0;
        res.push_back(c + 1000 * (b + 1000 * a));
    }
    return Big_Int(res);
}

// Комбинированный метод умножения:
Big_Int Big_Int::mult(const Big_Int &other) const {
// Выбор метода умножения:
    int len1 = (int) this->digits.size();
    int len2 = (int) other.digits.size();
    int temp = 3 * std::max(len1, len2);
    int pow = 1;
    while (pow < temp) pow *= 2;
    pow *= 2;
    int op1 = len1 * len2;
    int op2 = 3 * pow * std::log(pow) / std::log(2);
    return op1 >= 15 * op2 ? fast_mult(other) : slow_mult(other);
}

// Деление на короткое:
Big_Int &Big_Int::operator/=(const int num) {
    assert(num > 0);
    if (num >= BASE) {
        return *this /= Big_Int(num);
    }
    int64_t rem = 0;
    for (int j = (int) digits.size() - 1; j >= 0; --j) {
        (rem *= BASE) += digits[j];
        auto div = rem / num;
        digits[j] = div;
        rem -= div * num;
    }
    return this->normalize();
}

// Остаток от деления на короткое:
int operator%(const Big_Int &a, const int num) {
    assert(num > 0);
    int64_t rem = 0;
    for (int i = (int) a.digits.size() - 1; i >= 0; --i) {
        ((rem *= Big_Int::BASE) += a.digits[i]) %= num;
    }
    return rem;
}

// Целая часть и остаток от деления:
std::pair<Big_Int, Big_Int> Big_Int::div_mod(const Big_Int &other) const {
    if (other.digits.size() == 1u) {
        return {std::move(*this / other.digits[0]), *this % other.digits[0]};
    }
    const int norm = BASE / (other.digits.back() + 1);
    const Big_Int a = *this * norm;
    const Big_Int b = other * norm;
    const int a_size = (int) a.digits.size();
    const int b_size = (int) b.digits.size();
    Big_Int q, r;
    q.digits.resize(a_size);
    for (int i = a_size - 1; i >= 0; --i) {
        r *= BASE;
        r += a.digits[i];
        int s1 = (int) r.digits.size() <= b_size ? 0 : r.digits[b_size];
        int s2 = (int) r.digits.size() <= b_size - 1 ? 0 : r.digits[b_size - 1];
        int d = (1LL * BASE * s1 + s2) / b.digits.back();
        auto temp = b * d;
        while (r < temp) {
            r += b;
            --d;
        }
        r -= temp;
        q.digits[i] = d;
    }
    return {std::move(q.normalize()), std::move(r /= norm)};
}

// Сравнение: result < 0 (меньше), result == 0 (равно), result > 0 (больше)
int Big_Int::compare(const Big_Int &other) const {
    if (this->digits.size() > other.digits.size()) return 1;
    if (this->digits.size() < other.digits.size()) return -1;
    for (int i = (int) digits.size() - 1; i >= 0; --i) {
        if (this->digits[i] > other.digits[i]) return 1;
        if (this->digits[i] < other.digits[i]) return -1;
    }
    return 0;
}

// Операторы сравнения:
bool operator<(const Big_Int &a, const Big_Int &b) { return a.compare(b) < 0; }

bool operator>(const Big_Int &a, const Big_Int &b) { return a.compare(b) > 0; }

bool operator==(const Big_Int &a, const Big_Int &b) { return a.compare(b) == 0; }

bool operator<=(const Big_Int &a, const Big_Int &b) { return a.compare(b) <= 0; }

bool operator>=(const Big_Int &a, const Big_Int &b) { return a.compare(b) >= 0; }

bool operator!=(const Big_Int &a, const Big_Int &b) { return a.compare(b) != 0; }

// Ввод из потока:
std::istream &operator>>(std::istream &is, Big_Int &number) {
    std::string s;
    is >> s;
    number = Big_Int(s);
    return is;
}

// Вывод в поток:
std::ostream &operator<<(std::ostream &os, const Big_Int &number) {
    os << number.digits.back();
    for (int i = (int) number.digits.size() - 2; i >= 0; --i) {
        os << std::setw(Big_Int::WIDTH) << std::setfill('0') << number.digits[i];
    }
    return os << std::setfill(' ');
}

// Сумма:
Big_Int operator+(const Big_Int &a, const Big_Int &b) {
    return Big_Int(a) += b;
}

// Разность:
Big_Int operator-(const Big_Int &a, const Big_Int &b) {
    return Big_Int(a) -= b;
}

// Произведение:
Big_Int operator*(const Big_Int &a, const Big_Int &b) {
    return a.mult(b);
}

// Деление:
Big_Int operator/(const Big_Int &a, const Big_Int &b) {
    return a.div_mod(b).first;
}

// Взятие остатка:
Big_Int operator%(const Big_Int &a, const Big_Int &b) {
    return a.div_mod(b).second;
}

// Умножение:
Big_Int &Big_Int::operator*=(const Big_Int &other) {
    return *this = *this * other;
}

// Деление с присваиванием:
Big_Int &Big_Int::operator/=(const Big_Int &other) {
    return *this = *this / other;
}

// Взятие остатка с присваиванием:
Big_Int &Big_Int::operator%=(const Big_Int &other) {
    return *this = *this % other;
}

Big_Int operator+(const Big_Int &a, const int b) { return Big_Int(a) += b; }

Big_Int operator+(const int a, const Big_Int &b) { return b * a; }

Big_Int operator-(const Big_Int &a, const int b) { return Big_Int(a) -= b; }

Big_Int operator*(const Big_Int &a, const int b) { return Big_Int(a) *= b; }

Big_Int operator*(const int a, const Big_Int &b) { return b * a; }

Big_Int operator/(const Big_Int &a, const int b) { return Big_Int(a) /= b; }

Big_Int operator^(const Big_Int &a, const int n) { return pow(a, n); } // Возведение в степень

// Возведение в степень:
Big_Int pow(Big_Int a, int n) {
    Big_Int res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        n /= 2;
    }
    return res;
}

// Наибольший общий делитель:
Big_Int gcd(Big_Int a, Big_Int b) {
    while (b != 0) {
        auto rem = a % b;
        a = b;
        b = rem;
    }
    return a;
}

int main() {
    Big_Int p, g, K;
    cin >> p >> g >> K;
    string trash;
    string str;
    getline(cin, trash);
    getline(cin, str);
    string snum;
    Big_Int data = 0;
    Big_Int cnt = 1;
    for (char i: str) {
        data += cnt * char_to_number(i);
        cnt *= 64;
    }
    vector<Big_Int> num;
    while (data != 0) {
        num.push_back(data % p);
        data /= p;
    }
    int k = 2;
    for (const Big_Int &i: num) {
        Big_Int ans1 = (g ^ k) % p;
        Big_Int ans2 = ((K ^ k) * i) % p;
        cout << ans1 << ' ' << ans2 << '\n';
    }
}
