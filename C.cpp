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
    Big_Int &operator*=(const unsigned int num);     // Умножение на короткое
    Big_Int &operator*=(const Big_Int &other); // Умножение на длинное
    Big_Int &operator/=(const int num);     // Деление на короткое
    Big_Int &operator/=(const Big_Int &other); // Деление на длинное
    Big_Int &operator%=(const Big_Int &other); // Остаток от деления на длинное
};

std::istream &operator>>(std::istream &, Big_Int &); // Ввод из потока
std::ostream &operator<<(std::ostream &, const Big_Int &); // Вывод в поток

Big_Int pow(Big_Int, Big_Int); // Возведение в степень
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
// возведение в степень

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
Big_Int &Big_Int::operator*=(const unsigned int num) {
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
int operator%(const Big_Int &a, const unsigned int num) {
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

// Возведение в степень
Big_Int p;

// Возведение в степень:
Big_Int pow(Big_Int a, Big_Int n) {
    Big_Int res = 1;
    while (n > 0) {
        if (n % 2 != 0) res *= a;
        a *= a;
        a %= p;
        n /= 2;
        res %= p;
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

class Polynomial {
private:
    void modp() {
        for (Big_Int &i: coefficients) {
            i %= p;
        }
    }

private:
    void remove_zeros_from_end() {
        size_t len = 0;
        for (long long i = coefficients.size() - 1; i > -1; --i) {
            if (coefficients[i] % p!= Big_Int(0)) {
                len = i + 1;
                break;
            }
        }
        coefficients.resize(len);
    }

public:
    std::vector<Big_Int> coefficients;

    Polynomial(std::vector<Big_Int> data) : coefficients(data) {
        remove_zeros_from_end();
    }

    Polynomial(const Big_Int &data = Big_Int()) {
        coefficients.push_back(data);
        remove_zeros_from_end();
    }

    template<typename It>
    Polynomial(It begin, It end) : coefficients(begin, end) {
        remove_zeros_from_end();
    }

    bool operator==(const Polynomial &other) const {
        return coefficients == other.coefficients;
    }

    bool operator!=(const Polynomial &other) const {
        return coefficients != other.coefficients;
    }

    Polynomial &operator+=(const Polynomial &other) {
        coefficients.resize(std::max(coefficients.size(), other.coefficients.size()));
        for (long long i = 0; i != other.coefficients.size(); ++i) {
            coefficients[i] += other.coefficients[i];
        }
        remove_zeros_from_end();
        return *this;
    }

    Polynomial operator+(const Polynomial &B) const {
        Polynomial result = *this;
        result += B;
        return result;
    }

    Polynomial operator-( Polynomial &B) const {
        Polynomial result = *this;
        result.remove_zeros_from_end();
        B.remove_zeros_from_end();
        for (long long i = 0; i!= B.coefficients.size(); ++i) {
            result.coefficients[i] = (result.coefficients[i] + p - B.coefficients[i]) % p;
        }
        result.remove_zeros_from_end();
        return result;
    }

    const Big_Int operator[](size_t i) const {
        if (i > coefficients.size() - 1) {
            return 0;
        }
        return coefficients[i];
    }

    size_t Degree() const {
        return coefficients.size();
    }

    typename std::vector<Big_Int>::const_iterator begin() const {
        return coefficients.cbegin();
    }

    typename std::vector<Big_Int>::const_iterator end() const {
        return coefficients.cend();
    }

    Polynomial &operator*=(const Polynomial &other) {
        std::vector<Big_Int> res(coefficients.size() + other.coefficients.size() - 1);
        for (long long i = 0; i != coefficients.size(); ++i) {
            for (long long j = 0; j != other.coefficients.size(); ++j) {
                res[i + j] = (res[i + j]  + (coefficients[i] * other.coefficients[j]) % p) % p;
            }
        }
        coefficients = res;
        return *this;
    }

    Polynomial operator*(const Polynomial &other) {
        Polynomial res = *this;
        res *= other;
        return res;
    }

    Polynomial operator/(const Polynomial& other) {
        Polynomial p1 = *this;
        Polynomial p2;
        Polynomial res;
        vector<Big_Int> vres(p1.coefficients.size(), 0);
        if (p1.Degree() < other.Degree()) {
            res = Polynomial(Big_Int(1));
            return res;
        }
        Big_Int multiplier1;
        Big_Int multiplier2;
        Polynomial p3;
        while (p1.Degree() >= other.Degree()) {
            p2 = other;
            multiplier1 = pow(p2.coefficients.back(), p - 2);
            multiplier2 = (p1.coefficients.back() * multiplier1) % p;
            vres[p1.Degree() - p2.Degree()] = multiplier2;
            p3 = (p2 * Polynomial(multiplier2));
            vector<Big_Int> x(p1.Degree() - other.Degree() + 1);
            x.back() = 1;
            p3 = p3 * Polynomial(x);
            p1 = p1 - p3;
        }
        res.coefficients = vres;
        res.remove_zeros_from_end();
        return res;
    }

    Polynomial operator %(const Polynomial& other) {
        Polynomial p1 = *this;
        p1.remove_zeros_from_end();
        Polynomial ans;
        if (p1.Degree() < other.Degree()) {
            return  p1;
        }
        Polynomial p3 = (p1 / other) * other;
        ans = p1 - p3;
        ans.remove_zeros_from_end();
        return ans;
    }

    friend std::ostream &operator<<(std::ostream &out, const Polynomial &p);
};

std::ostream &operator<<(std::ostream &out, const Polynomial &p) {
    for (const Big_Int & coefficient : p.coefficients) {
        out << coefficient << ' ';
    }
    return out;
}

Big_Int char_to_number(char symbol) {
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

void string_to_coeffs (string& str, vector<Big_Int>& coeffs) {
    Big_Int num_in_10 = 0;
    bool is_minus = 0;
    for (char i : str) {
        if (i == '-') {
            is_minus = 1;
        } else if (i == ' ') {
            if (is_minus) {
                coeffs.push_back(p - num_in_10);
            } else {
                coeffs.push_back(num_in_10);
            }
            num_in_10 = 0;
            is_minus = 0;
        } else {
            num_in_10 *= 10;
            num_in_10 += Big_Int(i - '0');
        }
    }
    if (is_minus)
        coeffs.push_back(p - num_in_10);
    else
        coeffs.push_back(num_in_10);
}

int main() {
    string trash;
    cin >> p;
    getline(cin, trash);
    vector<Big_Int> f_poly;
    vector<Big_Int> g_poly;
    vector<Big_Int> k_poly;
    string f_coef;
    getline(cin, f_coef);
    string_to_coeffs(f_coef, f_poly);
    string g_coef;
    getline(cin, g_coef);
    string_to_coeffs(g_coef, g_poly);
    string k_coef;
    getline(cin, k_coef);
    string_to_coeffs(k_coef, k_poly);
    string message;
    getline(cin, message);
    Big_Int n = f_poly.size() - 1;
    Polynomial F(f_poly);
    Polynomial G(g_poly);
    Polynomial K(k_poly);
    Big_Int data = 0;
    Big_Int cnt = 1;
    for (char i: message) {
        data += cnt * Big_Int(char_to_number(i));
        cnt *= 64;
    }
    vector<Big_Int> num;
    while (data != 0) {
        num.push_back(data % p);
        data /= p;
    }
    vector<vector<Big_Int>> polynomials;
    for (unsigned long long i = 0; i != num.size(); ++i) {
        if (!(i % (f_poly.size() - 1))) {
            polynomials.emplace_back();
        }
        polynomials.back().push_back(num[i]);
    }
    for (const vector<Big_Int>& pol : polynomials) {
        Polynomial pp(pol);
        Polynomial ans1 = G * G;
        ans1 = ans1 % F;
        cout << ans1 << '\n';
        Polynomial ans2 = K * K;
        ans2 = ans2 * pp;
        ans2 = ans2 % F;
        cout << ans2 << '\n';
    }
}
