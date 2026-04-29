/*
 * Math.h
 *
 *  Created on: Jun 20, 2023
 *      Author: sergey
 */

#ifndef COMMON_MATH_H_
#define COMMON_MATH_H_

#include <cmath>
#include <utility>
#include <algorithm>
#include <limits>
#include <type_traits>

namespace Skasp {

//
//  Константы
//

#define M_TAU          0.61803398874989484820  /* tau */

//
//  Определения
//

template <typename T> using RangeT = std::pair<T, T>;
using Range  = RangeT<float>;
using RangeD = RangeT<double>;

//
//  Математические функции
//

//! Модуль числа
/**
 * @note Для тех, кто не хочет использовать std::abs()
 * @param __val число
 * @return значение по модулю
 */
template<class X> X ABS(const X& __val) { return __val < 0 ? -__val : __val; }

//! Знак числа в формате int
/**
 * @param val число
 * @return +1 если число положительное, -1 если число отрицательное, 0 если число равно нулю
 */
template <typename T>
int sgn(T val) { return (val == T(0)) ? 0 : val < T(0) ? -1 : 1; }

//! Квадрат числа
/**
 * @param a число
 * @return число в квадрате
 */
template<typename T>
inline T SQ(T a) { return a*a; }

//! Куб числа
/**
 * @param a число
 * @return число в кубе
 */
template<typename T>
inline T CUB(T a) { return a*a*a; }

//! Среднее значение двух переменных
/**
 * @details Переменнные должны поддерживать сложение и умножение на число (это могут быть, например, числа, вектора или матрицы)
 * @param a Первая переменная
 * @param b Вторая переменная
 * @return Среднее значение
 */
template <typename T>
T avg(T a, T b) { return 0.5*(a + b); }

//! Линейная интерполяция двух чисел
/**
 * @details a и b - вещественные числа, t должен быть внутри интервала [0, 1).
 * @param a Первое число
 * @param b Второе число
 * @param t Параметр интероляции
 * @return Результат интероляции
 */
template <typename T>
T Lerp( T a, T b, T t ) { return a + t*(b - a); }

//! Факториал для числа, ограниченного 32-битной переменной
/**
 * @param k Число
 * @return Факториал числа
 */
unsigned int factorial(unsigned int k);

//
// Сравнения
//

//! Шаблон для параметров вещественных чисел
template <typename T> struct fp_limits {};
template <> struct fp_limits<float>  { static float eps()  { return 1.0e-05f; } };
template <> struct fp_limits<double> { static double eps() { return 1.0e-14;  } };

//! Корректное сравнение двух вещественных чисел с учётом дискретного представления
/**
 * @note "Comparing Two Floating-Point Numbers" by Burkhard Stubert, 2019/08/26, https://embeddeduse.com/2019/08/26/qt-compare-two-floats
 * @param r1 Первое значение
 * @param r2 Второе значение
 * @return Результат сравнения
 */
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type fp_compare(T r1, T r2)
{
    if (r1 == r2)
        return true;

    T diff = std::abs(r1 - r2);
    return diff <= fp_limits<T>::eps() ? true :
           diff <= fp_limits<T>::eps()*std::max(std::abs(r1), std::abs(r2));
}

//! Прямое сравнение двух переменных
/**
 * @param __a Первое значение
 * @param __b Второе значение
 * @return Результат сравнения
 */
template <typename T>
struct exact_equal {
    static bool is_equal(T __a, T __b) { return __a == __b; }
};

//! Сравнение двух переменных через модуль разности
/**
 * @param __a Первое значение
 * @param __b Второе значение
 * @return Результат сравнения
 */
template <typename T, typename Enable = void>
struct inexact_equal {
    static bool is_equal(T __a, T __b) { return std::abs(__a - __b) < std::numeric_limits<T>::epsilon(); }
};

//! Сравнение двух переменных для float и double
/**
 * @param __a Первое значение
 * @param __b Второе значение
 * @return Результат сравнения
 */
template <typename T>
struct inexact_equal<T, typename std::enable_if<std::is_floating_point<T>::value>::type>  {
    static bool is_equal(T __a, T __b) { return fp_compare(__a, __b); }
};

//! Сравнение двух переменных c автоматическим выбором метода сравнения
/**
 * @param __a Первое значение
 * @param __b Второе значение
 * @return Результат сравнения
 */
template <typename T>
bool equal(T __a, T __b)
{
    return std::conditional<std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_exact,
           inexact_equal<T>, exact_equal<T> >::type::is_equal(__a, __b);
}

//! Сравнение переменной с нулём
/**
 * @param __a Значение
 * @return Результат сравнения
 */
template <typename T>
bool IsZero(T __a) {
    return equal(__a, T(0));
}

//! Проверка возможности деления на данное число
/**
 * @param x Число
 * @return Результат проверки
 */
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type allow_division(T x)
{
    return std::abs(x) > fp_limits<T>::eps();
}

template <typename T>
typename std::enable_if<std::is_integral<T>::value && !std::is_unsigned<T>::value, bool>::type allow_division(T x)
{
    return std::abs(x) > T(0);
}

template <typename T>
typename std::enable_if<std::is_integral<T>::value && std::is_unsigned<T>::value, bool>::type allow_division(T x)
{
    return x > T(0);
}

} // namespace Skasp {

#endif /* COMMON_MATH_H_ */
