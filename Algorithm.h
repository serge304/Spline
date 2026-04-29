/*
 * Algorithm.h
 *
 *  Created on: Jun 20, 2023
 *      Author: sergey
 */

#ifndef COMMON_ALGORITHM_H_
#define COMMON_ALGORITHM_H_

#include <cstdint>
#include <utility>
#include <functional>

#include <type_traits>
#include <algorithm>

namespace Skasp {

//
//  Численные функции
//

//! Фильтр ограничивает число неотрицательными значениями
/** Неположительные значения заменяются на ноль
 * @param val Число для фильтрации
 * @return Результат фильтра
 */
template <typename T>
T positive_filter(T val) { return val > T(0) ? val : T(0); }

//! Проверяет, что число находится в заданном полуоткрытом интервале
/**
 * \f$x \in [a, b)\f$
 * @param x Число для проверки
 * @param a Левая граница интервала
 * @param b Правая граница интервала
 * @return Результат проверки
 */
template <typename T, typename U, typename W>
bool in_range(T x, U a, W b) { return x >= static_cast<T>(a) && x < static_cast<T>(b); }

//! Проверяет, что число не находится в заданном полуоткрытом интервале
/**
 * \f$x \notin [a, b)\f$
 * @param x Число для проверки
 * @param a Левая граница интервала
 * @param b Правая граница интервала
 * @return Результат проверки
 */
template <typename T, typename U, typename W>
bool not_in_range(T x, U a, W b) { return !in_range(x, a, b); }

//! Проверяет, что число находится в заданном полуоткрытом интервале
/**
 * \f$x \in r\f$
 * @param x Число для проверки
 * @param r Интервал в формате std::pair
 * @return Результат проверки
 */
template <typename T, typename U>
bool in_range(T x, const std::pair<U, U>& r) { return in_range(x, r.first, r.second); }

//! Проверяет, что число не находится в заданном полуоткрытом интервале
/**
 * \f$x \notin r\f$
 * @param x Число для проверки
 * @param r Интервал в формате std::pair
 * @return Результат проверки
 */
template <typename T, typename U>
bool not_in_range(T x, const std::pair<U, U>& r) { return !in_range(x, r); }

//! Проверяет, что число находится в заданном закрытом интервале
/**
 * \f$x \in [a, b]\f$
 * @param x Число для проверки
 * @param a Левая граница интервала
 * @param b Правая граница интервала
 * @return Результат проверки
 */
template <typename T, typename U, typename W>
bool in_range_inc(T x, U a, W b) { return x >= static_cast<T>(a) && x <= static_cast<T>(b); }

//! Ограничивает число сверху
/**
 * Если число выше максимального значения, возвращается максимальное значение
 * @param x Число для проверки
 * @param b Максимальное значение
 * @return Результат ограничения
 */
template <typename T, typename W>
T limit(T x, W b) { return x > b ? b : x; }

//! Ограничивает число заданным закрытым интервалом
/**
 * Если число вне интервала, возвращается соответствующая граница интервала
 * @param x Число для проверки
 * @param a Левая граница интервала
 * @param b Правая граница интервала
 * @return Результат ограничения
 */
template <typename T, typename U, typename W>
T crop(T x, U a, W b) { return x < a ? a : x > b ? b : x; }  // std::clamp

//! Ограничивает число заданным закрытым интервалом
/**
 * @param x Переменная для проверки
 * @param r Интервал в формате std::pair
 * @return Результат ограничения
 */
template <typename T, typename U>
T crop(T x, const std::pair<U, U>& r) { return crop(x, r.first, r.second); }

//! Округление вещественного числа до одной цифры после запятой
/**
 * @param x Исходное значение
 * @return Результат округления
 */
double round_one_dp(double x);

//! Округление вещественного числа до двух цифр после запятой
/**
 * @param x Исходное значение
 * @return Результат округления
 */
double round_two_dp(double x);

//! Гарантированное преобразования целого числа в bool
inline bool to_bool(uint8_t val) { return val > 0; }

//! Гарантированное преобразования bool в целое число
inline uint8_t from_bool(bool val) { return val ? 1 : 0; }

//
// Функторы
//

template <typename T>
auto not_in_range_fn(T _min, T _max) {
    return [_min = std::move(_min), _max = std::move(_max)](const auto& val) {
        return not_in_range(val, _min, _max); 
    };
}

template <typename T>
auto not_in_range_fn(const std::pair<T, T>& _range) {
    return not_in_range_fn(_range.first, _range.second);
}

//
// Контейнеры
//

//! Добавление элементов одного контейнера в другой контейнер
/**
 * @param A Контейнер - получатель элементов
 * @param B Контейнер - источник элементов
 */
template<typename A, typename B>
void append(A& receiver, const B& provider)
{
  receiver.insert(receiver.end(), provider.begin(), provider.end());
}

//
// Получение значения из карты
//

template<class M>
typename M::mapped_type FindInMap(
    const M& m,
    const typename M::key_type key1,
    typename M::mapped_type defaultValue)
{
  typename M::const_iterator it1 = m.find(key1);
  return it1 != m.end() ? it1->second : defaultValue;
}

template<class M>
bool TryFindInMap(
    const M& m,
    const typename M::key_type key1,
    typename M::mapped_type& value)
{
  typename M::const_iterator it1 = m.find(key1);
  return it1 != m.end() ? (value = it1->second, true) : false;
}

template<class M>
typename M::mapped_type::mapped_type FindInMap(
    const M& m,
    const typename M::key_type key2,
    const typename M::mapped_type::key_type key1,
    typename M::mapped_type::mapped_type defaultValue)
{
  typename M::const_iterator it = m.find(key2);
  return it != m.end() ? FindInMap(it->second, key1, defaultValue) : defaultValue;
}

template<class M>
bool TryFindInMap(
    const M& m,
    const typename M::key_type key2,
    const typename M::mapped_type::key_type key1,
    typename M::mapped_type::mapped_type& value)
{
  typename M::const_iterator it = m.find(key2);
  return it != m.end() ? TryFindInMap(it->second, key1, value) : false;
}

template<class M>
typename M::mapped_type::mapped_type::mapped_type FindInMap(
    const M& m,
    const typename M::key_type key3,
    const typename M::mapped_type::key_type key2,
    const typename M::mapped_type::mapped_type::key_type key1,
    typename M::mapped_type::mapped_type::mapped_type defaultValue)
{
  typename M::const_iterator it = m.find(key3);
  return it != m.end() ? FindInMap(it->second, key2, key1, defaultValue) : defaultValue;
}

//
// Пары
//

//! Инкремент обоих элментов пары
/**
 * @param p Пара
 */
template <typename first_type, typename second_type>
void operator++(std::pair<first_type, second_type> &p)
{
  ++p.first;
  ++p.second;
}

//! Декремент обоих элментов пары
/**
 * @param p Пара
 */
template <typename first_type, typename second_type>
void operator--(std::pair<first_type, second_type> &p)
{
  --p.first;
  --p.second;
}

//
//  Утилиты
//

#define STATIC_ASSERT(x) static_assert(x, #x)

//
// Преобразование типов
//

template<typename T, typename F>
struct alias_cast_t
{
    union
    {
        F raw;
        T data;
    };
};

template<typename T, typename F>
T alias_cast(F raw_data)
{
	STATIC_ASSERT(sizeof(T) == sizeof(F));

    alias_cast_t<T, F> ac;
    ac.raw = raw_data;
    return ac.data;
}

} // namespace Skasp {

#endif /* COMMON_ALGORITHM_H_ */
