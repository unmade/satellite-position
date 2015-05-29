/** \file date_converters.h
 * Файл определяет различные функции для перехода от
 * одного представления даты к другому.
 */

#ifndef SATELLITE_PROPAGATION_DATE_CONVERTERS_H
#define SATELLITE_PROPAGATION_DATE_CONVERTERS_H

#endif //SATELLITE_PROPAGATION_DATE_CONVERTERS_H

/**
 * Переводит календарную дату в модифицированный юлианский день
 *
 * @param[in] nyear Год
 * @param[in] nmonth Месяц
 * @param[in] nday День
 * @param[in] nhour Часы
 * @param[in] nminute Минуты
 * @param[in] nsecond Секунды
 *
 * @return Модифицированный юлианский день
 */
long double utc_to_mjd(int nyear, int nmonth, int nday,
                  int nhour, int nminute, long double nsecond);


/**
 * Переводит модифицированный юлианский день в календарную дату
 *
 * @param[in] mjd Модифицированный юлианский день
 * @param[out] nyear Год
 * @param[out] nmonth Месяц
 * @param[out] nday День
 * @param[out] nhour Часы
 * @param[out] nminute Минуты
 * @param[out] nsecond Секунды
 *
 * @return 0
 */
void mjd_to_utc(long double mjd, int *nyear, int *nmonth, int *nday,
                int *nhour, int *nminute, long double *nsecond);


/**
 * Переводит календарную дату в Земное время
 *
 * @param[in] nyear Год
 * @param[in] nmonth Месяц
 * @param[in] nday День
 * @param[in] nhour Часы
 * @param[in] nminute Минуты
 * @param[in] nsecond Секунды
 *
 * @return Земное время в юлианских днях
 */
long double utc_to_tt(int nyear, int nmonth, int nday,
                 int nhour, int nminute, long double nsecond);


/**
 * Переводит календарную дату в Земное время
 *
 * @param[in] mjd Момент времени UTC, выраженный в модифицированных юлианских днях
 *
 * @return Земное время в юлианских днях
 */
long double mjd_to_tt(long double mjd_in_utc);


/**
 * Переводит Земное время в Барицентрическое динамическое
 *
 * @param[in] tt Земное время
 *
 * @return Барицентрическое динамическое время
 */
long double tt_to_tdb(long double tt);


/**
 * Переводит момент времени UTC, выраженный в модифицированных в юлианских днях
 * в гринвическое среднее звездное время
 *
 * @param[in] utc_in_mjd момент времени UTC, выраженный в модифицированных в юлианских днях
 * @param[out] delta_ut Поправка ∆T = UT1 - UTC
 *
 * @return Гринвическое среднее звездное время выраженное в радианах
 */
long double utc_to_gmst(long double utc_in_mjd, long double delta_ut);


/**
 * Переводит момент времени UTC, выраженный в модифицированных в юлианских днях
 * в гринвическое истинное звездное время
 *
 * @param[in] year Год
 * @param[in] month Месяц
 * @param[in] day День
 * @param[in] hour Часы
 * @param[in] minute Минуты
 * @param[in] second Секунды
 * @param[out] delta_ut Поправка ∆T = UT1 - UTC
 *
 * @return Гринвическое истинное звездное время выраженное в радианах
 */
long double utc_to_gast(int year, int month, int day, int hour,
                   int minute, long double second, long double delta_ut);


/**
 * Переводит момент времени UTC, выраженный в модифицированных в юлианских днях
 * в гринвическое истинное звездное время
 *
 * @param[in] utc_in_mjd момент времени UTC, выраженный в модифицированных юлианских днях
 * @param[out] delta_ut Поправка ∆T = UT1 - UTC
 *
 * @return Гринвическое среднее звездное время выраженное в радианах
 */
long double mjd_to_gast(long double utc_in_mjd, long double delta_ut);


/**
 * Извлекает из дробной части дня часы, минуты и секунды
 *
 * @param[in] days Дробная часть дня
 * @param[out] nhour Часы
 * @param[out] nminute Минуты
 * @param[out] nsecond Секунды
 *
 * @return 0
 */
void days_to_hms(long double days, int *nhour, int *nminute, long double *nsecond);


/**
 * Определяет ∆T - временную разницу
 * между земным временем (TT) и всемирным временем (UT).
 *
 * @param[in] nyear Год, на который узнается поправка
 * @param[in] nmonth Месяц, на который узнается поправка
 *
 * @return 0
 */
long double get_deltaT(int nyear, int nmonth);