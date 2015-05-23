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
double utc_to_mjd(int nyear, int nmonth, int nday,
                  int nhour, int nminute, double nsecond);


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
void mjd_to_utc(double mjd, int *nyear, int *nmonth, int *nday,
                int *nhour, int *nminute, double *nsecond);


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
double utc_to_tt(int nyear, int nmonth, int nday,
                 int nhour, int nminute, double nsecond);


/**
 * Переводит Земное время в Барицентрическое динамическое
 *
 * @param[in] tt Земное время
 *
 * @return Барицентрическое динамическое время
 */
double tt_to_tdb(double tt);


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
void days_to_hms(double days, int *nhour, int *nminute, double *nsecond);


/**
 * Определяет ∆T - временную разницу
 * между земным временем (TT) и всемирным временем (UT).
 *
 * @param[in] nyear Год, на который узнается поправка
 * @param[in] nmonth Месяц, на который узнается поправка
 *
 * @return 0
 */
double get_deltaT(int nyear, int nmonth);
