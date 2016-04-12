//
// Created by Леша on 12.06.15.
//

#include <constants.h>

long double FM = 398600.4415L;   ///< Геоцентрическая гравитационная постоянная [км^3 / с^2]
long double R0 = 6378.1363L;     ///< Экваториальный радиус Земли [км]
long double A0 = 1.0 / 298.257L; ///< Полярное сжатие Земли

long double FM_M = 4.902799e03L;   ///< Гравитационная постоянная на массу Луны [км^3 / с^2]
long double FM_S = 1.32712438e11L; ///< Гравитационная постоянная на массу Солнца [км^3 / с^2]
long double P0 = 4.5606e-06L;      ///< Давление солнечного света на среднем расстоянии Земли от Солнца [н / м^2]

long double AS = 0.0L;    ///< Площадь поперечного сечения [м]
long double MS = 0.0L;    ///< Масса [кг]
long double SB = 0.0L;    ///< Баллистический коэффициент
long double KR = 0.0L;    ///< Эмпирический коэффициент отражения
long double CREFL = (10e-3 * 0.0000045606 * 0.0 * (0.0 / 1e7)); ///< Коэффициент эффективного отражения [км/с^2]

int USE_EARTH_FORCE = 1; ///< Признак задействования ускорений, обусловленных Геопотенциалом Земли
int USE_MOON_FORCE = 1;  ///< Признак задействования ускорений, обусловленных действием Луны
int USE_SUN_FORCE = 1;   ///< Признак задействования ускорений, обусловленных действием Солнца
int USE_SOLAR_FORCE = 1; ///< Признак задействования ускорений, обусловленных действием Солнечного давления
int USE_ATMOSPHERE_FORCE = 1; ///< Признак задействования ускорений, обусловленных торможением атмосферы


static void set_SB(void);
static void set_CREFL(void);

void set_default(void)
{
    FM = 398600.4415L;
    R0 = 6378.1363L;
    A0 = 1.0 / 298.257L;

    FM_M = 4.902799e03L;
    FM_S = 1.32712438e11L;
    P0 = 4.5606e-06L;

    AS = 0.0L;
    MS = 0.0L;
    SB = 0.0L;
    KR = 0.0L;
    CREFL = (10e-3 * 0.0000045606 * 0.0 * (0.0 / 1e7));

    USE_EARTH_FORCE = 1;
    USE_MOON_FORCE = 1;
    USE_SUN_FORCE = 1;
    USE_SOLAR_FORCE = 1;
    USE_ATMOSPHERE_FORCE = 1;
}

void set_FM(long double fm)
{
    FM = fm;
    return;
}
void set_R0(long double r0)
{
    R0 = r0;
    return;
}
void set_A0(long double a0)
{
    A0 = a0;
    return;
}


void set_FM_M(long double fm_m)
{
    FM_M = fm_m;
    return;
}
void set_FM_S(long double fm_s)
{
    FM_S = fm_s;
    return;
}

void set_P0(long double p0)
{
    P0 = p0;
    set_CREFL();
    return;
}


void set_SB(void)
{
    SB = (MS > 0) ? AS/MS : 0.0L;
    return;
}

void set_CREFL(void)
{
    CREFL = (MS > 0) ? (10e-3 * P0 * KR * (AS/MS)) : 0.0L;
    return;
}

void set_AS(long double as)
{
    AS = as;
    set_SB();
    set_CREFL();
    return;
}

void set_MS(long double ms)
{
    MS = ms;
    set_CREFL();
    return;
}

void set_KR(long double kr)
{
    KR = kr;
    set_CREFL();
    return;
}


void set_USE_EARTH_FORCE(int use)
{
    USE_EARTH_FORCE = (use > 0) ? 1 : 0;
    return;
}

void set_USE_MOON_FORCE(int use)
{
    USE_MOON_FORCE = (use > 0) ? 1 : 0;
    return;
}

void set_USE_SUN_FORCE(int use)
{
    USE_SUN_FORCE = (use > 0) ? 1 : 0;
    return;
}

void set_USE_SOLAR_FORCE(int use)
{
    USE_SOLAR_FORCE = (use > 0) ? 1 : 0;
    return;
}

void set_USE_ATMOSPHERE_FORCE(int use)
{
    USE_ATMOSPHERE_FORCE = (use > 0) ? 1 : 0;
    return;
}
