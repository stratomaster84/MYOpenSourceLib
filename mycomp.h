
#ifndef mycompH
#define mycompH
//===========================================================================
#define _USE_MATH_DEFINES
#include <math.h>
//===========================================================================
class Complex                           // класс комплексного числа
{
public:
    double  real;                   // действительная часть
    double  imaginary;              // мнимая часть

    Complex();                                 // конструктор по-умолчанию
    Complex(const Complex &z);                 // конструктор копирования
    Complex(const double &x);                  // конструктор real=x, imaginary=0.0
    Complex(const double &x,const double &y);  // конструктор real=x, imaginary=y

           Complex  operator-();               // оператор перемены знака
                    operator double();         // оператор приведения типа к double
                                               // - даёт только действительную часть

// УНАРНЫЕ ОПЕРАТОРЫ (операторы действия с присваиванием
           Complex  operator+=(const Complex &z);
           Complex  operator+=(const double  &x);

           Complex  operator-=(const Complex &z);
           Complex  operator-=(const double  &x);

           Complex  operator*=(const Complex &z);
           Complex  operator*=(const double  &x);

           Complex  operator/=(const Complex &z);
           Complex  operator/=(const double  &x);
};
//===========================================================================
const Complex Im(0.0,1.0);  // мнимая единица = sqrt(-1)
const Complex _Zero(0.0);   // комплексный нуль

//===========================================================================
//============================БИНАРНЫЕ ОПЕРАТОРЫ=============================
bool     operator==(const Complex &z1,const Complex &z2);      // операторы сравнения
bool     operator==(const Complex &z ,const double  &x );
bool     operator==(const double  &x ,const Complex &z );

bool     operator!=(const Complex &z1,const Complex &z2);
bool     operator!=(const Complex &z ,const double  &x );
bool     operator!=(const double  &x ,const Complex &z );

Complex  operator+(const Complex &z1,const Complex &z2);       // операторы сложения
Complex  operator+(const Complex &z ,const double  &x );
Complex  operator+(const double  &x ,const Complex &z );

Complex  operator-(const Complex &z1,const Complex &z2);       // операторы вычитания
Complex  operator-(const Complex &z ,const double  &x );
Complex  operator-(const double  &x ,const Complex &z );

Complex  operator*(const Complex &z1,const Complex &z2);       // операторы умножения
Complex  operator*(const Complex &z ,const double  &x );
Complex  operator*(const double  &x ,const Complex &z );

Complex  operator/(const Complex &z1,const Complex &z2);       // операторы деления
Complex  operator/(const Complex &z ,const double  &x );
Complex  operator/(const double  &x ,const Complex &z );
//===========================================================================
#endif
