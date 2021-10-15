//===========================================================================
#include "mycomp.h"
//===========================================================================
Complex::Complex()
{
}
//===========================================================================
Complex::Complex(const Complex &z)
{
    real = z.real;
    imaginary = z.imaginary;
}
//===========================================================================
Complex::Complex(const double &x)
{
    real = x;
    imaginary = 0.0;
}
//===========================================================================
Complex::Complex(const double &x,const double &y)
{
    real = x;
    imaginary = y;
}
//===========================================================================
Complex Complex::operator-()
{
    return Complex(-real,-imaginary);
}
//===========================================================================
Complex::operator double(){
        return real;
}
//===========================================================================
Complex Complex::operator+=(const Complex &z)
{
    real += z.real;
    imaginary += z.imaginary;
    return *this;
}
//--------------------------------------------
Complex Complex::operator+=(const double &x)
{
    real += x;
    return *this;
}
//===========================================================================
Complex Complex::operator-=(const Complex &z)
{
    real -= z.real;
    imaginary -= z.imaginary;
    return *this;
}
//--------------------------------------------
Complex Complex::operator-=(const double &x)
{
    real -= x;
    return *this;
}
//===========================================================================
Complex Complex::operator*=(const Complex &z)
{
    const double tmp = real*z.real - imaginary*z.imaginary;
    imaginary = real*z.imaginary + imaginary*z.real;
    real = tmp;
    return *this;
}
//--------------------------------------------
Complex Complex::operator*=(const double &x)
{
    real *= x;
    imaginary *= x;
    return *this;
}
//===========================================================================
Complex Complex::operator/=(const Complex &z)
{
    const double _mod = z.real*z.real + z.imaginary*z.imaginary;
    if(_mod == 0.0)
        real = imaginary = HUGE_VAL;  // при делении на нуль - в обе части числа записывается HUGE_VAL
    else{
        const double tmp = (real*z.real + imaginary*z.imaginary)/_mod;
        imaginary = (-real*z.imaginary + imaginary*z.real)/_mod;
        real = tmp;
    }
    return *this;
}
//--------------------------------------------
Complex Complex::operator/=(const double &x)
{
    if(x == 0.0)
        real = imaginary = HUGE_VAL;  // при делении на нуль - в обе части числа записывается HUGE_VAL
    else{
        real /= x;
        imaginary /= x;
    }
    return *this;
}
//===========================================================================
bool operator==(const Complex &z1,const Complex &z2)
{
    if(z1.real != z2.real || z1.imaginary != z2.imaginary)
        return false;
    else
        return true;
}
//--------------------------------------------
bool operator==(const Complex &z ,const double  &x )
{
    if(z.real != x || z.imaginary != 0.0)
        return false;
    else
        return true;
}
//--------------------------------------------
bool operator==(const double  &x ,const Complex &z )
{
    if(z.real != x || z.imaginary != 0.0)
        return false;
    else
        return true;
}
//===========================================================================
bool operator!=(const Complex &z1,const Complex &z2)
{
    if(z1.real != z2.real || z1.imaginary != z2.imaginary)
        return true;
    else
        return false;
}
//--------------------------------------------
bool operator!=(const Complex &z ,const double  &x )
{
    if(z.real != x || z.imaginary != 0.0)
        return true;
    else
        return false;
}
//--------------------------------------------
bool operator!=(const double  &x ,const Complex &z )
{
    if(z.real != x || z.imaginary != 0.0)
        return true;
    else
        return false;
}
//===========================================================================
Complex operator+(const Complex &z1,const Complex &z2)
{
    return Complex(z1.real + z2.real,
                   z1.imaginary + z2.imaginary);
}
//--------------------------------------------
Complex operator+(const Complex &z ,const double  &x )
{
    return Complex(z.real+x,
                   z.imaginary);
}
//--------------------------------------------
Complex operator+(const double  &x ,const Complex &z )
{
    return Complex(z.real+x,
                   z.imaginary);
}
//===========================================================================
Complex operator-(const Complex &z1,const Complex &z2)
{
    return Complex(z1.real - z2.real,
                   z1.imaginary - z2.imaginary);
}
//--------------------------------------------
Complex operator-(const Complex &z ,const double  &x )
{
    return Complex(z.real-x,
                   z.imaginary);
}
//--------------------------------------------
Complex operator-(const double  &x ,const Complex &z )
{ 
    return Complex(x - z.real,
                   -z.imaginary);
}
//===========================================================================
Complex operator*(const Complex &z1,const Complex &z2)
{
    return Complex(z1.real*z2.real - z1.imaginary*z2.imaginary,
                   z1.real*z2.imaginary + z1.imaginary*z2.real);
}
//--------------------------------------------
Complex operator*(const Complex &z ,const double  &x )
{
    return Complex(z.real*x,
                   z.imaginary*x);
}
//--------------------------------------------
Complex operator*(const double  &x ,const Complex &z )
{
    return Complex(z.real*x,
                   z.imaginary*x);
}
//===========================================================================
Complex operator/(const Complex &z1,const Complex &z2)
{
    const double _mod = z2.real*z2.real + z2.imaginary*z2.imaginary;
    if(_mod == 0.0)
        return Complex (HUGE_VAL,HUGE_VAL);  // при делении на нуль - в обе части числа записывается HUGE_VAL
    else
        return Complex((z1.real*z2.real + z1.imaginary*z2.imaginary)/_mod,
                       (-z1.real*z2.imaginary + z1.imaginary*z2.real)/_mod);
}
//--------------------------------------------
Complex operator/(const Complex &z ,const double  &x )
{
    if(x == 0.0)
        return Complex (HUGE_VAL,HUGE_VAL);  // при делении на нуль - в обе части числа записывается HUGE_VAL
    else
        return Complex(z.real/x,
                       z.imaginary/x);
}
//--------------------------------------------
Complex operator/(const double  &x ,const Complex &z )
{
    const double _mod = z.real*z.real + z.imaginary*z.imaginary;
    if(_mod == 0.0)
        return Complex (HUGE_VAL,HUGE_VAL);  // при делении на нуль - в обе части числа записывается HUGE_VAL
    else
        return Complex(x*z.real/_mod,
                       -x*z.imaginary/_mod);
}
//===========================================================================
