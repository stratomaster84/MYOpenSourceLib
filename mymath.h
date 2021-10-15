
#ifndef mymathH
#define mymathH
//---------------------------------------------------------------------------
#include "mycomp.h"
#define _2PI 6.283185307179586476926            // число 2*Пи
#define _ZERO_EPSILON 1e-20                     // минимальное квадратичное отклонение от нуля при проверке на близость к НУЛЮ
#define _ZERO_FUNC(Z,_Value) if(C_SquareAbs(Z) < _ZERO_EPSILON) Z = _Value  // функция проверки на близость к НУЛЮ
#define _IS_ZERO(Z) ((C_SquareAbs(Z) < _ZERO_EPSILON)?true:false)         // функция - является ли "Z" - нулём
//---------------------------------------------------------------------------
double  C_SquareAbs(const Complex &z);                 // вычисляет квадрат модуля "z"
double  C_Arg(const Complex &z, int n = 0);            // вычисляет фазу числа "z" (n-ый "корень")
Complex C_Conj(const Complex &z);                      // вычисляет комплексно-сопряжённое число

Complex C_DefinePolar(const double &wAbs,const double &wArg); // возвращает число wAbs*exp(Im*wArg)
Complex C_Exp(const Complex &z);                       // вычисляет exp(z)
Complex C_Ln(const Complex &z, int n = 0);             // вычисляет натуральный логорифм "z" (n-ый "корень")
Complex C_Pow(const Complex &z1, const Complex &z2,
                                 int n = 0);           // вычисляет z1^z2 (n-ый "корень")

Complex C_Sin(const Complex &z);                       // вычисляет синус "z"
Complex C_Cos(const Complex &z);                       // вычисляет косинус "z"
Complex C_Tg(const Complex &z);                        // вычисляет тангенс "z"
Complex C_Sh(const Complex &z) ;                       // вычисляет гиперболический синус "z"
Complex C_Ch(const Complex &z);                        // вычисляет гиперболический косинус "z"
Complex C_Th(const Complex &z);                        // вычисляет гиперболический тангенс "z"

Complex C_Ash(const Complex &z, int n = 0);            // вычисляет гиперболический ареасинус "z" (n-ый "корень")
Complex C_Ach(const Complex &z, int n = 0);            // вычисляет гиперболический ареакосинус "z" (n-ый "корень")
Complex C_Ath(const Complex &z, int n = 0);            // вычисляет гиперболический ареатангенс "z" (n-ый "корень")
Complex C_Asin(const Complex &z, int n = 0);           // вычисляет арксинус "z" (n-ый "корень")
Complex C_Acos(const Complex &z, int n = 0);           // вычисляет арккосинус "z" (n-ый "корень")
Complex C_Atg(const Complex &z, int n = 0);            // вычисляет арктангенс "z" (n-ый "корень")

void    C_FFT(Complex *z, int T, int reverse);
                                                // вычисляет БПФ последовательности "z" размерности 2^T,
                                                /* reverse: 0 - прямое преобразование Фурье,
                                                            1 - обратное преобразование Фурье БЕЗ деления всех элементов на 2^T
                                                            2 - обратное преобразование Фурье С делением всех элементов на 2^T */
int     myfloor(const double &x);                      // аналог "floor"
int     myceil(const double &x);                       // аналог "ceil"
int     compare(const double &x, const double &y);            // правильно сравнивает "double"
                                                // возвращает "1" - если "x>y", "-1" - если "x<y" и "0" - если "x==y"
#endif
