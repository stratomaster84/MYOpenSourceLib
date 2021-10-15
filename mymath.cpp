//===========================================================================
#include "mymath.h"
//===========================================================================
double C_SquareAbs(const Complex &z)
{
    return z.real*z.real + z.imaginary*z.imaginary;
}
//--------------------------------------------
double C_Arg(const Complex &z, int n)
{
    if(z.real != 0.0 || z.imaginary != 0.0){
        if(n == 0)
            return atan2(z.imaginary,z.real);
        else
            return atan2(z.imaginary,z.real) + double(n)*_2PI;
    }
    else
        return 0.0;
}
//--------------------------------------------
Complex C_Conj(const Complex &z)
{
    return Complex(z.real,
                   -z.imaginary);
}
//===========================================================================
Complex C_DefinePolar(const double &wAbs,const double &wArg)
{
    return Complex(wAbs*cos(wArg),
                   wAbs*sin(wArg));
}
//--------------------------------------------
Complex C_Exp(const Complex &z)
{
    const double _exp = exp(z.real);
    return Complex(_exp*cos(z.imaginary),
                   _exp*sin(z.imaginary));
}
//--------------------------------------------
Complex C_Ln(const Complex &z, int n)
{
    if(z.real == 0.0 && z.imaginary == 0.0)
        return Complex(-HUGE_VAL,0.0);
    double _arg = atan2(z.imaginary,z.real);
    if(n != 0)
        _arg += double(n)*_2PI;
    return Complex(0.5*log(z.real*z.real + z.imaginary*z.imaginary),
                   _arg);
}
//--------------------------------------------
Complex C_Pow(const Complex &z1, const Complex &z2, int n)
{
    if((z1.real == 0.0 || z1.real == 1.0) && z1.imaginary == 0.0)
        return z1;
    const double _log = 0.5*log(z1.real*z1.real + z1.imaginary*z1.imaginary);
    double _atg = atan2(z1.imaginary,z1.real);
    if(n != 0)
        _atg += double(n)*_2PI;
    const double _exp = exp(z2.real*_log - _atg*z2.imaginary);
    const double _arg = z2.imaginary*_log + _atg*z2.real;
    return Complex(_exp*cos(_arg),
                   _exp*sin(_arg));
}
//===========================================================================
Complex C_Sin(const Complex &z)
{
    return Complex(sin(z.real)*cosh(z.imaginary),
                   cos(z.real)*sinh(z.imaginary));
}
//--------------------------------------------
Complex C_Cos(const Complex &z)
{
    return Complex(cos(z.real)*cosh(z.imaginary),
                   -sin(z.real)*sinh(z.imaginary));
}
//--------------------------------------------
Complex C_Tg(const Complex &z)
{
    const double _s = sin(z.real);
    const double _c = cos(z.real);
    const double _sh = sinh(z.imaginary);
    const double _ch = cosh(z.imaginary);
    const double _mod = _c*_c*_ch*_ch+_s*_s*_sh*_sh;
    if(_mod == 0.0)
        return Complex(HUGE_VAL,HUGE_VAL);
    else
        return Complex(_s*_c/_mod,
                       _ch*_sh/_mod);
}
//--------------------------------------------
Complex C_Sh(const Complex &z)
{
    return Complex(cos(z.imaginary)*sinh(z.real),
                   sin(z.imaginary)*cosh(z.real));
}
//--------------------------------------------
Complex C_Ch(const Complex &z)
{
    return Complex(cos(z.imaginary)*cosh(z.real),
                   sin(z.imaginary)*sinh(z.real));
}
//--------------------------------------------
Complex C_Th(const Complex &z)
{
    const double _s = sin(z.imaginary);
    const double _c = cos(z.imaginary);
    const double _sh = sinh(z.real);
    const double _ch = cosh(z.real);
    const double _mod = _c*_c*_ch*_ch+_s*_s*_sh*_sh;
    if(_mod == 0.0)
        return Complex(HUGE_VAL,HUGE_VAL);
    else
        return Complex(_sh*_ch/_mod,
                       _c*_s/_mod);
}
//===========================================================================
Complex C_Ash(const Complex &z, int n)
{
    return C_Ln(z+C_Pow(z*z+1.0 , 0.5), n);
}
//--------------------------------------------
Complex C_Ach(const Complex &z, int n)
{
    return C_Ln(z+C_Pow(z*z-1.0 , 0.5), n);
}
//--------------------------------------------
Complex C_Ath(const Complex &z, int n)
{
    return Complex(0.5)*C_Ln((1.0+z)/(1.0-z), n);
}
//--------------------------------------------
Complex C_Asin(const Complex &z, int n)
{
    Complex _iz(-z.imaginary,z.real);
    return Complex(0.0,-1.0)*C_Ln(_iz+C_Pow(_iz*_iz+1.0 , 0.5), n);
}
//--------------------------------------------
Complex C_Acos(const Complex &z, int n)
{
    Complex _iz(-z.imaginary,z.real);
    return (M_PI_2 + Complex(0.0,1.0)*C_Ln(_iz+C_Pow(_iz*_iz+1.0 , 0.5), n));
}
//--------------------------------------------
Complex C_Atg(const Complex &z, int n)
{
    Complex _iz(-z.imaginary,z.real);
    return Complex(0.0,0.5)*C_Ln((1.0-_iz)/(1.0+_iz), n);
}
//===========================================================================
//===============Вспомогательные массивы для C_FFT===========================
//===========================================================================
// Этот массив содержит числа от 0 до 255 с обратным порядком битов
unsigned char reverse256[256]= {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0,
    0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8,
    0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4,
    0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC,
    0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2,
    0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA,
    0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6,
    0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE,
    0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1,
    0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9,
    0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5,
    0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED,
    0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3,
    0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB,
    0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7,
    0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF,
    0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};
//--------------------------------------------
//Этот массив содержит числа C_Exp(-Im*_2PI/2^n) для n= 1,...,32
Complex W2n[32]={
    Complex(-1.0                              ,  0.0                               ), // n = 1
    Complex(0.0                               , -1.0                               ), // n = 2
    Complex(0.70710678118654752440084436210485, -0.70710678118654752440084436210485), // n = 3
    Complex(0.92387953251128675612818318939679, -0.38268343236508977172845998403040), // n = 4
    Complex(0.98078528040323044912618223613424, -0.19509032201612826784828486847702), // n = 5
    Complex(0.99518472667219688624483695310948, -9.80171403295606019941955638886e-2), // n = 6
    Complex(0.99879545620517239271477160475910, -4.90676743274180142549549769426e-2), // n = 7
    Complex(0.99969881869620422011576564966617, -2.45412285229122880317345294592e-2), // n = 8
    Complex(0.99992470183914454092164649119638, -1.22715382857199260794082619510e-2), // n = 9
    Complex(0.99998117528260114265699043772857, -6.13588464915447535964023459037e-3), // n = 10
    Complex(0.99999529380957617151158012570012, -3.06795676296597627014536549091e-3), // n = 11
    Complex(0.99999882345170190992902571017153, -1.53398018628476561230369715026e-3), // n = 12
    Complex(0.99999970586288221916022821773877, -7.66990318742704526938568357948e-4), // n = 13
    Complex(0.99999992646571785114473148070739, -3.83495187571395589072461681181e-4), // n = 14
    Complex(0.99999998161642929380834691540291, -1.91747597310703307439909561989e-4), // n = 15
    Complex(0.99999999540410731289097193313961, -9.58737990959773458705172109764e-5), // n = 16
    Complex(0.99999999885102682756267330779455, -4.79368996030668845490039904946e-5), // n = 17
    Complex(0.99999999971275670684941397221864, -2.39684498084182187291865771650e-5), // n = 18
    Complex(0.99999999992818917670977509588385, -1.19842249050697064215215615969e-5), // n = 19
    Complex(0.99999999998204729417728262414778, -5.99211245264242784287971180889e-6), // n = 20
    Complex(0.99999999999551182354431058417300, -2.99605622633466075045481280835e-6), // n = 21
    Complex(0.99999999999887795588607701655175, -1.49802811316901122885427884615e-6), // n = 22
    Complex(0.99999999999971948897151921479472, -7.49014056584715721130498566730e-7), // n = 23
    Complex(0.99999999999992987224287980123973, -3.74507028292384123903169179084e-7), // n = 24
    Complex(0.99999999999998246806071995015625, -1.87253514146195344868824576593e-7), // n = 25
    Complex(0.99999999999999561701517998752946, -9.36267570730980827990672866808e-8), // n = 26
    Complex(0.99999999999999890425379499688176, -4.68133785365490926951155181385e-8), // n = 27
    Complex(0.99999999999999972606344874922040, -2.34066892682745527595054934190e-8), // n = 28
    Complex(0.99999999999999993151586218730510, -1.17033446341372771812462135032e-8), // n = 29
    Complex(0.99999999999999998287896554682627, -5.85167231706863869080979010083e-9), // n = 30
    Complex(0.99999999999999999571974138670657, -2.92583615853431935792823046906e-9), // n = 31
    Complex(0.99999999999999999892993534667664, -1.46291807926715968052953216186e-9)  // n = 32
};
//==================================C_FFT====================================
// вычисляет БПФ последовательности "z" размерности 2^T,
/* reverse: 0 - прямое преобразование Фурье,
            1 - обратное преобразование Фурье БЕЗ деления всех элементов на 2^T
            2 - обратное преобразование Фурье С делением всех элементов на 2^T */

void C_FFT(Complex *z, int T, int reverse)
{
    unsigned int i, j, N, Nd2, k, m, mpNd2, Skew,Nmax;
    unsigned char *Ic = (unsigned char*) &i;
    unsigned char *Jc = (unsigned char*) &j;
    Complex tmp;
    Complex *Wstore, *Warray;
    Complex WN, W, Temp, *pWN;

    Nmax = 1<<T;

    //первая перестановка
    for(i = 1; i < Nmax - 1; i++){
        Jc[0] = reverse256[Ic[3]];
        Jc[1] = reverse256[Ic[2]];
        Jc[2] = reverse256[Ic[1]];
        Jc[3] = reverse256[Ic[0]];
        j >>= (32 - T);
        if (i < j){
            tmp = z[i];
            z[i] = z[j];
            z[j] = tmp;
        } //if(i<j)
    }// for(i)

    //выделяем память под поворачивающие множители
    Wstore = new Complex[ Skew = (Nmax >> 1) ];
    Wstore[0] = Complex(1.0);

    //Главный цикл
    for(N = 2, Nd2 = 1, pWN = W2n ; N <= Nmax; Nd2 = N, N += N, pWN++, Skew >>= 1){
        //WN = W(1, N) = C_Exp(-Im*_2PI/N)
        WN= *pWN;
        if (reverse)
            WN.imaginary = -WN.imaginary;
        for(Warray = Wstore, k = 0; k < Nd2; k++, Warray += Skew){
            if (k & 0x1){
                W *= WN;
                *Warray = W;
            }
            else
                W = *Warray;

            for(m = k; m < Nmax; m += N){
                mpNd2 = m + Nd2;
                Temp = W;
                Temp *= z[mpNd2];
                z[mpNd2] = z[m];
                z[mpNd2] -= Temp;
                z[m] += Temp;
            }//for(m)
        }//for(Warray)
    }//for(N = 2, Nd2 = 1, pWN = W2n)

    delete [] Wstore;

    if (reverse==2)
        for( i = 0; i < Nmax; i++ )
            z[i] /= Nmax;
}
//===========================================================================
int myfloor(const double &x){
    const double _tmp = floor(x);
    const double _d = x-_tmp-1.0;
    if(_IS_ZERO(_d))
        return (int)_tmp + 1;
    else
        return (int)_tmp;
}
//--------------------------------------------
int myceil(const double &x){
    const double _tmp = ceil(x);
    const double _d = _tmp-x-1.0;
    if(_IS_ZERO(_d))
        return (int)_tmp - 1;
    else
        return (int)_tmp;
}
//--------------------------------------------
int compare(const double &x, const double &y){
    const double _d = x-y;
    if(_IS_ZERO(_d))
        return 0;
    if(_d>0.0)
        return 1;
    else
        return -1;
}
//===========================================================================
