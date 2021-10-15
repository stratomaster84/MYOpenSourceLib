//===========================================================================
#include "mypuchok.h"
//===========================================================================
//=================КОНСТРУКТОРЫ и ДЕСТРУКТОРЫ================================
//===========================================================================
Puchok::Puchok(){                   // конструктор по-умолчанию - просто создаёт НЕИНИЦИАЛИЗИРОВАННЫЙ объект
    _inited = false;
}
//---------------------------------------------------------------------------
Puchok::Puchok(const Puchok &p){    // конструктор копирования
    _inited = false;                // создаёт НЕИНИЦИАЛИЗИРОВАННЫЙ объект
    Copy(p);                        // копирует в него "p"
}
//---------------------------------------------------------------------------
Puchok::~Puchok(){                  // деструктор - очищает и уничтожает
        Clear();
}
//===========================================================================
//=================ФУНКЦИИ ИНИЦИАЛИЗАЦИИ и КОПИРОВАНИЯ=======================
//===========================================================================
void Puchok::Clear(){
    if(_inited){
        delete [] _E;           // если инициализирован - удалить все узлы поля и разинициализировать
        _inited = false;
    }
}
//---------------------------------------------------------------------------
bool Puchok::Init(double D, double F, double Lambda, int rang, double X, int Z, double sinX, int dZ){
    if(rang<3 || !(rang % 2) || D <= 0.0 || Lambda <= 0.0 || sinX <= -1.0 || sinX >= 1.0 || dZ<-1 || dZ>1 || dZ == 0)
                                // проверка "rang" >=3 и нечётное и остальные проверки
        return true;
    _N=rang*rang;
    if(!_inited)                // если не был инициализирован
        _E = new Complex [_N];
    else if(_rang!=rang){       // если был инициализирован, но при этом количество узлов было другим
        delete [] _E;
        _E = new Complex [_N];
    }
// копируем остальные переменные (кроме "_E" и "_N")
    _rang = rang;
    _D = D;
    _F = F;
    _Lambda = Lambda;
    _H = _D/(double)(_rang+1);  // ВСЕГДА НУЛЕВЫЕ края тоже считаются и входят в "_D"

    _centerX = X;
    _centerZ = Z;
    _sinX = sinX;
    _dZ = dZ;

    _inited = true;
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::Copy(const Puchok &p){
    if(!p._inited)      // проверка инициализирован ли "p" и не ИСПОЛЬЗУЕМ ли мы "самого себя"
        return true;
    if(this==(&p))
        return false;   // если используем самого себя - ничего не делаем

    _N=p._N;
    if(!_inited)                    // если не был инициализирован
        _E = new Complex [_N];
    else if(_rang!=p._rang){        // если был инициализирован, но при этом количество узлов было другим
        delete [] _E;
        _E = new Complex [_N];
    }
    // копируем остальные переменные (кроме "_N")
    _rang = p._rang;
    _D = p._D;
    _F = p._F;
    _Lambda = p._Lambda;
    _H = p._H;

    _centerX = p._centerX;
    _centerZ = p._centerZ;
    _sinX = p._sinX;
    _dZ = p._dZ;

    for(int i=_N;i--;)
        _E[i] = p._E[i];        // копируем каждый элемент

    _inited = true;
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::Transform(const Puchok &p, double NewD, int NewRang){
    if(!p._inited || NewD <= (double)PUCHOK_MIN_DELTAS_FOR_TRANSFORM*p._H || NewRang<3 || !(NewRang % 2))
        return true;
    if(p._rang == NewRang && compare(p._D,NewD)==0){
        return Copy(p); // если новые параметры совпадают со старыми - просто копировать
    }

    int nx,ny,Cn,CPN;       // переменные для номеров точек, адреса центральной точки по столбцу(строке) и по всему пучку для НОВЫХ координат
    int oldnx,oldny,oldCn,oldCPN;   // тоже самое для СТАРЫХ координат
    int _e1,_e2,_e3,_e4,_a1,_a2,_a3,_a4;    // переменные для адресации
    double deltaX,deltaY,deltaX2,deltaY2;   // переменные для хранения коэффициентов интерполяции
    double _dst_Amp1,_dst_Amp2,_dst_Amp3,_dst_Amp4; // переменные для временного хранения целевых амплитуды и фазы
    double _dst_Fi1,_dst_Fi2,_dst_Fi3,_dst_Fi4;
    double _tmp_koef;                       // произведение всех коэффициентов интерполяции для конкретного узла
    bool Y_int,X_int;                       // для определения необходимости интерполяции по координатам

    double *_Amp,*_Fi;                      // здесь будут храниться значения амплитуды и фазы исходного пучка "p"
    _Amp = new double [p._N];
    _Fi = new double [p._N];
    for(nx = p._N;nx--;){
        _Amp[nx] = C_SquareAbs(p._E[nx]);
        _Fi[nx] = C_Arg(p._E[nx]);
    }

    _N = NewRang*NewRang;           // пишем количество точек в НОВОМ пучке
    if(!_inited)                    // если не был инициализирован
        _E = new Complex [_N];
    else if(_rang!=NewRang){        // если был инициализирован, но при этом количество узлов было другим
        delete [] _E;
        _E = new Complex [_N];
    }

    Cn = NewRang>>1;        // так как "_rang" - нечётное по определению, то "Cn" = адресу центральной координаты
    oldCn = p._rang>>1;
    CPN = Cn*(NewRang+1);   // пишем адрес центральной точки по всему пучку
    oldCPN = oldCn*(p._rang+1);

    double _newH = NewD/(double)(NewRang+1);// расстояние между узлами в НОВОМ пучке
    double _koef = _newH/p._H;              // отношение этих расстояний для НОВОГО и СТАРОГО пучков

    _E[CPN] = C_DefinePolar(sqrt(_Amp[oldCPN]),_Fi[oldCPN]);// записываем центральную координату без изменеий (интерполяция не нужна)
    int forNull = myceil(0.5*p._D/_newH);           // при увеличении размеров, с какого номера узла заполнять нулями
    if(forNull > Cn)                                // если этот номер узла больше номера центральной точки
        forNull = Cn+1;                             // то заполнять ничего не надо и пишем сюда предел для основного цикла
    else for(ny=forNull;ny<=Cn;ny++){               // иначе обнуляем
        CPN -= Cn;
        for(nx=NewRang;nx--;){                      // обнуляем всё что...
            _E[CPN+ny*NewRang+nx] = _Zero;          // ...выше "forNULL" и...
            _E[CPN-ny*NewRang+nx] = _Zero;          // ...ниже "-forNULL"
        }
        CPN += Cn;
        for(nx=forNull;nx--;){
            _E[CPN+nx*NewRang+ny] = _Zero;          // всё что ниже "forNULL" и правее "forNULL"
            _E[CPN+nx*NewRang-ny] = _Zero;          // всё что ниже "forNULL" и левее "-forNULL"
            _E[CPN-nx*NewRang+ny] = _Zero;          // всё что выше "-forNULL" и правее "forNULL"
            _E[CPN-nx*NewRang-ny] = _Zero;          // всё что выше "-forNULL" и левее "-forNULL"
        }
    }// else for(ny)
// ОСНОВНОЙ ЦИКЛ
    for(ny=1;ny<forNull;ny++){          // для всех НОВЫХ координат от "1" до "forNULL" по "y"
        deltaY = _koef*(double)ny;      // дробная СТАРАЯ координата
        oldny = myfloor(deltaY);        // целая СТАРАЯ координата

        deltaY -= (double)oldny;        // ПЕРВЫЙ коэффициент интерполяции по "Y"
        if(oldny == oldCn || _IS_ZERO(deltaY))  // определение необходимости интерполяции по "Y"
            Y_int = false;
        else
            Y_int = true;
        deltaY2 = 1.0 - deltaY;         // ВТОРОЙ коэффициент интерполяции по "Y"

        // подготовка адресов для интерполяции "КРЕСТА" (центр уже записан ВЫШЕ)
        _e1 = CPN + ny*NewRang;         // верхней ветви
        _e2 = CPN - ny*NewRang;         // нижней ветви
        _e3 = CPN + ny;                 // правой ветви
        _e4 = CPN - ny;                 // левой ветви
        _a1 = oldCPN + oldny*p._rang;   // то же для СТАРЫХ координат
        _a2 = oldCPN - oldny*p._rang;
        _a3 = oldCPN + oldny;
        _a4 = oldCPN - oldny;

        // пишем значение ОСНОВНОГО узла
            _dst_Amp1 = _Amp[_a1]*deltaY2;
            _dst_Fi1  = _Fi [_a1]*deltaY2;
            _dst_Amp2 = _Amp[_a2]*deltaY2;
            _dst_Fi2  = _Fi [_a2]*deltaY2;
            _dst_Amp3 = _Amp[_a3]*deltaY2;
            _dst_Fi3  = _Fi [_a3]*deltaY2;
            _dst_Amp4 = _Amp[_a4]*deltaY2;
            _dst_Fi4  = _Fi [_a4]*deltaY2;
        // если НЕ попали между краем и последней значащей точкой, то в интерполяции участвует и следующая (более удалённая от центра) точка
        if(Y_int){
            _dst_Amp1 += _Amp[_a1+p._rang]*deltaY;  // следующая точка по верхней ветви
            _dst_Fi1  += _Fi [_a1+p._rang]*deltaY;
            _dst_Amp2 += _Amp[_a2-p._rang]*deltaY;  // следующая точка по нижней ветви
            _dst_Fi2  += _Fi [_a2-p._rang]*deltaY;
            _dst_Amp3 += _Amp[_a3+1]      *deltaY;  // следующая точка по правой ветви
            _dst_Fi3  += _Fi [_a3+1]      *deltaY;
            _dst_Amp4 += _Amp[_a4-1]      *deltaY;  // следующая точка по левой ветви
            _dst_Fi4  += _Fi [_a4-1]      *deltaY;
        }//if(oldny!=oldCn)
        _E[_e1] = C_DefinePolar(sqrt(_dst_Amp1),_dst_Fi1);      // запись параметров в поле "Complex *_E"
        _E[_e2] = C_DefinePolar(sqrt(_dst_Amp2),_dst_Fi2);
        _E[_e3] = C_DefinePolar(sqrt(_dst_Amp3),_dst_Fi3);
        _E[_e4] = C_DefinePolar(sqrt(_dst_Amp4),_dst_Fi4);

        _e3 = _e1;_e4 = _e2;    // начальная подготовка адресов для интерполяции по квадрантам ( "КРЕСТ" и центр уже записаны ВЫШЕ)
        _a3 = _a1;_a4 = _a2;    // тоже для СТАРЫХ координат
        for(nx=1;nx<forNull;nx++){          // для всех НОВЫХ координат от "1" до "forNULL" по "X"
            deltaX = _koef*(double)nx;      // дробная СТАРАЯ координата
            oldnx = myfloor(deltaX);        // целая СТАРАЯ координата

            deltaX -= (double)oldnx;        // ПЕРВЫЙ коэффициент интерполяции по "X"
            if(oldnx == oldCn || _IS_ZERO(deltaX))  // определение необходимости интерполяции по "X"
                X_int = false;
            else
                X_int = true;
            deltaX2 = 1.0 - deltaX;         // ВТОРОЙ коэффициент интерполяции по "X"

            _e1++;_e2++;_e3--;_e4--;        // окончательная подготовка адресов для интерполяции по квадрантам
            _a1 += oldnx;_a2 += oldnx;_a3 -= oldnx;_a4 -= oldnx;    // СТАРЫЕ координаты будут не всегда последовательными
                                                                    // поэтому в конце тела цикла "for(nx)" они будут возвращены обратно
            // пишем значение ОСНОВНОГО узла
                _tmp_koef = deltaY2*deltaX2;            // запись произведения всех коэффициентов интерполяции
                _dst_Amp1 = _Amp[_a1]*_tmp_koef;
                _dst_Fi1  = _Fi [_a1]*_tmp_koef;
                _dst_Amp2 = _Amp[_a2]*_tmp_koef;
                _dst_Fi2  = _Fi [_a2]*_tmp_koef;
                _dst_Amp3 = _Amp[_a3]*_tmp_koef;
                _dst_Fi3  = _Fi [_a3]*_tmp_koef;
                _dst_Amp4 = _Amp[_a4]*_tmp_koef;
                _dst_Fi4  = _Fi [_a4]*_tmp_koef;
            if(X_int){              // если НЕ попали между краем и последней значащей точкой по "X"
                                        // в интерполяции участвует и следующая точка вдоль оси "X" (более удалённая от центра)
                _tmp_koef = deltaY2*deltaX;
                _dst_Amp1 += _Amp[_a1+1]*_tmp_koef;     // следующая точка по оси "X" в 1-ом квадранте
                _dst_Fi1  += _Fi [_a1+1]*_tmp_koef;
                _dst_Amp2 += _Amp[_a2+1]*_tmp_koef;     // следующая точка по оси "X" в 4-ом квадранте
                _dst_Fi2  += _Fi [_a2+1]*_tmp_koef;
                _dst_Amp3 += _Amp[_a3-1]*_tmp_koef;     // следующая точка по оси "X" в 2-ом квадранте
                _dst_Fi3  += _Fi [_a3-1]*_tmp_koef;
                _dst_Amp4 += _Amp[_a4-1]*_tmp_koef;     // следующая точка по оси "X" в 3-ем квадранте
                _dst_Fi4  += _Fi [_a4-1]*_tmp_koef;
            }//if(X_int)
            if(Y_int){              // если НЕ попали между краем и последней значащей точкой по "Y"
                                    // в интерполяции участвует и следующая точка вдоль оси "Y" (более удалённая от центра)
                _tmp_koef = deltaX2*deltaY;
                _dst_Amp1 += _Amp[_a1+p._rang]*_tmp_koef;       // следующая точка по оси "Y" в 1-ом квадранте
                _dst_Fi1  += _Fi [_a1+p._rang]*_tmp_koef;
                _dst_Amp2 += _Amp[_a2-p._rang]*_tmp_koef;       // следующая точка по оси "Y" в 4-ом квадранте
                _dst_Fi2  += _Fi [_a2-p._rang]*_tmp_koef;
                _dst_Amp3 += _Amp[_a3+p._rang]*_tmp_koef;       // следующая точка по оси "Y" в 2-ом квадранте
                _dst_Fi3  += _Fi [_a3+p._rang]*_tmp_koef;
                _dst_Amp4 += _Amp[_a4-p._rang]*_tmp_koef;       // следующая точка по оси "Y" в 3-ем квадранте
                _dst_Fi4  += _Fi [_a4-p._rang]*_tmp_koef;
            }//if(Y_int)
            if(X_int && Y_int){     // если НЕ попали между краем и последней значащей точкой по "X" и по "Y"
                                    // в интерполяции участвует и следующая точка вдоль осей "X" и "Y" (более удалённая от центра)
                _tmp_koef = deltaX*deltaY;
                _dst_Amp1 += _Amp[_a1+p._rang+1]*_tmp_koef;     // следующая точка по оси "X" и "Y" в 1-ом квадранте
                _dst_Fi1  += _Fi [_a1+p._rang+1]*_tmp_koef;
                _dst_Amp2 += _Amp[_a2-p._rang+1]*_tmp_koef;     // следующая точка по оси "X" и "Y" в 4-ом квадранте
                _dst_Fi2  += _Fi [_a2-p._rang+1]*_tmp_koef;
                _dst_Amp3 += _Amp[_a3+p._rang-1]*_tmp_koef;     // следующая точка по оси "X" и "Y" в 2-ом квадранте
                _dst_Fi3  += _Fi [_a3+p._rang-1]*_tmp_koef;
                _dst_Amp4 += _Amp[_a4-p._rang-1]*_tmp_koef;     // следующая точка по оси "X" и "Y" в 3-ем квадранте
                _dst_Fi4  += _Fi [_a4-p._rang-1]*_tmp_koef;
            }//if(X_int && Y_int)

            _E[_e1] = C_DefinePolar(sqrt(_dst_Amp1),_dst_Fi1);  // запись параметров в поле "Complex *_E"
            _E[_e2] = C_DefinePolar(sqrt(_dst_Amp2),_dst_Fi2);
            _E[_e3] = C_DefinePolar(sqrt(_dst_Amp3),_dst_Fi3);
            _E[_e4] = C_DefinePolar(sqrt(_dst_Amp4),_dst_Fi4);

            _a1 -= oldnx;_a2 -= oldnx;_a3 += oldnx;_a4 += oldnx; // как и обещали в конце тела цикла "for(nx)" возвращаем обратно адреса СТАРЫХ точек
        }// for(nx)
    }// for(ny)
    delete [] _Amp;
    delete [] _Fi;

// запись/копирование остальных параметров (кроме "_E" и "_N")
    _rang = NewRang;
    _D = NewD;
    _F = p._F;
    _Lambda = p._Lambda;
    _H = _D/(double)(_rang+1);

    _centerX = p._centerX;
    _centerZ = p._centerZ;
    _sinX = p._sinX;
    _dZ = p._dZ;

    _inited = true;
    return false;
}
//===========================================================================
//=================ФУНКЦИИ ПОЛУЧЕНИЯ ХАРАКТЕРИСТИК ПУЧКА=====================
//===========================================================================
double Puchok::FullPower() const{
    if(!_inited)
        return 0.0;                     // если не инициализирован - возвращаем мощность = 0
    double _pow = 0.0;
    for(int i = _N;i--;)
        _pow += C_SquareAbs(_E[i]);     // вычисляем сумму интенсивностей в каждом узле
    _pow *= _H*_H;                      // и умножаем их на площадь квадрата ограниченного четырьмя соседними узлами
    return _pow;
}
//---------------------------------------------------------------------------
double Puchok::DiamOfCircleWithPower(double Power) const{
    if(!_inited || Power<=0.0)
        return 0.0;
    if(Power > FullPower())
        return 0.0;

    double _M;                         // сюда запишем результат GetAngleWithW
    double _W = Power/_H/_H;           // вычисляем требуемый поток
    _M = GetAngleWithW(_E,_W,_rang);   // вычисляем угол (в узлах), в котором содержится данный поток
    if(_M == 0.0)                      // возвратить "0" мог только в случае, если требуемого потока нету в пучке
        return HUGE_VAL;               // поэтому возвращаем "Бесконечность"

    return 2.0*_M*_H;                  // диаметр равен двойному углу умноженному на расстояние между узлами
}
//---------------------------------------------------------------------------
double Puchok::DivergenceAngle(double Ppart, double FULLANGLE, int rang) const{
    if(!_inited || Ppart<=0.0 || FULLANGLE <= 0.0 || FULLANGLE >= M_PI || rang<3 || !(rang % 2) )
        return 0.0;         // если не инициализирован или "Ppart" <=0 - возвращаем расходимость = 0
    if(Ppart>=1.0)
        return M_PI;        // если по уровню полной мощности - угол равен "ПИ"
    double _result;
    Puchok _far;
    _far.WriteFarField(*this,FULLANGLE,rang);    //расчёт дальней зоны
    _result = _far.DiamOfCircleWithPower(Ppart*FullPower());
    _far.Clear();
    return _result;
}
//===========================================================================
//=================ФУНКЦИИ РАЗЛИЧНЫХ ПРЕОБРАЗОВАНИЙ ПУЧКА====================
//===========================================================================
bool Puchok::CircleDiafragm(double Diam, bool Hole, double _delta){
    if(!_inited)
        return true;
    double _R = 0.5*Diam;
    if(_R >= M_SQRT1_2*_D  || _R <= _H || _delta < 0.0 || _delta > _R)
        return true;

    double _Mmax = (_R+_delta)/_H;
    double _Mmin = (_R-_delta)/_H;
    double _rad,_par,_XX;

    int nx,ny,Cn = _rang>>1;// так как "_rang" - нечётное по определению, то "Cn" = центральной координате
    int CPN = Cn*(_rang+1); // пишем адрес центральной точки по всему пучку

    if(Hole){                       // для вырезки центральной части (дырки)
        _E[CPN] = _Zero;            // стираем осевой узел - обе координаты равны нулю
        for(ny=1;ny<=Cn;ny++){
            _rad = (double)ny;
            if(_rad < _Mmin){                       // стираем "КРЕСТ" - одна из координат равна нулю
                _E[CPN+ny*_rang] = _Zero;           // ниже "Mmin"
                _E[CPN-ny*_rang] = _Zero;           // выше "-Mmin"
                _E[CPN+ny]       = _Zero;           // левее "-Mmin"
                _E[CPN-ny]       = _Zero;           // правее "Mmin"
            }//if(_rad < _Mmin)
            else if(_rad < _Mmax){
                _par = (_rad-_Mmin)/(_Mmax-_Mmin);
                _XX = 0.5*(1.0 + sin(M_PI*(_par-0.5)));
                _E[CPN+ny*_rang] *= _XX;            // ниже "Mmax", но выше "Mmin"
                _E[CPN-ny*_rang] *= _XX;            // выше "-Mmax", но ниже "-Mmin"
                _E[CPN+ny]       *= _XX;            // левее "-Mmax", но правее "-Mmin"
                _E[CPN-ny]       *= _XX;            // правее "Mmax", но левее "Mmin"
            }//else if(_rad < _Mmax)
            for(nx=1;nx<=Cn;nx++){
                _rad = sqrt(double(ny*ny + nx*nx));
                if(_rad < _Mmin){                   // стираем всё остальное
                    _E[CPN+ny*_rang+nx] = _Zero;    // первый квадрант
                    _E[CPN-ny*_rang+nx] = _Zero;    // четвёртый квадрант
                    _E[CPN+ny*_rang-nx] = _Zero;    // второй квадрант
                    _E[CPN-ny*_rang-nx] = _Zero;    // третий квадрант
                }//if(_rad < _Mmin)
                else if(_rad < _Mmax){
                    _par = (_rad-_Mmin)/(_Mmax-_Mmin);
                    _XX = 0.5*(1.0 + sin(M_PI*(_par-0.5)));
                    _E[CPN+ny*_rang+nx] *= _XX;     // первый квадрант
                    _E[CPN-ny*_rang+nx] *= _XX;     // четвёртый квадрант
                    _E[CPN+ny*_rang-nx] *= _XX;     // второй квадрант
                    _E[CPN-ny*_rang-nx] *= _XX;     // третий квадрант
                }//else if(_rad < _Mmax)
            }//for(nx)
        }//for(ny)
    }//if(Hole)
    else{
        for(ny=1;ny<=Cn;ny++){
            _rad = (double)ny;
            if(_rad > _Mmax){                   // стираем "КРЕСТ" - одна из координат равна нулю
                _E[CPN+ny*_rang] = _Zero;       // выше "Mmax"
                _E[CPN-ny*_rang] = _Zero;       // ниже "-Mmax"
                _E[CPN+ny]       = _Zero;       // правее "Mmax"
                _E[CPN-ny]       = _Zero;       // левее "-Mmax"
            }//if(_rad > _Mmax)
            else if(_rad > _Mmin){
                _par = (_rad-_Mmin)/(_Mmax-_Mmin);
                _XX = 0.5*(1.0 - sin(M_PI*(_par-0.5)));
                _E[CPN+ny*_rang] *= _XX;        // выше "Mmin", но ниже "Mmax"
                _E[CPN-ny*_rang] *= _XX;        // ниже "-Mmin", но выше "-Mmax"
                _E[CPN+ny]       *= _XX;        // правее "Mmin", но левее "Mmax"
                _E[CPN-ny]       *= _XX;        // левее "-Mmin", но правее "-Mmax"
            }//else if(_rad > _Mmin)
            for(nx=1;nx<=Cn;nx++){
                _rad = sqrt(double(ny*ny + nx*nx));
                if(_rad > _Mmax){                // стираем всё остальное
                    _E[CPN+ny*_rang+nx] = _Zero; // первый квадрант
                    _E[CPN-ny*_rang+nx] = _Zero; // четвёртый квадрант
                    _E[CPN+ny*_rang-nx] = _Zero; // второй квадрант
                    _E[CPN-ny*_rang-nx] = _Zero; // третий квадрант
                }//if(_rad > _Mmax)
                else if(_rad > _Mmin){
                    _par = (_rad-_Mmin)/(_Mmax-_Mmin);
                    _XX = 0.5*(1.0 - sin(M_PI*(_par-0.5)));               // стираем всё остальное
                    _E[CPN+ny*_rang+nx] *= _XX;   // первый квадрант
                    _E[CPN-ny*_rang+nx] *= _XX;   // четвёртый квадрант
                    _E[CPN+ny*_rang-nx] *= _XX;   // второй квадрант
                    _E[CPN-ny*_rang-nx] *= _XX;   // третий квадрант
                }//else if(_rad > _Mmin)
            }//for(nx)
        }//for(ny)
    } // else(Hole)
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::YDiafragm(double Y,bool Under, double _delta){
    if(!_inited || _delta < 0.0)
        return true;
    double _Mmax = (Y + _delta)/_H + (double)(_rang>>1);
    double _Mmin = (Y - _delta)/_H + (double)(_rang>>1);
    double _par,_XX;

    int nx,ny;
    int mmin = myfloor(_Mmin);
    int mmax = myfloor(_Mmax);

    if(Under){
        ny = mmax;
        if(ny>=_rang)
            ny = _rang;
        while(ny>=0 && ny > mmin){
            _par = ((double)ny-_Mmin)/(_Mmax-_Mmin);
            _XX = 0.5*(1.0 + sin(M_PI*(_par-0.5)));
            for(nx = _rang;nx--;)
                _E[ny*_rang + nx] *= _XX;
            ny--;
        }
        while(ny>=0){
            for(nx = _rang;nx--;)
                _E[ny*_rang + nx] = _Zero;
            ny--;
        }
    }
    else{
        ny = _rang-1;
        while(ny>=0 && ny > mmax){
            for(nx = _rang;nx--;)
                _E[ny*_rang + nx] = _Zero;
            ny--;
        }
        while(ny>=0 && ny > mmin){
            _par = ((double)ny-_Mmin)/(_Mmax-_Mmin);
            _XX = 0.5*(1.0 - sin(M_PI*(_par-0.5)));
            for(nx = _rang;nx--;)
                _E[ny*_rang + nx] *= _XX;
            ny--;
        }
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::XDiafragm(double X,bool Left, double _delta){
    if(!_inited || _delta < 0.0)
        return true;
    double _Mmax = (X + _delta)/_H + (double)(_rang>>1);
    double _Mmin = (X - _delta)/_H + (double)(_rang>>1);
    double _par,_XX;

    int nx,ny;
    int mmin = myfloor(_Mmin);
    int mmax = myfloor(_Mmax);

    if(Left){
        nx = mmax;
        if(nx>=_rang)
            nx = _rang;
        while(nx>=0 && nx > mmin){
            _par = ((double)nx-_Mmin)/(_Mmax-_Mmin);
            _XX = 0.5*(1.0 + sin(M_PI*(_par-0.5)));
            for(ny = _rang;ny--;)
                _E[ny*_rang + nx] *= _XX;
            nx--;
        }
        while(nx>=0){
            for(ny = _rang;ny--;)
                _E[ny*_rang + nx] = _Zero;
            nx--;
        }
    }
    else{
        nx = _rang-1;
        while(nx>=0 && nx > mmax){
            for(ny = _rang;ny--;)
                _E[ny*_rang + nx] = _Zero;
            nx--;
        }
        while(nx>=0 && nx > mmin){
            _par = ((double)nx-_Mmin)/(_Mmax-_Mmin);
            _XX = 0.5*(1.0 - sin(M_PI*(_par-0.5)));
            for(ny = _rang;ny--;)
                _E[ny*_rang + nx] *= _XX;
            nx--;
        }
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::Sub(const Puchok &p){
    if(!_inited || !p._inited || _rang != p._rang)
        return true;
                
    for(int i=_N;i--;){
        _E[i] -= p._E[i];               // из каждого элемента вычитается соответствующий элемент пучка "p"
#if _ZERO_CHECK
        _ZERO_FUNC(_E[i],_Zero);        // проверка на НУЛЬ
#endif
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::InverseX(){
    if(!_inited)
        return true;
    int _t1,i;
    Complex _tmp;
    for(_t1 = _N + (_rang>>1);(_t1-=_rang)>0;){
        for(i=(_rang>>1)+1;--i;){
            _tmp = _E[_t1 + i];             // перестановка
            _E[_t1+i] = _E[_t1-i];
            _E[_t1-i] = _tmp;
        }
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::InverseY(){
    if(!_inited)
        return true;
    int _t1,i,j;
    Complex _tmp;

    _t1 = (_rang>>1)*_rang;
    for(j=(_rang>>1)+1;--j;){
        for(i=_rang;i--;){
            _tmp = _E[_t1+j*_rang+i];       // перестановка
            _E[_t1+j*_rang+i] = _E[_t1-j*_rang+i];
            _E[_t1-j*_rang+i] = _tmp;
        }
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::TransformFToPhase(double n){
    if(!_inited || _F == 0.0)       // если плоский фронт - ничего делать не надо
        return true;

    int i,j,_t1;
    double _K,_F2,_H2,_fabsF,_amp;  // вспомогательные переменные для учёта радиуса кривизны ВФ

    _K = _2PI/_Lambda*n;            // волновое число
    _fabsF = _F;                    // модуль "_F" - пока равен "_F"
    if(_F < 0.0){                   // если фокус отрицательный - знаки волнового числа и "_fabsF" поменять
        _K = -_K;
        _fabsF = -_F;
    }
    _F2 = _F*_F;                    // "_F^2"

    _amp = sqrt(_F2 - 0.25*_D*_D)/_fabsF;   // коэффициент уменьшения интенсивности
    _D /= _amp;                     // диаметр увеличивается при проецировании на экран
    _H = _D/(double)(_rang+1);      // и расстояние между узлами тоже

    _H2 = _H*_H;                    // "_H^2"

    int Cn = (_rang>>1);            // центральная точка
    _t1 = _N;
    for(j=_rang;j--;)
        for(i=_rang;i--;)
            _E[--_t1] *= C_DefinePolar(_amp,_K*(sqrt(_F2+_H2*(double)((i-Cn)*(i-Cn)+(j-Cn)*(j-Cn))) - _fabsF));
    _F = 0.0;                       // теперь можно обнулить _F
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::DefineEAsAverage(const Complex *E, double K){
    if(!_inited || K <= 0.0 || K >1.0)            // проверка на ошибки
        return true;
    double _amp,_arg;       // переменные для хранения значений амплитуды и фазы в текущем узле
    double _K2 = 1.0 - K;   // весовой коэффициент для первоначальных значений пучка
    for(int i=_N;i--;){
        _amp = sqrt(  _K2*C_SquareAbs(_E[i]) + K*C_SquareAbs(E[i])  );  // амплитуда - среднее квадратическое с весами
        _arg = _K2*C_Arg(_E[i]) + K*C_Arg(E[i]);                        // фаза - среднее арифметическое
        _E[i] = C_DefinePolar(_amp,_arg);
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::Mul(Complex Z){
    if(!_inited)
        return true;
    for(int i=_N;i--;)
        _E[i] *= Z;
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::PhaseAddition(const double *Arg){
    if(!_inited)
        return true;
    for(int i=_N;i--;)
        _E[i] *= Complex(cos(Arg[i]),sin(Arg[i]));
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::InvertInt()
{
    if(!_inited)
        return true;
    double _koef = 0.0;
    double _temp = 0.0;
    int i;
    for(i=_N;i--;)           // получаем максимальное значение интенсивности
        if((_temp=C_SquareAbs(_E[i])) > _koef)
            _koef = _temp;
    for(i=_N;i--;){
        _temp = _koef - C_SquareAbs(_E[i]);
        if(_temp<0.0)
            _E[i] = C_DefinePolar(0.0,C_Arg(_E[i]));
        _E[i] = C_DefinePolar(sqrt(_temp),C_Arg(_E[i]));
    }
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::Shade(const Puchok &shadow)
{
    if(!_inited || !shadow._inited || _rang != shadow._rang)
        return true;
    for(int i=_N;i--;)
        _E[i] *= shadow._E[i];
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::CircleTurn(double Angle)
{
    if(!_inited)
        return true;
    if(_IS_ZERO(Angle))
        return false;

    int i,j,_t;
    double _tmp_i,_tmp_j;   // временные переменные для прибавок
    double *_Amp,*_Fi;      // здесь будут храниться значения амплитуды и фазы исходного пучка
    double iold,jold;       // старые дробные координаты
    int fli,flj;            // крайняя левая нижняя координата для интерполяции
    double i_koef,j_koef;
    double _amp,_fi;

    _tmp_i = (double)(_rang>>1)*(1.0 - cos(Angle) + sin(Angle));
    _tmp_j = (double)(_rang>>1)*(1.0 - sin(Angle) - cos(Angle));

    _Amp = new double [_N];
    _Fi = new double [_N];
    for(i = _N;i--;){
        _Amp[i] = C_SquareAbs(_E[i]);
        _Fi[i] = C_Arg(_E[i]);
    }
    _t = _N-1;
    for(j=_rang;j--;)
        for(i=_rang;i--;_t--){
            iold = (double)i*cos(Angle) - (double)j*sin(Angle) + _tmp_i;
            jold = (double)i*sin(Angle) + (double)j*cos(Angle) + _tmp_j;
            fli = myfloor(iold);
            flj = myfloor(jold);
            i_koef = iold - (double)fli;
            j_koef = jold - (double)flj;
            _amp = 0.0;
            _fi = 0.0;
            if(!(fli<0 || fli>=_rang || flj<0 || flj>=_rang)){
                _amp += (1.0 - i_koef)*(1.0 - j_koef)*_Amp[flj*_rang + fli];
                _fi  += (1.0 - i_koef)*(1.0 - j_koef)*_Fi[flj* _rang + fli];
            }
            fli++;
            if(!(fli<0 || fli>=_rang || flj<0 || flj>=_rang)){
                _amp += i_koef*(1.0 - j_koef)*_Amp[flj*_rang + fli];
                _fi  += i_koef*(1.0 - j_koef)*_Fi[flj* _rang + fli];
            }
            fli--;
            flj++;
            if(!(fli<0 || fli>=_rang || flj<0 || flj>=_rang)){
                _amp += (1.0 - i_koef)*j_koef*_Amp[flj*_rang + fli];
                _fi  += (1.0 - i_koef)*j_koef*_Fi[flj* _rang + fli];
            }
            fli++;
            if(!(fli<0 || fli>=_rang || flj<0 || flj>=_rang)){
                _amp += i_koef*j_koef*_Amp[flj*_rang + fli];
                _fi  += i_koef*j_koef*_Fi[flj* _rang + fli];
            }
            _E[_t] = C_DefinePolar(sqrt(_amp),_fi);
        }
    delete [] _Amp;
    delete [] _Fi;
    return false;
}
//---------------------------------------------------------------------------
double Puchok::Normalize()
{
    if(!_inited)
        return 0.0;

    double _koef = 0.0;
    double _temp = 0.0;
    int i;
    for(i=_N;i--;)                      // получаем максимальное значение интенсивности
        if((_temp=C_SquareAbs(_E[i])) > _koef)
            _koef = _temp;
    if(_koef != 0.0 && _koef != 1.0){
        _koef = sqrt(_koef);            // получаем максимальное значение интенсивности
        for(i=_N;i--;){                 // нормируем
            _E[i] /= _koef;
#if _ZERO_CHECK
            _ZERO_FUNC(_E[i],_Zero);    // проверка на нуль
#endif
        }
    }
    return _koef;
}
//===========================================================================
//=================ФУНКЦИИ РАСПРОСТРАНЕНИЯ===================================
//===========================================================================
bool Puchok::ShiftNearRange(const Puchok &p,const Complex *s, double Hz){
    if(!p._inited || Hz <= 0.0)
        return true;
    if(p._F < 0.0 && Hz >= (-p._F)) // проверка на переход через фокус
        return true;
    _N=p._N;
    if(!_inited)                    // если не был инициализирован
        _E = new Complex [_N];
    else if(_rang!=p._rang){        // если был инициализирован, но при этом количество узлов было другим
        delete [] _E;
        _E = new Complex [_N];
    }
// Копирование и запись...
    _rang = p._rang;
    _D = p._D;
    _F = p._F;
    _Lambda = p._Lambda;

    _sinX = p._sinX;
    _dZ = p._dZ;

    _centerX = p._centerX + _sinX*Hz;
    _centerZ = p._centerZ + _dZ;

// Вычисление...
    int i,j,_t1;

    double _Nsred=0.0;
    for(i=_N;i--;)
        _Nsred += s[i].real;
    _Nsred /= (double)_N;   // средний показатель преломления в параметрах "s"

    double _koef = 0.0;     // доп. константа для сферического фронта
    double _Hxy = p._H;     // здесь будет записано расстояние между узлами поля на среднем (половинном) сечении

    if(_F != 0.0){          // если сферический фронт
        _koef = _Hxy/_F;
        _Hxy *= 1.0 + 0.5*Hz/_F;    // расстояние между узлами на среднем (половинном) сечении увеличено по сравнению с начальным
    }
    double K_Hxy = _Hxy/_Lambda*_2PI;   // (волновое число)*(_Hxy)
    double KN_Hxy = K_Hxy*_Nsred;       // (волновое число)*(_Hxy)*(средний показатель преломления)
    double _Hxy_Hz = _Hxy/Hz;

    Complex _paramC1,_paramF1;          // параметры для прогонки по каждой строке
    Complex _paramC2,_paramF2;          // параметры для прогонке по каждому столбцу
    if(_F != 0.0){
        _paramC1 = -2.0 + Complex(0.0, 4.0*_Hxy_Hz*KN_Hxy+2.0*_koef*KN_Hxy) - KN_Hxy*KN_Hxy;
        _paramF1 = -2.0 + Complex(0.0,-4.0*_Hxy_Hz*KN_Hxy);
        _paramC2 = -2.0 + Complex(0.0, 4.0*_Hxy_Hz*KN_Hxy);
        _paramF2 = -2.0 + Complex(0.0,-4.0*_Hxy_Hz*KN_Hxy+2.0*_koef*KN_Hxy) - KN_Hxy*KN_Hxy;
    }
    else{
        _paramC1 = -2.0 + Complex(0.0, 4.0*_Hxy_Hz*KN_Hxy) - KN_Hxy*KN_Hxy;
        _paramF1 = -2.0 + Complex(0.0,-4.0*_Hxy_Hz*KN_Hxy);
        _paramC2 = -2.0 + Complex(0.0, 4.0*_Hxy_Hz*KN_Hxy);
        _paramF2 = -2.0 + Complex(0.0,-4.0*_Hxy_Hz*KN_Hxy) - KN_Hxy*KN_Hxy;
    }

    Complex *ParamSred, *SredSloy;
    SredSloy = new Complex[_N];     // сюда будут записываться результаты прогонки по каждой строке
    ParamSred = new Complex[_N];    // сюда будут записываться коэффициенты в волновом уравнении стоящие перед мгновенным значением функции
    Complex _tmpC,_tmpF;
    Complex *_M,*_L;
    _M = new Complex[_rang+1];      // вспомогательные массивы для записи коэффициентов прямого хода прогонки
    _L = new Complex[_rang+1];
    _M[_rang] = _Zero;
    _L[_rang] = _Zero;
    _t1=_N;

    for(j=_rang;j--;){              // для каждой строки вычисляем сначала прямой ход прогонки, потом обратный
        for(i=_rang;i--;){          // прямой ход (но в обратном направлении от "_rang-1" до "0")
            _t1--;

            double _tmp = K_Hxy*s[_t1].real;    //(волновое число)*(_Hxy)*(показатель преломления в узле)
            ParamSred[_t1] = Complex(_tmp*_tmp,-_tmp*s[_t1].imaginary*_Hxy);
                    // _tmp^2 - Im*(волновое число)*(показ. преломления)*(коэф. усиления)*(Hxy^2)
                    // здесь НЕ УЧИТЫВАЕТСЯ "-0.25*(Hxy^2)*(коэф. усиления)^2" ввиду своей малости относительно (_tmp^2)
                    // для его учёта можно РАЗкоментить следующие две строки и ЗАкоментить предыдущие две

             //Complex _tmp = Complex(K_Hxy*s[_t1].real,-0.5*_Hxy*s[_t1].imaginary);
             //ParamSred[_t1] = _tmp*_tmp;

            _tmpC = _paramC1 + ParamSred[_t1];      // коэффициент при значении функции в данном узле
            _tmpF = _paramF1*p._E[_t1];             // коэффициент в правой части уравнения (перенесённый в левую)
            if(j)
                _tmpF += p._E[_t1-_rang];           // прибавляем к нему значения в верхнем...
            if(j!=_rang-1)
                _tmpF += p._E[_t1+_rang];           // ...и нижнем узлах
            _L[i] = -1.0/(_tmpC + _L[i+1]);         // считаем элементы массивов прямого хода прогонки
            _M[i] = _L[i]*(_M[i+1]+_tmpF);
        }//for(i)
        SredSloy[_t1] = _M[0];          // обратный ход (но в прямом направлении от "0" до "_rang-1")
        for(i=1;i<_rang;i++){
            _t1++;
            SredSloy[_t1] = _M[i] + _L[i]*SredSloy[_t1-1];
        }//for(i)
                _t1-= _rang-1;
    }//for(j)
    for(i=_rang;i--;){                  // для каждого столбца вычисляем сначала прямой ход прогонки, потом обратный
        _t1 = _N + i;
        for(j=_rang;j--;){              // прямой ход (но в обратном направлении от "_rang-1" до "0")
            _t1-=_rang;

            _tmpC = _paramC2;           // коэффициент при значении функции в данном узле
            _tmpF = (_paramF2 + ParamSred[_t1])*SredSloy[_t1];
                                        // коэффициент в правой части уравнения  (перенесённый в левую)
            if(i)
                _tmpF += SredSloy[_t1-1];   // прибавляем к нему значения в левом...
            if(i!=_rang-1)
                _tmpF += SredSloy[_t1+1];   // ...и правом узлах
            _L[j] = -1.0/(_tmpC + _L[j+1]); // считаем элементы массивов прямого хода прогонки
            _M[j] = _L[j]*(_M[j+1]+_tmpF);
        }//for(j)
        _E[_t1] = _M[0];                    // обратный ход (но в прямом направлении от "0" до "_rang-1")
        for(j=1;j<_rang;j++){
            _t1+=_rang;
            _E[_t1] = _M[j] + _L[j]*_E[_t1-_rang];
        }//for(j)
    }//for(i)
#if _ZERO_CHECK
    for(_t1=_N;_t1--;)
        _ZERO_FUNC(_E[_t1],_Zero);          // проверка на нуль
#endif

    delete [] _M;
    delete [] _L;
    delete [] SredSloy;
    delete [] ParamSred;
// поправка скопированных переменных и запись незаписанных...

    if(_F != 0.0){              // если сферический фронт:
        _D *= 1.0 + Hz/_F;      // увеличиваем диаметр пучка
        _F += Hz;               // увеличиваем радиус кривизны пучка
    }
    _H = _D/(double)(_rang+1);  // считаем расстояние между соседними узлами
    _inited = true;
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::ShiftMidRange(const Puchok &p, Complex defparam, double Hz,int rang){
    if(!p._inited || rang<3 || !(rang % 2) || Hz <= 0.0)
        return true;
    if(p._F < 0.0 && Hz >= (-p._F))                 // проверка на переход через фокус
        return true;

                                                    // ограничение на "Hz" чтобы пучки не перекрывались
    double _Apovt = Hz*p._Lambda/defparam.real/p._H;// координата с которой распределение будет повторяться (т.е. где фазовый множитель отличается на "_2PI")
    double _Dnew;                                   // предполагаемый размер нового поля = геометрический + "_DEFRAG_N" зон Френеля
    if(p._F == 0.0)
        _Dnew = p._D +  sqrt(_DEFRAG_N*p._Lambda/defparam.real*Hz);
    else
        _Dnew = p._D*(Hz+p._F)/p._F  + sqrt( _DEFRAG_N*p._Lambda/defparam.real*Hz*p._F/(Hz+p._F));
    if(_Apovt < _Dnew)                              // если повторный пучок перекрывается с основным - ошибка
        return true;

// Копирование и запись...
    _H = p._H;
    _F = p._F;
    _Lambda = p._Lambda;

    _sinX = p._sinX;
    _dZ = p._dZ;

    _centerX = p._centerX + _sinX*Hz;
    _centerZ = p._centerZ + _dZ;

// Вычисление
    int i,j,_t1;
    double _Ro_H,_g;
    Complex _tmp;

    int _fin = rang + p._rang - 1;  // количество различных комбинаций фазового множителя в принципе Гюйгенса-Френеля

    int _period=1,T=1;              // вычислим период - _period = 2^T >= _fin
    while(_period<_fin){
        _period <<= 1;
        T++;
    }
    T--;

    Complex B(0.5*_H*defparam.imaginary,_2PI*_H/_Lambda*defparam.real);     // константа учитывающая характеристики среды
    double C = Hz/_H;                                                       // отношение шагов по "XY" и по "Z"
    double A = C*C;
    double _par1 = _2PI/_period*(p._rang-1);                                // константа для поправки (сдвига) преобразования Фурье (учёт отрицательных (x2-x1) )
    Complex _par2(0.0,-0.5*_H/_Lambda*defparam.real/_period/_period);       // общий множитель для всех "_E"
    double _par3 = 1.0, _par4 = 0.0;                                                     // дополнительные константы для неплоских ВФ
    if(_F != 0.0){
        _par3 = 1.0 + Hz/_F;
        _par4 = 0.5*_H/_F;
    }
    int _Fi = (rang>>1) - (p._rang>>1) + p._rang - 1;                       // разность адресов середин пучков (нового и старого)...
                                                                            // ...плюс учёт отрицательных (x2-x1)

    Complex **U,**P;                                                        // массивы для записи амплитуды поля и множителей
    Complex *koef;                                                          // вспомогательный массив для поправки (сдвига) преобразования Фурье

    koef = new Complex [_period];
    U = new Complex* [_period];
    P = new Complex* [_period];
    for(j=_period;j--;){
        koef[j] = C_DefinePolar(1.0,_par1*j);       // определяем поправки и обнуляем "U" и "P"
        U[j] = new Complex[_period];
        P[j] = new Complex[_period];
        for(i=_period;i--;){
            U[j][i] = _Zero;
            P[j][i] = _Zero;
        }
    }

    _t1 = p._N;
    for(j=p._rang;j--;){
        for(i=p._rang;i--;)
            U[j][i] = p._E[--_t1];
        C_FFT(U[j],T,0);
    }                                               // после этого в "U" преобразование Фурье поля "_E" по "x"
    for(j=_fin;j--;){
        if(_F == 0.0) for(i=_fin;i--;){
            _Ro_H = sqrt(A+(double)(   (i-_Fi)*(i-_Fi)  +  (j-_Fi)*(j-_Fi)  ));
            P[j][i] = C_Exp(B*_Ro_H)/_Ro_H*(1.0+C/_Ro_H);
        }
        else for(i=_fin;i--;){
            _g = _par3*(double)(   (i-_Fi)*(i-_Fi)  +  (j-_Fi)*(j-_Fi)  );
            _Ro_H = sqrt(A+_g);
            P[j][i] = C_Exp(B*_Ro_H)/_Ro_H*(1.0+(C-_par4*_g)/_Ro_H);
        }

        C_FFT(P[j],T,0);
        for(i=_period;i--;)
            P[j][i] *= koef[i];
    }                                               // после этого в "P" преобразование Фурье (сдвинутое) множителей по "X"

// U[i][j] = U[j][i].....
// P[i][j] = P[j][i].....
    for(j=_period;--j;){
        for(i=j;i--;){
            _tmp=U[j][i];
            U[j][i] = U[i][j];
            U[i][j] = _tmp;
            _tmp=P[j][i];
            P[j][i] = P[i][j];
            P[i][j] = _tmp;
        } // for(i)
    } // for(j)                                 // поменяли координаты "X" и "Y" для удобства преобразования

    for(i=_period;i--;){
        C_FFT(U[i],T,0);                        // преобразования Фурье "U" и "P" по "Y"
        C_FFT(P[i],T,0);
        for(j=_period;j--;)
            U[i][j] *= P[i][j]*koef[j];         // перемножение преобразованных "U" и "P" (с учётом сдвига "P" по "Y")

        delete [] P[i];

        C_FFT(U[i],T,1);                        // обратное преобразование "U" по "Y"
    }
    delete [] P;
    delete [] koef;

// U[j][i] = U[i][j]...
    for(j=_period;--j;){
        for(i=j;i--;){
            _tmp=U[j][i];
            U[j][i] = U[i][j];
            U[i][j] = _tmp;
        } // for(i)
    } // for(j)                                 // снова меняем координаты "X" и "Y" для удобства преобразования

    _t1 = _N = rang*rang;
    if(!_inited)                                // если не был инициализирован
        _E = new Complex [_N];
    else if(_rang!=rang){                       // если был инициализирован, но при этом количество узлов было другим
        delete [] _E;
        _E = new Complex [_N];
    }

    for(j=rang;j--;){
        C_FFT(U[j],T,1);                        // обратное преобразование "U" по "X"
        for(i=rang;i--;){
            _E[--_t1] = U[j][i]*_par2;          // запись конечного поля
#if _ZERO_CHECK
            _ZERO_FUNC(_E[_t1],_Zero);
#endif
        }
        delete [] U[j];
    }
    delete [] U;

// поправка скопированных переменных и запись незаписанных...
    if(_F != 0.0){              // если сферический фронт:
        _F += Hz;               // увеличиваем радиус кривизны пучка
        _H *= _par3;            // и расстояние между узлами
    }
    _rang = rang;
    _D = _H*(double)(_rang+1);
    _inited = true;
    return false;
}
//===========================================================================
//=================ФУНКЦИИ РАСЧЁТОВ ДАЛЬНЕЙ ЗОНЫ=============================
//===========================================================================
bool Puchok::ShiftFarRange(const Puchok &p, double Hz, int rang, double D){
    if(!p._inited || Hz <= 0.0 || D <= 0.0 || rang<3 || !(rang % 2) )       // проверка на ошибки
        return true;

    Complex *field;                     // массив для записи нового поля
    double NewH,_delta;                 // переменные для записи расстояния между узлами
    field = new Complex [rang*rang];

    NewH = D/(double)(rang+1);          // расстояние между узлами в мм
    _delta = NewH/Hz;                   // угловое расстояние

    p.GetFarField(field,_delta,rang);   // получим распределение в дальней зоне
    if(_inited)                         // если "this" инициализирован...
        delete [] _E;                   // ...удалить значения интенсивностей
    _E = field;

    _N=rang*rang;
    for(int i=_N;i--;){                 // делаем из углового потока - привычное поле
        _E[i] /= Hz;
#if _ZERO_CHECK
        _ZERO_FUNC(_E[i],_Zero);        // проверка на нуль
#endif
    }
// копирование и запись...
    _D = D;                             // не будем уменьшать диаметр на "sqrt(_F^2-D^2)/|_F|"- т.к. "D << |_F|"
    _Lambda = p._Lambda;
    _H = NewH;
    _rang = rang;
    _sinX = p._sinX;
    _dZ = p._dZ;
    _centerZ = p._centerZ + _dZ;
    _centerX = p._centerX + _sinX*Hz;
    _F = Hz;                            // радиус кривизны ВФ теперь равен "Hz"

    _inited = true;
    return false;
}
//---------------------------------------------------------------------------
bool Puchok::WriteFarField(const Puchok &p, double Angle, int rang){
    if(!p._inited  || Angle <= 0.0 || Angle >= M_PI || rang<3 || !(rang % 2) )
        return true;

    double _delta = Angle/(double)(rang-1); // используем "rang-1", а не "rang+1" как везде, чтобы были значащие края
    Complex *field;
    field = new Complex [(_N=rang*rang)];
    p.GetFarField(field,_delta,rang);       // получаем распределение в дальней зоне

    if(_inited)                             // если "this" инициализирован...
    delete [] _E;                           // ...удалить значения интенсивностей
    _E = field;

    _inited = true;

// копирование и запись...
    _rang = rang;
    _D = _delta*(_rang+1);
    _Lambda = p._Lambda;
    _H = _delta;
    _F = 0.0;

    _centerX = 0.0;
    _centerZ = 0;
    _sinX = 0.0;
    _dZ = 1;

    return false;
}
//===========================================================================
//=================ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ===================================
//===========================================================================
bool Puchok::GetFarField(Complex *result, double Delta, int rang) const{
    if(!_inited || rang<3 || !(rang % 2) || Delta <= 0.0)   // проверка на ошибки
        return true;

    int ind1,ind2,ind3;     // переменные для записи вспомогательных индексов
    int i,j,a,b;            // переменные для организации циклов
    int Cn1 = _rang>>1;     // индекс центрального узла в строке (столбце) в исходном пучке
    int Cn2 = rang>>1;      // индекс центрального узла в строке (столбце) в результирующем пучке

    Complex *W,*U,*NewE;    // массивы для записи фазовых множителей и промежуточных значений плотности
    Complex _tmp;
    W = new Complex [rang*_rang];
    U = new Complex [rang*_rang];

    double A;
    Complex B;

    if(_F == 0.0){                      // ЕСЛИ плоский ВФ...
        NewE = _E;                      // ...используем исходное поле
        A = -_H*Delta*_2PI/_Lambda;     // вспомогательная константа для фазового множителя
        B = Complex(0.0,-_H*_H/_Lambda);// вспомогательная константа для коррекции значения плотности
    }
    else{                               // ИНАЧЕ...
        double _K,_F2,_H2,_fabsF,_amp;  // вспомогательные переменные для учёта радиуса кривизны ВФ
        _fabsF = _F;                    // модуль "_F" - пока равен "_F"
        _K = _2PI/_Lambda;              // волновое число
        _F2 = _F*_F;                    // "_F^2"
        _amp = sqrt((_F2 - 0.25*_D*_D)/_F2);   // коэффициент уменьшения интенсивности
        _H2 = _H*_H/_amp/_amp;          // расстояние между узлами увеличилось

        A = -sqrt(_H2)*Delta*_K;
        B = Complex(0.0,-_H2/_Lambda);

        if(_F < 0.0){                   // если фокус отрицательный - поменять знаки "_K" и "_fabsF"
            _K = -_K;
            _fabsF = -_F;
        }

        NewE = new Complex [_N];        // создаём новое поле и...
        ind1 = _N-1;
        for(j=_rang;j--;)
            for(i=_rang;i--;ind1--)     // ...учитываем в нём радиус кривизны ВФ
                NewE[ind1] = _E[ind1]*C_DefinePolar(_amp,_K*(sqrt(_F2+_H2*(double)((i-Cn1)*(i-Cn1)+(j-Cn1)*(j-Cn1))) - _fabsF));
    }
    ind1 = _rang*rang;
    for(a=rang;a--;)                    // запись фазовых множителей без сферической составляющей
        for(i=_rang;i--;)
            W[--ind1] = C_DefinePolar(1.0,A*(double)(a-Cn2)*(double)(i-Cn1));

    ind1 = _rang*rang;
    for(a=rang;a--;){                   // запись промежуточных значений плотности
        ind2 = _N;                      // Uaj = SUMi (_Eij*Wai)
        for(j=_rang;j--;){
            _tmp = _Zero;
            ind3 = (a+1)*_rang;
            for(i=_rang;i--;)
                _tmp += NewE[--ind2]*W[--ind3];
            U[--ind1] = _tmp;
        }
    }

    ind1 = rang*rang;                   // запись основных значений плотности
    for(b=rang;b--;){                   // Fab = B*SUMj (_Uaj*Wbj)
        ind2 = _rang*rang;
        for(a=rang;a--;){
            _tmp = _Zero;
            ind3 = (b+1)*_rang;
            for(j=_rang;j--;)
                _tmp += U[--ind2]*W[--ind3];
            result[--ind1] = _tmp*B;    // проверки на нуль нет, так как это вспомогательная функция
        }
    }
    if(_F != 0.0)
        delete [] NewE;                 // удалять новое поле, только если было создано в теле функции
    delete [] W;                        // очистка вспомогательных массивов
    delete [] U;
    return false;
}
//---------------------------------------------------------------------------
double Puchok::GetAngleWithW(const Complex *field, double W,int rang) const{
    if(W<=0.0 || rang<3 || !(rang % 2))     // если полный заданный поток меньше "0" - значит ошибка
        return 0.0;     // поэтому возвращаем нуль
    int i,j;                // индексы
    int Cn = rang>>1;       // центральная точка по строке (столбцу)
    int CPN = Cn*(rang+1);  // пишем адрес центральной точки по всему полю

    int M,M2,MM,MmMm;       // переменные для хранения
    double sum1,sum2,tmpSUM;// записывается основные и резервная (для восстановления) копии суммы интенсивностей
    int _tmp;               // дополнительная переменная для хранения суммы квадратов координат
    int Mmax;               // переменная максимального радиуса


                    // записываем в обе СУММЫ значения ближайших пяти центральных узлов
    tmpSUM = sum1 = C_SquareAbs(field[CPN])         +
                    C_SquareAbs(field[CPN + rang])  +
                    C_SquareAbs(field[CPN - rang])  +
                    C_SquareAbs(field[CPN + 1])     +
                    C_SquareAbs(field[CPN - 1]);
    if(sum1 >= W)                               // если уже больше...
        return sqrt(W/sum1);                    // ...возвращаем радиус между "0" и "1"

// ПЕРВЫЙ ОСНОВНОЙ цикл
    MM = 1;                                     // на момент входа в ПЕРВЫЙ основной цикл MM = (M-1)^2
    Mmax = ceil(M_SQRT2*(double)Cn);
    for(M=2;M<=Mmax;M++){
        if(M<=Cn){      // если не вышли за границы пучка - прибавляем "КРЕСТ" с координатами "M"
            sum1 += C_SquareAbs(field[CPN + M*rang]) +
                    C_SquareAbs(field[CPN - M*rang]) +
                    C_SquareAbs(field[CPN + M])      +
                    C_SquareAbs(field[CPN - M]);
        }// if(M<=Cn)
        MmMm = MM;                              // MmMm = (M-1)^2
        MM += (M<<1)-1;                         // MM  = M^2
        for(i=1;i<M && i<=Cn;i++){              // для всех возможных "i" и "j" меньших "M" провести проверку
            for(j=1;j<M && j<=Cn;j++){
                _tmp = i*i + j*j;               // сумма квадратов двух координат
                if(_tmp <= MM && _tmp > MmMm){  // если в диапазоне от "(M-1)^2" до "M^2" прибавляем этот узел со всех квадрантов
                    sum1 += C_SquareAbs(field[CPN + j*rang + i]) +
                            C_SquareAbs(field[CPN - j*rang + i]) +
                            C_SquareAbs(field[CPN + j*rang - i]) +
                            C_SquareAbs(field[CPN - j*rang - i]);
                }// if(_tmp)
            }// for(j)
        }// for(i)
        if(sum1>=W)         // если получившийся поток больше заданного...
            break;          // ...завершить первый основной цикл
        tmpSUM = sum1;      // сохранить текущее значение
    }// for(M)
    if(M>Mmax)      // если вышли из цикла не через "break" - значит заданного потока не содержится в максимально возможном круге
        return 0.0; // поэтому возвращаем нуль
                    // если же вышли через "break" - значит поток находится в круге угловым радиусом между "(M-1)" и "M"
    sum2 = tmpSUM;  // восстанавливаем сумму и количество точек соответсвующее угловому радиусу "(M-1)"
//ВТОРОЙ основной цикл
    MM--;           // MM = M^2 - 1
    for(M2=MmMm;M2<MM;M2++){// перебираем все возможные варианты угловых радиусов от "sqrt((M-1)^2)" до "sqrt(M^2-1)"
        for(i=1;i<M && i<=Cn;i++){              // для всех возможных "i" и "j" меньших "M" провести проверку
            for(j=1;j<M && j<=Cn;j++){
                _tmp = i*i + j*j;               // сумма квадратов двух координат
                if(_tmp <= M2+1 && _tmp>M2){    // если в диапазоне от "M2" до "M2+1" прибавляем этот узел со всех квадрантов
                    sum2 += C_SquareAbs(field[CPN + j*rang + i]) +
                            C_SquareAbs(field[CPN - j*rang + i]) +
                            C_SquareAbs(field[CPN + j*rang - i]) +
                            C_SquareAbs(field[CPN - j*rang - i]);
                }// if(_tmp)
            }//for(j)
        }//for(i)
        if(sum2>=W)     // если получившийся поток больше заданного...
            break;      // ...завершить второй основной цикл
        tmpSUM = sum2;  // сохранить текущее значение
    }// for(M2)
    if(M2 == MM)        // если вышли из цикла не через "break" - значит заданный поток содержится в круге с угловым радиусом между "sqrt(M^2-1)" и "M"
        return sqrt((double)MM + (W-sum2)/(sum1-sum2));
                        // если же вышли через "break" - значит поток находится в круге с угловым радиусом между "sqrt(M2)" и "sqrt(M2+1)"
    return sqrt((double)M2 + (W-tmpSUM)/(sum2-tmpSUM));
}
//===========================================================================
//===========================================================================
//===========================================================================

