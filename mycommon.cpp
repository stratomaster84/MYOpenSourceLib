//===========================================================================
#include <iomanip>
#include <fstream>
#include <windows.h>
#include "mycommon.h"

//===========================================================================
//===========ФУНКЦИИ СОХРАНЕНИЯ/ЗАГРУЗКИ ЛИСТА В ПОТОК И ОБРАТНО=============
//===========================================================================
bool SaveListToStream(const std::vector<std::string>& List, std::ostream& stream){
    if(List.empty() || !stream)
        return true;
    bool _first = true;
    for(auto str:List){
        if(_first){
            _first = false;
        }
        else{
            stream << '\n';
        }
        stream << str;
    }
    stream.flush();
    return false;
}
//===========================================================================
bool LoadListFromStream(std::istream& stream, std::vector<std::string>& List){
    if(!stream)
        return true;
    std::string str;
    while(std::getline(stream,str)){
        List.push_back(str);
    }
    return false;
}

//===========================================================================
//===============ФУНКЦИЯ СОХРАНЕНИЯ ПУЧКА В ТЕКСТОВЫЙ ФАЙЛ===================
//===========================================================================
bool PrintPuchok(const Puchok &p, const std::string& filename, int param, int mode, int Precision){
    if(!p._inited) // проверка на ошибки
            return true;

    std::ofstream Output(filename);
    if(!Output.is_open()){
        //std::ostringstream _ex("");
        //_ex << "ERROR: Can't create/rewrite file \"" << filename << "\"!" << '\n' << "May be the directory is not exists!" << std::endl;
        //throw _ex.str();
        return true;
    }

    Output << std::scientific;
    int _temp_prec = 1;
    if(Precision > 9)
        _temp_prec = 9;
    else if(Precision > 1)
        _temp_prec = Precision;
    Output << std::setprecision(_temp_prec);

    int i,j,_t1;                    // переменные индексации
    bool _first = true;             // первый элемент или нет - выводить табуляцию или нет

    if(mode == PUCHOK_3D){          // если выбрана запись ВСЕХ узлов
        _t1 = p._N - p._rang;   // начинаем с верхнего левого узла
        for(j=p._rang;j--;_t1 -= (p._rang<<1)){ // и для всех строк СВЕРХУ ВНИЗ...
            if(!_first)
                Output << '\n';
            if(param == PUCHOK_INT){                 // если выбрана "ИНТЕНСИВНОСТЬ"
                _first = true;
                for(i=0;i<p._rang;i++,_t1++){    // записать (по заданному формату) значения интенсивности во всех узлах j-ой строки через знак табуляции
                    if(_first){
                        _first = false;
                    }
                    else{
                        Output << '\t';
                    }
                    Output << C_SquareAbs(p._E[_t1]);
                }
            }//end of if(param == PUCHOK_INT)

            else{ // param == PUCHOK_ARG             // если выбрана "ФАЗА"
                _first = true;
                for(i=0;i<p._rang;i++,_t1++){    // записать (по заданному формату) значения интенсивности во всех узлах j-ой строки через знак табуляции
                    if(_first){
                        _first = false;
                    }
                    else{
                        Output << '\t';
                    }
                    Output << C_Arg(p._E[_t1]);
                }
            }//end of else(param == PUCHOK_ARG) 
        }//end of for(j=p._rang;j--;_t1 -= (p._rang<<1))
    }//end of if(mode == PUCHOK_3D)

    else{   //if(mode != PUCHOK_3D)
        if(mode == PUCHOK_VERT){        // если выбрана запись узлов на главной вертикали (сверху вниз)
            i = p._N - (p._rang>>1) - 1;    // выбрана самая верхняя точка по главной вертикали
            _t1 = -p._rang;                 // чтобы перейти на следующую - нужно опуститься строкой ниже
        }
        else if(mode == PUCHOK_HORIZ){  // если выбрана запись узлов на главной горизонтали (слева направо)
            i = p._rang*(p._rang>>1);       // выбрана самая левая точка по главной горизонтали
            _t1 = 1;                        // чтобы перейти на следующую - нужно передвинуться на следующий столбец
        }
        else if(mode == PUCHOK_DIAG13){ // если выбрана запись узлов по диагонали, проходящей в 1-ом и 3-ем квадранте (слева направо)
            i = 0;                          // выбрана левая нижняя точка
            _t1 = p._rang+1;                // чтобы перейти на следующую - нужно передвинуться на следующий столбец и подняться строкой выше
        }
        else{ //else (mode == PUCHOK_DIAG24)// если выбрана запись узлов по диагонали, проходящей в 2-ом и 4-ом квадранте (слева направо)
            i = p._N - p._rang;             // выбрана левая верхняя точка
            _t1 = -p._rang+1;               // чтобы перейти на следующую - нужно передвинуться на следующий столбец и опуститься строкой ниже
        }

        if(param == PUCHOK_INT){                // если выбрана "ИНТЕНСИВНОСТЬ"
            for(j=p._rang;j--;i+=_t1){
                if(_first){
                    _first = false;
                }
                else{
                    Output << '\t';
                }
                Output << C_SquareAbs(p._E[i]);
            }// end of for(j=p._rang;j--;i+=_t1)
        } // end of if(param == PUCHOK_INT)
        else{ // else if(param == PUCHOK_ARG)    // если выбрана "ФАЗА"
            for(j=p._rang;j--;i+=_t1){
                if(_first){
                    _first = false;
                }
                else{
                    Output << '\t';
                }
                Output << C_Arg(p._E[i]);
            }// end of for(j=p._rang;j--;i+=_t1)
        }// end of else if(param == PUCHOK_ARG)
    }//end of else(mode == PUCHOK_3D)

    Output.close();// СОХРАНЕНИЕ ЛИСТА в файл
    return false;
}

//===========================================================================
//=========== ФУНКЦИИ ПРОВЕРКИ НА СУЩЕСТВОВАНИЕ ФАЙЛА/ПАПКИ =================
//===========================================================================
int FileOrFolder(const std::string& path)
{
    DWORD ftyp = GetFileAttributesA(path.c_str());
    if (ftyp == INVALID_FILE_ATTRIBUTES)
        return 0;  //something is wrong with path!

    if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
        return 2;   // this is a directory!

    return 1;       // this is file!
}
//===========================================================================
//===========================================================================
//===========================================================================
