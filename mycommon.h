
#ifndef mycommonH
#define mycommonH
//============================================================================
#include <iostream>
#include <vector>
#include "mypuchok.h"

// переменные для "param" в "PrintPuchok"
#define PUCHOK_INT 0
#define PUCHOK_ARG 1
// переменные для "mode" в "PrintPuchok"
#define PUCHOK_3D 0
#define PUCHOK_VERT 1
#define PUCHOK_HORIZ 2
#define PUCHOK_DIAG13 3
#define PUCHOK_DIAG24 4

//===========ФУНКЦИИ СОХРАНЕНИЯ/ЗАГРУЗКИ ЛИСТА В ПОТОК И ОБРАТНО=============
bool SaveListToStream(const std::vector<std::string>& List,std::ostream& stream);
	// false - удачно
bool LoadListFromStream(std::istream& stream,std::vector<std::string>& List);
	// false - удачно
//============================================================================

//================ФУНКЦИЯ СОХРАНЕНИЯ ПУЧКА В ТЕКСТОВЫЙ ФАЙЛ==================
bool PrintPuchok(const Puchok &p, const std::string &filename, int param, int mode, int Precision);
                        // Результирующее число начинается с минусом, если оно отрицательное и одна цифра всегда предшествует децимальной точке.
                        // Количество цифр после запятой задаётся параметром "Precision". Должно быть между 1 и 9 включительно, иначе будет округлено до 1 или 9
                        // "param" - параметр который нужно сохранить
                        //              ="PUCHOK_INT" для интенсивности
                        //              ="PUCHOK_ARG"  для фазы
                        // "mode" - определяет какие узлы сохранять
                        //              ="PUCHOK_3D" для всех узлов
                        //              ="PUCHOK_VERT" для всех узлов по главной вертикали (сверху вниз)
                        //              ="PUCHOK_HORIZ" для всех узлов по главной горизонтали (слева направо)
                        //              ="PUCHOK_DIAG13" для всех узлов по диагонали, проходящей в 1-ом и 3-ем квадранте (слева направо)
                        //
                        // ВОЗВРАЩАЕТ: true наличии ошибки, false - если ошибок нет
//============================================================================

//==ФУНКЦИИ ПРОВЕРКИ НА СУЩЕСТВОВАНИЕ ФАЙЛА/ПАПКИ==
int FileOrFolder(const std::string& path);
// 0 - неверный путь
// 1 - файл
// 2 - папка
//=================================================

#endif
