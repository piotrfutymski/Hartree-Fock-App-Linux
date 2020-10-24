#pragma once
#include <cmath>
#include <string>

#define PI 3.14159265358979323846
#define M_SQRTPI 1.77245385090551602729816748334
#define BOHR_TO_ANGSTROM 0.52917721
#define HARTREE_TO_ELECTRONOVOLTS 27.211396641308
#define ERROR 1.0e-6


constexpr int NUCLEON_COUNT = 26;
const std::string NUCLEON_TAB[] ={
    "X","H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","X","Ti","X","X","X","Fe"
};