#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <string>
#include "Logger.h"
#include "Mol.h"
#include "InputParser.h"

class OutputCreator
{
public:

    OutputCreator(){};
    ~OutputCreator(){}

    void createOutputFile(Logger & logger, Mol & mol, InputParser & parser);

    void createOutputFileCompare(Logger & logger, Mol & molA, Mol & B, InputParser & parser);

    void createMOPlainFiles(Mol & mol, InputParser & parser);

    void createMOPlainFilesCompare(Mol & molA, Mol & molB, InputParser & parser);

    void createDPlainFile(Mol & mol, InputParser & parser);

    void createDPlainFileCompare(Mol & molA,Mol & molB, InputParser & parser);

    void createScriptOutputFile(Mol & mol, InputParser & parser);

    void createScriptOutputFileCompare(Mol & molA, Mol & molB, InputParser & parser);

private:

    int _xPoints, _yPoints;    
    double _xMin, _yMin;


};