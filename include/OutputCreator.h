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

    void createMOPlainFiles(Mol & mol, InputParser & parser);

    void createDPlainFile(Mol & mol, InputParser & parser);

    void createScriptOutputFile(Mol & mol, InputParser & parser);

private:

    int _xPoints, _yPoints;    
    double _xMin, _yMin;


};