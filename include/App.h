#pragma once
#include <string>
#include "InputParser.h"
#include "Logger.h"
#include "OutputCreator.h"
#include "Mol.h"

class App
{
public:
    App(){};
    ~App(){};

    void run(const std::string & inputFile, const std::string & outputFile);

private:

    InputParser _parser;

    Logger _logger;

    OutputCreator _outputCreator;

private:

    Mol singlePointComputation();

    void opt2Computation();

};