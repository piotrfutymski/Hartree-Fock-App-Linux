#include "App.h"

void App::run(const std::string & inputFile, const std::string & outputFile)
{
    _parser.parse(inputFile);

    _logger.open(outputFile);

    Mol mol = std::move(this->singlePointComputation());

    _outputCreator.createOutputFile(_logger, mol, _parser);
    _logger.close();


    std::cout<<"Creating pictures, it can take few minutes...\n";
    if(_parser.getMOGraphFlag())
        _outputCreator.createMOPlainFiles(mol, _parser);
    if(_parser.getDGraphFlag())
        _outputCreator.createDPlainFile(mol, _parser);

    _outputCreator.createScriptOutputFile(mol, _parser);
    
}


Mol App::singlePointComputation()
{
    Mol molecule{_parser.getNucleons()};

    if(_parser.getMOCount() != -1)
        molecule.setMOcount(_parser.getMOCount());

    if(_parser.getApproxFlag())
        molecule.setApproxIntegrals(true);

    molecule.HFComputation(_logger);

    return molecule;
}

void App::opt2Computation()
{
    
}