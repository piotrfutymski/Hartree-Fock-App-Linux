#include "App.h"

void App::run(const std::string & inputFile, const std::string & outputFile)
{
    _parser.parse(inputFile);

    _logger.open(outputFile);

    Mol mol{_parser.getNucleons()};

    if(_parser.getMOCount() != -1)
        mol.setMOcount(_parser.getMOCount());

    if(_parser.getApproxFlag())
        mol.setApproxIntegrals(true);

    mol.HFComputation(_logger);
    _logger.say("Drawing plots and finishing (for big molecules up to few minutes)\n");
    _outputCreator.createOutputFile(_logger, mol, _parser);

    if(_parser.getMOGraphFlag())
        _outputCreator.createMOPlainFiles(mol, _parser);
    if(_parser.getDGraphFlag())
        _outputCreator.createDPlainFile(mol, _parser);

    _outputCreator.createScriptOutputFile(mol, _parser);

    _logger.close();
    
}


void App::singlePointComputation()
{

}

void App::opt2Computation()
{
    
}