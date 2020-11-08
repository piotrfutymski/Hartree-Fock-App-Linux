#include "App.h"

void App::run(const std::string & inputFile, const std::string & outputFile)
{
    _parser.parse(inputFile);
    if(_parser.getCompareFlag())
        compareBasisComputation(outputFile, _parser);
    else
        singlePointComputation(outputFile, _parser);
}


void App::singlePointComputation(const std::string & outputFile, const InputParser & parser)
{
     _logger.open(outputFile);

    Mol mol{_parser.getNucleons()};

    if(_parser.getMOCount() != -1)
        mol.setMOcount(_parser.getMOCount());

    mol.HFComputation(_logger, _parser.getOptFlag());
    if(_parser.getDGraphFlag() || _parser.getMOGraphFlag())
        _logger.say("Drawing plots and finishing (for big molecules up to few minutes)\n");
    _outputCreator.createOutputFile(_logger, mol, _parser);
    _logger.close();
    if(_parser.getMOGraphFlag())
        _outputCreator.createMOPlainFiles(mol, _parser);
    if(_parser.getDGraphFlag())
        _outputCreator.createDPlainFile(mol, _parser);

    _outputCreator.createScriptOutputFile(mol, _parser);
}

void App::compareBasisComputation(const std::string & outputFile, const InputParser & parser)
{
     _logger.open(outputFile);

    Mol molA{_parser.getNucleons()};
    Mol molB{_parser.getNucleons()};

    if(_parser.getMOCount() != -1)
    {
        molA.setMOcount(_parser.getMOCount());
        molB.setMOcount(_parser.getMOCount());
    }       

    molA.HFComputation(_logger, true);
    molB.HFComputation(_logger, false);

    _logger.say("Drawing compare plot and finishing (for big molecules up to few minutes)\n");
    _outputCreator.createOutputFileCompare(_logger, molA, molB, _parser);
    _logger.close();

    _outputCreator.createMOPlainFilesCompare(molA, molB, _parser);
    _outputCreator.createDPlainFileCompare(molA, molB, _parser);

    _outputCreator.createScriptOutputFileCompare(molA, molB, _parser);
}