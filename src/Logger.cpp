#include "Logger.h"

void Logger::open(const std::string & filename)
{
    _outputFile.open(filename, std::ios::trunc | std::ios::out);

    if(!_outputFile.is_open())
        throw OutputOpenException();


}
void Logger::close()
{
    _outputFile.close();
}

void Logger::log(const std::string & msg, int lvl)
{
    _outputFile<<msg;
}
void Logger::say(const std::string & msg, int lvl)
{
    std::cout<<msg;
}
void Logger::sayAndLog(const std::string & msg, int lvl)
{
    this->log(msg);
    this->say(msg);
}