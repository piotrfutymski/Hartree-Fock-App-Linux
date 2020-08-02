#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <string>

class OutputOpenException : public std::exception{
    const char * what () const throw () {
      return "Unable to create output file";
   }
};



class Logger
{
public:
    Logger(){};
    ~Logger(){};

    void open(const std::string & filename);
    void close();

    void log(const std::string & msg, int lvl = 1);
    void say(const std::string & msg, int lvl = 1);
    void sayAndLog(const std::string & msg, int lvl = 1);

private:

    // LOG LVL: 0 - log nothing, 1-log base info, 2 -log all info

    int _logLVL {2};

    std::ofstream _outputFile;

};