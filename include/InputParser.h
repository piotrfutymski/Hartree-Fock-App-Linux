#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include "Nucleon.h"

class InputOpenException : public std::exception{
    const char * what () const throw () {
      return "Unable to open input file";
   }
};

class InputFileException : public std::exception{
    const char * what () const throw () {
      return "Unable to read input file";
   }
};


class InputParser
{
public:

    enum class Plane{
    XY, XZ, YZ
    };

    InputParser(){}
    ~InputParser(){}

    void parse(const std::string & input);

    std::vector<Nucleon> getNucleons()const;
    int getMOCount()const;
    Plane getPlane()const;
    double getShift()const;
    std::string getMOName()const;
    double getCornerX()const;
    double getCornerY()const;
    double getWidth()const;
    double getHeight()const;
    double getSameplDensity()const;

    std::string getInput()const;

    bool getApproxFlag()const;
    bool getOptFlag()const;
    bool getMOGraphFlag()const;
    bool getDGraphFlag()const;

private:

    // ALL IN ANGSTROMS HERE

    std::vector<Nucleon> _nucleons;
    int _MOCount{-1};
    Plane _plane{Plane::XY};
    double _shift{0.0};
    std::string _MOName{"Molecule"};
    double _cornerX{0};
    double _cornerY{0};
    double _width{0};
    double _height{0};
    double _sampleDencity{0.01};

    std::string _input;
    
    //flags

    bool _approxFarIntegrals{false};
    bool _HFOPT2{false};
    bool _MOGRAPHS{false};
    bool _DGRAPH{false};

private:

    std::string readHeader( std::ifstream & file);
    std::string readAtoms( std::ifstream & file);
    std::string readMOCount( std::ifstream & file);
    std::string readPlain(std::ifstream & file);
    std::string readShift(std::ifstream & file);
    std::string readMOName( std::ifstream & file);
    std::string readArea(std::ifstream & file);
    std::string readDensity( std::ifstream & file);

    std::string read(std::ifstream & file);

    int findInNucleonTab(const std::string & a);

};