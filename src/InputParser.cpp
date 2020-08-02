#include "InputParser.h"

void InputParser::parse(const std::string & input)
{
    std::ifstream file(input, std::ios::in);

    if(!file.is_open())
    {
        throw InputOpenException();
        return;
    }

    std::string tmp = this->readHeader(file);

    while(tmp.size() > 0)
    {
        if(tmp[0] == '*')
            tmp = tmp.substr(1, tmp.size()-1);

        if(tmp == "xyz")
            tmp = this->readAtoms(file);
        else if(tmp == "plain")
            tmp = this->readPlain(file);
        else if(tmp == "shift")
            tmp = this->readShift(file);
        else if(tmp == "mo_count")
            tmp = this->readMOCount(file);
        else if(tmp == "name")
            tmp = this->readMOName(file);
        else if(tmp == "area")
            tmp = this->readArea(file);
        else if(tmp == "density")
            tmp = this->readDensity(file);    
        else
            tmp = this->read(file);           
    }

    file.close();
    file.open(input, std::ios::in);
    std::string line;
    while(std::getline(file, line))
    {
        _input += line;
        _input += "\n";
    }
    file.close();
}

std::vector<Nucleon> InputParser::getNucleons()const
{
    return _nucleons;
}
int InputParser::getMOCount()const
{
    return _MOCount;
}
InputParser::Plane InputParser::getPlane()const
{
    return _plane;
}
double InputParser::getShift()const
{
    return _shift;
}
std::string InputParser::getMOName()const
{
    return _MOName;
}
double InputParser::getCornerX()const
{
    return _cornerX;
}
double InputParser::getCornerY()const
{
    return _cornerY;
}
double InputParser::getWidth()const
{
    return _width;
}
double InputParser::getHeight()const
{
    return _height;
}
double InputParser::getSameplDensity()const
{
    return _sampleDencity;
}
std::string InputParser::getInput()const
{
    return _input;
}

bool InputParser::getApproxFlag()const
{
    return _approxFarIntegrals;
}
bool InputParser::getOptFlag()const
{
    return _HFOPT2;
}
bool InputParser::getMOGraphFlag()const
{
    return _MOGRAPHS;
}
bool InputParser::getDGraphFlag()const
{
    return _DGRAPH;
}


std::string InputParser::readHeader(std::ifstream & file)
{
    auto tmp = this->read(file);
    if(tmp == "HF")
        ;
    else if(tmp == "HF_OPT2")
        _HFOPT2 = true;
    else
        throw InputFileException();
    
    while(tmp[0] != '*')
    {
        tmp = this->read(file);
        if(tmp == "MO_GRAPHS")
            _MOGRAPHS = true;
        else if(tmp == "DENSITY_GRAPH")
            _DGRAPH = true;
        else if(tmp == "APPROX")
            _approxFarIntegrals = true;
    }
    return tmp;
}
std::string InputParser::readAtoms(std::ifstream & file)
{
    auto tmp = this->read(file);
    while(this->findInNucleonTab(tmp) == -1)
        tmp = this->read(file);
    while(tmp[0] != '*')
    {
        Nucleon n;
        n.charge = this->findInNucleonTab(tmp);
        if(n.charge == -1)
            throw InputFileException();
        try{
            tmp = this->read(file);
            n.p.x = std::stod(tmp);
            tmp = this->read(file);
            n.p.y = std::stod(tmp);
            tmp = this->read(file);
            n.p.z = std::stod(tmp);
        }catch (const std::exception & e){
            throw InputFileException();
        }
        n.p.reloadSpherical();
        _nucleons.push_back(n);
        tmp = this->read(file);
    }
    return tmp;
}
std::string InputParser::readMOCount(std::ifstream & file)
{
    auto tmp = this->read(file);
    try{
        _MOCount = std::stod(tmp);
    }catch (const std::exception & e){
        throw InputFileException();
    }
    return this->read(file);
}

std::string InputParser::readPlain(std::ifstream & file)
{
    auto tmp = this->read(file);
    if(tmp == "xy")
        _plane = Plane::XY;
    else if (tmp == "xz")
        _plane = Plane::XZ;
    else if(tmp == "yz")
        _plane = Plane::YZ;
    else
        throw InputFileException();
      
    return this->read(file);
}
std::string InputParser::readShift(std::ifstream & file)
{
    auto tmp = this->read(file);
    try{
        _shift = std::stod(tmp);
    }catch (const std::exception & e){
        throw InputFileException();
    }
    return this->read(file);
}
std::string InputParser::readMOName(std::ifstream & file)
{
    auto tmp = this->read(file);
    _MOName = tmp;
    return this->read(file);
}
std::string InputParser::readArea(std::ifstream & file)
{
    auto tmp = this->read(file);
    try{
        _cornerX = std::stod(tmp);
        tmp = this->read(file);
         _cornerY = std::stod(tmp);
        tmp = this->read(file);
         _width = std::stod(tmp);
        tmp = this->read(file);
         _height = std::stod(tmp);
        tmp = this->read(file);

    }catch (const std::exception & e){
        throw InputFileException();
    }
    return this->read(file);
}
std::string InputParser::readDensity(std::ifstream & file)
{
    auto tmp = this->read(file);
    try{
        _sampleDencity = std::stod(tmp);
    }catch (const std::exception & e){
        throw InputFileException();
    }
    return this->read(file);
}

std::string InputParser::read(std::ifstream & file)
{
    std::string res;
    file >> res;
    return res;
}

int InputParser::findInNucleonTab(const std::string & a)
{
    for (int i = 0; i < NUCLEON_COUNT + 1; i++)
    {
        if(NUCLEON_TAB[i] == a)
        return i;
    }
    return -1;
}