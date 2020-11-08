#include "OutputCreator.h"


void OutputCreator::createOutputFile(Logger & logger, Mol & mol, InputParser & parser)
{
    logger.log("--------------------------------------------------------------\n");
    logger.log("-------------------- Results for molecule --------------------\n");
    logger.log("--------------------------------------------------------------\n");
    logger.log("-------------------------- Input -----------------------------\n");
    logger.log(parser.getInput());
    logger.log("--------------------------------------------------------------\n");
    logger.log("Molecule name:\t");
    logger.log(parser.getMOName());
    logger.log("\n");
    logger.log("Number of filled orbitals:\t");
    logger.log(std::to_string(mol.getMOcount()));
    logger.log("\n");
    logger.log("Computation done in:\t");
    logger.log(std::to_string(mol.getTotalTime()));
    logger.log("s\n");
    if(mol.doesNotConv())
    {
        logger.log("Energy does not converges!\n");
        return;
    }
    logger.log("\n");
    logger.log("MOLECULE TOTAL ENERGY:\t");
    logger.log(std::to_string(mol.getTotalEnergy()));
    logger.log("\n");
    logger.log("Nucleon repulsion energy:\t");
    logger.log(std::to_string(mol.getRepulsionEnergy()));
    logger.log("\n");
    logger.log("Electron energy:\t");
    logger.log(std::to_string(mol.getElectronicEnergy()));
    logger.log("\n");
    logger.log("Orbital energies:\t");
    for(int i = 0; i < mol.getMOcount(); ++i)
    {
        logger.log(std::to_string(mol.getOribtalEnergy(i)));
        logger.log("\t");
    }
    logger.log("\n");
    logger.log("--------------------------------------------------------------\n");
}

void OutputCreator::createOutputFileCompare(Logger & logger, Mol & molA, Mol & B, InputParser & parser)
{
    logger.log("--------------------------------------------------------------\n");
    logger.log("-------------------- Results for molecule --------------------\n");
    logger.log("--------------------------------------------------------------\n");
    logger.log("-------------------------- Input -----------------------------\n");
    logger.log(parser.getInput());
    logger.log("--------------------------------------------------------------\n");
    logger.log("MOLECULE TOTAL ENERGY A:\t");
    logger.log(std::to_string(molA.getTotalEnergy()));
    logger.log("\n");
    logger.log("MOLECULE TOTAL ENERGY B:\t");
    logger.log(std::to_string(B.getTotalEnergy()));
    logger.log("\n");
}

void OutputCreator::createMOPlainFiles(Mol & mol, InputParser & parser)
{
    double minX =1000;
    double minY =1000;
    double maxX =-1000;
    double maxY =-1000;

    for (int i = 0; i < parser.getNucleons().size(); i++)
    {
        if(parser.getPlane()==InputParser::Plane::XY || parser.getPlane()==InputParser::Plane::XZ)
        {
            if(parser.getNucleons()[i].p.x < minX)
                minX = parser.getNucleons()[i].p.x;
            if(parser.getNucleons()[i].p.x > maxX)
                maxX = parser.getNucleons()[i].p.x;
        }
        else
        {
            if(parser.getNucleons()[i].p.y < minX)
                minX = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxX)
                maxX = parser.getNucleons()[i].p.y;
        }
        
        if(parser.getPlane()==InputParser::Plane::XY)
        {
            if(parser.getNucleons()[i].p.y < minY)
                minY = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxY)
                maxY = parser.getNucleons()[i].p.y;
        }
        else
        {
            if(parser.getNucleons()[i].p.z < minY)
                minY = parser.getNucleons()[i].p.z;
            if(parser.getNucleons()[i].p.z > maxY)
                maxY = parser.getNucleons()[i].p.z;
        }   

    }
    
    minX -=1.0;
    maxX +=1.0;
    minY -=1.0;
    maxY +=1.0;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
    {
        _xMin = minX - (_yPoints-_xPoints)*parser.getSameplDensity()/2;
        minX = _xMin;
        _xPoints = _yPoints;
        _yMin = minY;
    }       
    else
    {   
        _yMin = minY - (_xPoints-_yPoints)*parser.getSameplDensity()/2;
        minY = _yMin;
        _yPoints = _xPoints;
        _xMin = minX;
    }

    for (int i = 0; i < mol.getMOcount(); i++)
	{
		std::fstream file("int/plain_"+ std::to_string(i) + ".txt", std::ios::trunc | std::ios::out);
		for(int j = 0; j < _xPoints; j++)
		{
			for(int k = 0; k < _yPoints; k++)
			{
                double x,y,z;
                if(parser.getPlane()==InputParser::Plane::XY)
                {
                    x = minX + j * parser.getSameplDensity();
                    y = minY + k * parser.getSameplDensity();
                    z = parser.getShift();
                }
                else if(parser.getPlane()==InputParser::Plane::XZ)
                {
                    x = minX + j * parser.getSameplDensity();
                    z = minY + k * parser.getSameplDensity();
                    y = parser.getShift();
                }
                else if(parser.getPlane()==InputParser::Plane::YZ)
                {
                    y = minX + j * parser.getSameplDensity();
                    z = minY + k * parser.getSameplDensity();
                    x = parser.getShift();
                }

				auto v = mol.countMolecularFunction(i, { x,y,z });
				file << v <<"\n";
			}
		}
		file.close();
	}
}

void OutputCreator::createMOPlainFilesCompare(Mol & molA, Mol & molB, InputParser & parser)
{
    double minX =1000;
    double minY =1000;
    double maxX =-1000;
    double maxY =-1000;

    for (int i = 0; i < parser.getNucleons().size(); i++)
    {
        if(parser.getPlane()==InputParser::Plane::XY || parser.getPlane()==InputParser::Plane::XZ)
        {
            if(parser.getNucleons()[i].p.x < minX)
                minX = parser.getNucleons()[i].p.x;
            if(parser.getNucleons()[i].p.x > maxX)
                maxX = parser.getNucleons()[i].p.x;
        }
        else
        {
            if(parser.getNucleons()[i].p.y < minX)
                minX = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxX)
                maxX = parser.getNucleons()[i].p.y;
        }
        
        if(parser.getPlane()==InputParser::Plane::XY)
        {
            if(parser.getNucleons()[i].p.y < minY)
                minY = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxY)
                maxY = parser.getNucleons()[i].p.y;
        }
        else
        {
            if(parser.getNucleons()[i].p.z < minY)
                minY = parser.getNucleons()[i].p.z;
            if(parser.getNucleons()[i].p.z > maxY)
                maxY = parser.getNucleons()[i].p.z;
        }   

    }
    
    minX -=1.0;
    maxX +=1.0;
    minY -=1.0;
    maxY +=1.0;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
    {
        _xMin = minX - (_yPoints-_xPoints)*parser.getSameplDensity()/2;
        minX = _xMin;
        _xPoints = _yPoints;
        _yMin = minY;
    }       
    else
    {   
        _yMin = minY - (_xPoints-_yPoints)*parser.getSameplDensity()/2;
        minY = _yMin;
        _yPoints = _xPoints;
        _xMin = minX;
    }

    for (int i = 0; i < molA.getMOcount(); i++)
	{
		std::fstream file("int/plain_"+ std::to_string(i) + ".txt", std::ios::trunc | std::ios::out);
		for(int j = 0; j < _xPoints; j++)
		{
			for(int k = 0; k < _yPoints; k++)
			{
                double x,y,z;
                if(parser.getPlane()==InputParser::Plane::XY)
                {
                    x = minX + j * parser.getSameplDensity();
                    y = minY + k * parser.getSameplDensity();
                    z = parser.getShift();
                }
                else if(parser.getPlane()==InputParser::Plane::XZ)
                {
                    x = minX + j * parser.getSameplDensity();
                    z = minY + k * parser.getSameplDensity();
                    y = parser.getShift();
                }
                else if(parser.getPlane()==InputParser::Plane::YZ)
                {
                    y = minX + j * parser.getSameplDensity();
                    z = minY + k * parser.getSameplDensity();
                    x = parser.getShift();
                }

				auto v = molA.countMolecularFunction(i, { x,y,z }) - molB.countMolecularFunction(i, { x,y,z });
				file << v <<"\n";
			}
		}
		file.close();
	}
}

void OutputCreator::createDPlainFile(Mol & mol, InputParser & parser)
{
    double minX =1000;
    double minY =1000;
    double maxX =-1000;
    double maxY =-1000;

    for (int i = 0; i < parser.getNucleons().size(); i++)
    {
        if(parser.getPlane()==InputParser::Plane::XY || parser.getPlane()==InputParser::Plane::XZ)
        {
            if(parser.getNucleons()[i].p.x < minX)
                minX = parser.getNucleons()[i].p.x;
            if(parser.getNucleons()[i].p.x > maxX)
                maxX = parser.getNucleons()[i].p.x;
        }
        else
        {
            if(parser.getNucleons()[i].p.y < minX)
                minX = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxX)
                maxX = parser.getNucleons()[i].p.y;
        }
        
        if(parser.getPlane()==InputParser::Plane::XY)
        {
            if(parser.getNucleons()[i].p.y < minY)
                minY = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxY)
                maxY = parser.getNucleons()[i].p.y;
        }
        else
        {
            if(parser.getNucleons()[i].p.z < minY)
                minY = parser.getNucleons()[i].p.z;
            if(parser.getNucleons()[i].p.z > maxY)
                maxY = parser.getNucleons()[i].p.z;
        }   

    }
    
    minX -=1.0;
    maxX +=1.0;
    minY -=1.0;
    maxY +=1.0;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
    {
        _xMin = minX - (_yPoints-_xPoints)*parser.getSameplDensity()/2;
        minX = _xMin;
        _xPoints = _yPoints;
        _yMin = minY;
    }       
    else
    {   
        _yMin = minY - (_xPoints-_yPoints)*parser.getSameplDensity()/2;
        minY = _yMin;
        _yPoints = _xPoints;
        _xMin = minX;
    }



	std::fstream file("int/dplain.txt", std::ios::trunc | std::ios::out);
	for(int j = 0; j < _xPoints; j++)
	{
		for(int k = 0; k < _yPoints; k++)
		{
            double x,y,z;
            if(parser.getPlane()==InputParser::Plane::XY)
            {
                x = minX + j * parser.getSameplDensity();
                y = minY + k * parser.getSameplDensity();
                z = parser.getShift();
            }
            else if(parser.getPlane()==InputParser::Plane::XZ)
            {
                x = minX + j * parser.getSameplDensity();
                z = minY + k * parser.getSameplDensity();
                y = parser.getShift();
            }
            else if(parser.getPlane()==InputParser::Plane::YZ)
            {
                y = minX + j * parser.getSameplDensity();
                z = minY + k * parser.getSameplDensity();
                x = parser.getShift();
            }

			auto v = mol.countDencity({ x,y,z });
			file << v <<"\n";
		}

	}
    file.close();
}

void OutputCreator::createDPlainFileCompare(Mol & molA,Mol & molB, InputParser & parser)
{
    double minX =1000;
    double minY =1000;
    double maxX =-1000;
    double maxY =-1000;

    for (int i = 0; i < parser.getNucleons().size(); i++)
    {
        if(parser.getPlane()==InputParser::Plane::XY || parser.getPlane()==InputParser::Plane::XZ)
        {
            if(parser.getNucleons()[i].p.x < minX)
                minX = parser.getNucleons()[i].p.x;
            if(parser.getNucleons()[i].p.x > maxX)
                maxX = parser.getNucleons()[i].p.x;
        }
        else
        {
            if(parser.getNucleons()[i].p.y < minX)
                minX = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxX)
                maxX = parser.getNucleons()[i].p.y;
        }
        
        if(parser.getPlane()==InputParser::Plane::XY)
        {
            if(parser.getNucleons()[i].p.y < minY)
                minY = parser.getNucleons()[i].p.y;
            if(parser.getNucleons()[i].p.y > maxY)
                maxY = parser.getNucleons()[i].p.y;
        }
        else
        {
            if(parser.getNucleons()[i].p.z < minY)
                minY = parser.getNucleons()[i].p.z;
            if(parser.getNucleons()[i].p.z > maxY)
                maxY = parser.getNucleons()[i].p.z;
        }   

    }
    
    minX -=1.0;
    maxX +=1.0;
    minY -=1.0;
    maxY +=1.0;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
    {
        _xMin = minX - (_yPoints-_xPoints)*parser.getSameplDensity()/2;
        minX = _xMin;
        _xPoints = _yPoints;
        _yMin = minY;
    }       
    else
    {   
        _yMin = minY - (_xPoints-_yPoints)*parser.getSameplDensity()/2;
        minY = _yMin;
        _yPoints = _xPoints;
        _xMin = minX;
    }



	std::fstream file("int/dplain.txt", std::ios::trunc | std::ios::out);
	for(int j = 0; j < _xPoints; j++)
	{
		for(int k = 0; k < _yPoints; k++)
		{
            double x,y,z;
            if(parser.getPlane()==InputParser::Plane::XY)
            {
                x = minX + j * parser.getSameplDensity();
                y = minY + k * parser.getSameplDensity();
                z = parser.getShift();
            }
            else if(parser.getPlane()==InputParser::Plane::XZ)
            {
                x = minX + j * parser.getSameplDensity();
                z = minY + k * parser.getSameplDensity();
                y = parser.getShift();
            }
            else if(parser.getPlane()==InputParser::Plane::YZ)
            {
                y = minX + j * parser.getSameplDensity();
                z = minY + k * parser.getSameplDensity();
                x = parser.getShift();
            }

			auto v = molA.countDencity({ x,y,z }) - molB.countDencity({ x,y,z });
			file << v <<"\n";
		}

	}
    file.close();
}

void OutputCreator::createScriptOutputFile(Mol & mol, InputParser & parser)
{
    std::fstream file("int/pyinp.txt", std::ios::trunc | std::ios::out);
    if(parser.getMOGraphFlag() && parser.getDGraphFlag())
    {
        file << "t\n";
    }
    else if(parser.getDGraphFlag())
    {
        file << "d\n";
    }
    else if(parser.getMOGraphFlag())
    {
        file << "p\n";
    }
    else
    {
        file << "n\n";
    }
    file << mol.getMOcount()<<"\n";
    file << _xMin<<"\n"<<_yMin<<"\n"<<_xPoints<<"\n"<<_yPoints<<"\n"<<parser.getSameplDensity()<<"\n";
    file << parser.getMOName()<<"\n";

    
    file.close();
}

void OutputCreator::createScriptOutputFileCompare(Mol & molA, Mol & molB, InputParser & parser)
{
    std::fstream file("int/pyinp.txt", std::ios::trunc | std::ios::out);
    file << "t\n";
    file << molA.getMOcount()<<"\n";
    file << _xMin<<"\n"<<_yMin<<"\n"<<_xPoints<<"\n"<<_yPoints<<"\n"<<parser.getSameplDensity()<<"\n";
    file << parser.getMOName()<<"\n";  
    file.close();
}