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
    
    minX -=0.8;
    maxX +=0.8;
    minY -=0.8;
    maxY +=0.8;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
        _xPoints = _yPoints;
    else
        _yPoints = _xPoints;
    

    _xMin = minX;
    _yMin = minY;

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
    
    minX -=0.8;
    maxX +=0.8;
    minY -=0.8;
    maxY +=0.8;

    _xPoints = (maxX-minX)/parser.getSameplDensity();
    _yPoints = (maxY-minY)/parser.getSameplDensity();
    if(_xPoints < _yPoints)
        _xPoints = _yPoints;
    else
        _yPoints = _xPoints;
    _xMin = minX;
    _yMin = minY;


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

void OutputCreator::createScriptOutputFile(Mol & mol, InputParser & parser)
{
    std::fstream file("int/pyinp.txt", std::ios::trunc | std::ios::out);
    file << mol.getMOcount()<<"\n";
    file << _xMin<<"\n"<<_yMin<<"\n"<<_xPoints<<"\n"<<_yPoints<<"\n"<<parser.getSameplDensity()<<"\n";
    file << parser.getMOName()<<"\n";
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
    
    file.close();
}