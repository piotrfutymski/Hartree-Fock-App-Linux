#include <iostream>
#include "App.h"

int main(int argc, char** argv) {
    
    App app;
    try{
        app.run(argv[1], argv[2]);
    }catch (const std::exception & e){
        if(argc != 3)
            std::cout<< "Incorrect arguments number: Give input file name and output file name"<<std::endl;
        std::cout<< "Exception was thrown: "<<e.what()<<std::endl;
    }


    return 0;
}
