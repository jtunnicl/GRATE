/*******************
 *
 *
 *  GRATE 8
 *
 *  Command Line Interface
 *
 *
 *
*********************/
#include "model.h"
#include <iostream>
#include <string>
#include <stdexcept>


int main(int argc, char** argv) {
    // number of steps can be passed as a argument, otherwise default to 800
    int nsteps = 800;
    if (argc == 2) {
        try {
            nsteps = std::stoi(argv[1]);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Invalid argument: " << ia.what() << std::endl;
            std::cerr << "First argument should be integer number of steps" << std::endl;
            return 1;
        }
    }
    std::cout << "Running for " << nsteps << " steps" << std::endl;

    // model object to be populated when reading the input file
    Model *model;

    // default input file name
    std::string param_file = "Conway_Template.xml";

    // load xml input file
    std::cout << "Reading xml file: '" << param_file << "'" << std::endl;
    XMLDocument xml_params;
    if (xml_params.LoadFile(param_file.c_str()) != XML_SUCCESS) {
        std::cerr << "Error reading xml parameters:" << std::endl;
        std::cerr << xml_params.ErrorStr() << std::endl;
        return 1;
    }
    else {
        // file was loaded so initialise everything

        // get the root element of the XML document
        XMLElement *params_root = xml_params.FirstChildElement();
        if (params_root == NULL) {
            std::cerr << "Error getting root element from XML file" << std::endl;
            std::cerr << xml_params.ErrorStr() << std::endl;
            return 1;
        }
        else {
            // initialise components
            try {
                model = new Model(params_root, "GrateResults.txt");
            }
            catch (std::string msg) {
                std::cerr << "Error while initialising components: " << msg << std::endl;
                return 1;
            }
        }
    }

    // run the model
    std::cout << "Running model for " << nsteps << " steps..." << std::endl;
    for (int i = 0; i < nsteps; i++) {
        model->iteration();

        if (i % 100 == 0) {
            std::cout << "Step " << i << " (" << static_cast<double>(i) / nsteps * 100.0 << " %)" << std::endl;
        }
    }

    // free model object
    delete model;

    std::cout << "Finished!" << std::endl;
}
