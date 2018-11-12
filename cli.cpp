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


int main() {
    const int nsteps = 800;  // TODO: number of steps (or end time) should be parameter
    Model *model;

    // default input file name
    std::string param_file = "test_out.xml";

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
                model = new Model(params_root);
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
    }

    // free model object
    delete model;

    std::cout << "Finished!" << std::endl;
}
