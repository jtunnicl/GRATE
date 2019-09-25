#ifndef MODEL_H
#define MODEL_H

#include "riverprofile.h"
#include "hydro.h"
#include "sed.h"
#include "tinyxml2/tinyxml2.h"

using namespace tinyxml2;

class Model {
    public:
        Model(XMLElement* params_root, string out1);
        ~Model();
        void iteration();

        RiverProfile *rn;
        hydro *wl;
        sed *sd;
        int writeInterval;

    private:
        void stepTime();
        void writeResults(int count);
};

#endif
