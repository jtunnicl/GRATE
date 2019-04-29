// file to test the GrateTime object

#include "gratetime.h"
#include <iostream>


int main() {
    // create objects
    GrateTime t1 = GrateTime(2018, 8, 2, 11, 37, 21);
    GrateTime t2 = t1;
    
    // add some time
    int secsToAdd = 213;
    t2.addSecs(secsToAdd);

    // compute difference
    int diff = t1.secsTo(t2);
    if (diff != secsToAdd) {
        std::cerr << "Problem adding time and/or computing differences with GrateTime" << std::endl;
        return 1;
    }

    return 0;
}
