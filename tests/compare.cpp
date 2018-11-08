/*
 * Compare two output files
 */

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>


// Numbers a and b are the same if (from numpy.isclose):
//   abs(a -b) <= (atol + rtol * abs(b))
#define ATOL 1e-08
#define RTOL 1e-05


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: compare <REF_FILE> <NEW_FILE>" << std::endl;
        return 1;
    }
    std::string refFileName(argv[1]);
    std::string newFileName(argv[2]);

    std::cout << "Comparing Grate output files" << std::endl;

    // open reference file
    std::ifstream refFile(refFileName);
    if (refFile.fail()) {
        std::cerr << "Error opening file: " << refFileName << std::endl;
        return 1;
    }
    std::cout << "Opened ref file: " << refFileName << std::endl;

    // open new file
    std::ifstream newFile(newFileName);
    if (newFile.fail()) {
        std::cerr << "Error opening file: " << newFileName << std::endl;
        return 1;
    }
    std::cout << "Opened new file: " << newFileName << std::endl;

    // loop over header lines
    std::string newLine, refLine;
    bool gotHeader = false;
    int lineNumber = 0;
    while (not gotHeader) {
        lineNumber++; 

        // ref line
        if (not std::getline(refFile, refLine)) {
            std::cerr << "Error reading line in ref file line " << lineNumber << std::endl;
            return 1;
        }

        // check if the line starts with Count
        std::istringstream refiss(refLine);
        std::string token;
        refiss >> token;
        if (token.compare("Count:") == 0) {
            // header ends when we find a Count line
            gotHeader = true;
        }
        else {
            // new line
            if (not std::getline(newFile, newLine)) {
                std::cerr << "Error reading line in new file line " << lineNumber << std::endl;
                return 1;
            }

            // compare
            if (refLine.compare(newLine)) {
                std::cerr << "Lines differ:" << std::endl << refLine << std::endl << newLine << std::endl;
                return 1;
            }
        }
    }
    std::cout << "Header matches!" << std::endl;

    // loop over the rest of the file
    do {
        // we have refLine already; get newLine to match
        if (not std::getline(newFile, newLine)) {
            std::cerr << "Error reading line in new file " << lineNumber << std::endl;
            return 1;
        }

        std::istringstream refiss(refLine);
        std::istringstream newiss(newLine);
        std::string reftok;
        refiss >> reftok;
        if (reftok.compare("Count:") == 0) {
            // check count line matches directly
            if (refLine.compare(newLine)) {
                std::cerr << "Lines differ (L" << lineNumber << "):" << std::endl << refLine << std::endl << newLine << std::endl;
                return 1;
            }
        }
        else {
            // for other lines, we split into tokens and look at relative differences (allowed to be small difference)
            std::istringstream refstream(refLine);
            std::istringstream newstream(newLine);
            double newval, refval;
            while (refstream >> refval) {
                if (not (newstream >> newval)) {
                    std::cerr << "Ran out of tokens: " << lineNumber << std::endl;
                    return 1;
                }
                
                // formula from numpy.isclose
                if (abs(newval - refval) > (ATOL + RTOL * abs(refval))) {
                    std::cerr << "Values differ (L" << lineNumber << "): " << newval << " vs " << refval << std::endl;
                    return 1;
                }
            }
        }

        lineNumber++;
    }
    while (std::getline(refFile, refLine));
    std::cout << "Body matches!" << std::endl;

    // check if any extra lines left at the end of the new file
    if (std::getline(newFile, newLine)) {
        std::cerr << "New file has extra lines at the end" << std::endl;
        return 1;
    }

    return 0;
}
