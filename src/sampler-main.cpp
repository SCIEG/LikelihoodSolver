//==================================================================================================
// Name        : sampler-main.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include <iostream>
#include <string>

#include "Configuration.h"
#include "utils/FileReaderUtil.h"
#include "utils/StringUtil.h"

int main(int argc, char* argv[]) {
    const std::string prog = argv[0];
    if (argc != 4) {
        std::cerr << "Usage:" << std::endl
                  << "    " << prog << " input.csv N samples" << std::endl
                  << std::endl
                  << " - N is the number of unknowns in the numerator of the LR. Thus, N+1"
                  << std::endl
                  << "   is the number of unknowns in the denominator." << std::endl
                  << " - samples is the number of suspects to sample." << std::endl;
        return -1;
    }

    // The executable file name is 'sampler' but on windows is 'sampler.exe'
    int fileNameLength = 7;
    if (prog.find(".exe") != std::string::npos) {
        fileNameLength += 4;
    }
    const std::string executablePath = prog.substr(0, prog.length() - fileNameLength);

    const std::string inputFile = argv[1];

    int n;
    if (!LabRetriever::ToInt(argv[2], &n) || n < 0) {
      std::cerr << "N must be a non-negative integer!" << std::endl;
      return -1;
    }

    int samples;
    if (!LabRetriever::ToInt(argv[3], &samples) || samples <= 0) {
      std::cerr << "samples must be a positive integer!" << std::endl;
      return -1;
    }

    return 0;
}
