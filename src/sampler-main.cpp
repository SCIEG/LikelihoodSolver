//==================================================================================================
// Name        : sampler-main.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include <iostream>

int main(int argc, char* argv[]) {
  if (argc != 4) {
    const string usageString =
      "Usage: \n"
      "    " + argv[0] + " input.csv N samples\n"
      "\n"
      " - N is the number of unknowns in the numerator of the LR. Thus, N+1 is the number of\n"
      "   unknowns in the denominator.\n"
      " - samples is the number of suspects to sample.";
    std::cout << usageString << std::endl;
    return -1;
  }

  string executablePath = string(argv[0]);

  // The executable file name is 'sampler' but on windows is 'sampler.exe'
  int fileNameLength = 7;
  if (executablePath.find(".exe") != string::npos) {
    fileNameLength += 4;
  }

  // Parse N and samples as ints.
  
  return 0;
}
