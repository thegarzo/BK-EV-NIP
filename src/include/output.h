
#ifndef OUT_H
#define OUT_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

namespace Output{

    void printProgress2(double percentage1, double percentage2) {
        int val1 = (int) (percentage1 * 100);
        int lpad1 = (int) (percentage1 * PBWIDTH);
        int rpad1 = PBWIDTH - lpad1;

            int val2 = (int) (percentage2 * 100);
        int lpad2 = (int) (percentage2 * PBWIDTH);
        int rpad2 = PBWIDTH - lpad2;

            printf("\r%3d%% [%.*s%*s]", val1, lpad1, PBSTR, rpad1, "");
            printf("   %3d%% [%.*s%*s]", val2, lpad2, PBSTR, rpad2, "");
        fflush(stdout);
    }
}

#endif

