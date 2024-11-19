#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <fstream>

extern "C" {
    #include "grackle_macros.h"
}


extern "C" void FORTRAN_NAME(calc_tdust_1d_g)(double* tdust, double* tgas, double* nh, double* gasgr,
                                              double* gamma_isrfa, double* isrf, int* itmask, double* trad,
                                              int* in, int* is, int* ie, int* j, int* k);

TEST(FuncCalcTDust1DTest, CalcTDust1DTestJson) {

    std::ifstream f("/Users/bignamic/Development/SPH-EXA/GrackleWorks/datasets/calc_tdust_1d_g.json");
    nlohmann::json data = nlohmann::json::parse(f);

    double tdust[data["inputs"]["tdust"].size()];
    transform(data["inputs"]["tdust"].begin(), data["inputs"]["tdust"].end(), tdust, [](const double &x) { return x; });

    double tgas[data["inputs"]["tgas"].size()];
    transform(data["inputs"]["tgas"].begin(), data["inputs"]["tgas"].end(), tgas, [](const double &x) { return x; });

    double nh[data["inputs"]["nh"].size()];
    transform(data["inputs"]["nh"].begin(), data["inputs"]["nh"].end(), nh, [](const double &x) { return x; });

    double gasgr[data["inputs"]["gasgr"].size()];
    transform(data["inputs"]["gasgr"].begin(), data["inputs"]["gasgr"].end(), gasgr, [](const double &x) { return x; });

    double gamma_isrfa = data["inputs"]["gamma_isrfa"];

    double isrf[data["inputs"]["isrf"].size()];
    transform(data["inputs"]["isrf"].begin(), data["inputs"]["isrf"].end(), isrf, [](const double &x) { return x; });

    int itmask[data["inputs"]["itmask"].size()];
    transform(data["inputs"]["itmask"].begin(), data["inputs"]["itmask"].end(), itmask, [](const int &x) { return int(x); });

    double trad = data["inputs"]["trad"];
    int in = data["inputs"]["in"];
    int is = data["inputs"]["is"];
    int ie = data["inputs"]["ie"];

    int j = data["inputs"]["j"];
    int k = data["inputs"]["k"];

    FORTRAN_NAME(calc_tdust_1d_g)(tdust, tgas, nh, gasgr, &gamma_isrfa,
                                  isrf, itmask, &trad, &in, &is, &ie, &j, &k);

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(tdust[i], data["outputs"]["tdust"][i]);
    }

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(tgas[i], data["outputs"]["tgas"][i]);
    }

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(nh[i], data["outputs"]["nh"][i]);
    }

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(gasgr[i], data["outputs"]["gasgr"][i]);
    }

    EXPECT_DOUBLE_EQ(gamma_isrfa, data["outputs"]["gamma_isrfa"]);

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(isrf[i], data["outputs"]["isrf"][i]);
    }

    for(int i = 0; i < in; i++) {
        EXPECT_EQ(itmask[i], int(data["outputs"]["itmask"][i]));
    }

    for(int i = 0; i < in; i++) {
        EXPECT_DOUBLE_EQ(trad, data["outputs"]["trad"]);
    }

    EXPECT_EQ(in, data["outputs"]["in"]);
    EXPECT_EQ(is, data["outputs"]["is"]);
    EXPECT_EQ(ie, data["outputs"]["ie"]);
    EXPECT_EQ(j, data["outputs"]["j"]);
    EXPECT_EQ(k, data["outputs"]["k"]);

}
