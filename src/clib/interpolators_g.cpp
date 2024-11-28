#include <algorithm>
#include <vector>
#include <cstdint>
#include <cmath>

void interpolate_1D_g(
    double input1,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1,
    double dgridPar1,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1;
    double slope;

    // Calcolo dell'indice di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1) ));

    // Interpolazione sul parametro 1
    slope = (dataField[index1 + 1] - dataField[index1]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + dataField[index1];
}

void interpolate_2D_g(
    double input1, double input2,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, double dgridPar2,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1, index2, int_index;
    int q;
    double slope, value2[2];

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index2 = std::min(gridDim[1] - 2,
              std::max(int64_t(0),
              int64_t((input2 - gridPar2[0]) / dgridPar2)));

    for (q = 0; q < 2; ++q)
    {
        // Interpolazione sul parametro 2
        int_index = (index1 + q) * gridDim[1] + index2;

        slope = (dataField[int_index + 1] - dataField[int_index]) /
                (gridPar2[index2 + 1] - gridPar2[index2]);

        value2[q] = (input2 - gridPar2[index2]) * slope + dataField[int_index];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}


void interpolate_3D_g(
    double input1, double input2, double input3,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, double dgridPar2,
    const std::vector<double>& gridPar3, double dgridPar3,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1, index2, index3, int_index;
    int q, w;
    double slope, value3[2], value2[2];

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index2 = std::min(gridDim[1] - 2,
              std::max(int64_t(0),
              int64_t((input2 - gridPar2[0]) / dgridPar2)));
    index3 = std::min(gridDim[2] - 2,
              std::max(int64_t(0),
              int64_t((input3 - gridPar3[0]) / dgridPar3)));

    for (q = 0; q < 2; ++q)
    {
        for (w = 0; w < 2; ++w)
        {
            // Interpolazione sul parametro 3
            int_index = ((index1 + q) * gridDim[1] + (index2 + w)) * gridDim[2] + index3;

            slope = (dataField[int_index + 1] - dataField[int_index]) /
                    (gridPar3[index3 + 1] - gridPar3[index3]);

            value3[w] = (input3 - gridPar3[index3]) * slope + dataField[int_index];
        }

        // Interpolazione sul parametro 2
        slope = (value3[1] - value3[0]) /
                (gridPar2[index2 + 1] - gridPar2[index2]);

        value2[q] = (input2 - gridPar2[index2]) * slope + value3[0];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}

void interpolate_3Dz_g(
    double input1, double input2, double input3,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, int64_t index2,
    const std::vector<double>& gridPar3, double dgridPar3,
    int64_t dataSize,
    const std::vector<double>& dataField,
    bool end_int,
    double& value)
{
    // Variabili locali
    int64_t index1, index3, int_index;
    int q, w;
    double slope, value3[2], value2[2];

    if (end_int)
    {
        // Chiama interpolate_2Df3D_g (da implementare)
        interpolate_2Df3D_g(
            input1, input3, gridDim,
            gridPar1, dgridPar1,
            index2,
            gridPar3, dgridPar3,
            dataSize, dataField,
            value
        );
        return;
    }

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index3 = std::min(gridDim[2] - 2,
              std::max(int64_t(0),
              int64_t((input3 - gridPar3[0]) / dgridPar3)));

    for (q = 0; q < 2; ++q)
    {
        for (w = 0; w < 2; ++w)
        {
            // Interpolazione sul parametro 3
            int_index = ((index1 + q) * gridDim[1] + (index2 + w)) * gridDim[2] + index3;

            slope = (dataField[int_index + 1] - dataField[int_index]) /
                    (gridPar3[index3 + 1] - gridPar3[index3]);

            value3[w] = (input3 - gridPar3[index3]) * slope + dataField[int_index];
        }

        // Interpolazione sul parametro 2 (redshift)
        slope = (value3[1] - value3[0]) /
                std::log((1.0 + gridPar2[index2 + 1]) / (1.0 + gridPar2[index2]));

        value2[q] = std::log((1.0 + input2) / (1.0 + gridPar2[index2])) * slope + value3[0];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}

void interpolate_3Dz_g(
    double input1, double input2, double input3,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, int64_t index2,
    const std::vector<double>& gridPar3, double dgridPar3,
    int64_t dataSize,
    const std::vector<double>& dataField,
    bool end_int,
    double& value)
{
    // Variabili locali
    int64_t index1, index3, int_index;
    int q, w;
    double slope, value3[2], value2[2];

    if (end_int)
    {
        // Chiama interpolate_2Df3D_g (da implementare)
        interpolate_2Df3D_g(
            input1, input3, gridDim,
            gridPar1, dgridPar1,
            index2,
            gridPar3, dgridPar3,
            dataSize, dataField,
            value
        );
        return;
    }

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index3 = std::min(gridDim[2] - 2,
              std::max(int64_t(0),
              int64_t((input3 - gridPar3[0]) / dgridPar3)));

    for (q = 0; q < 2; ++q)
    {
        for (w = 0; w < 2; ++w)
        {
            // Interpolazione sul parametro 3
            int_index = ((index1 + q) * gridDim[1] + (index2 + w)) * gridDim[2] + index3;

            slope = (dataField[int_index + 1] - dataField[int_index]) /
                    (gridPar3[index3 + 1] - gridPar3[index3]);

            value3[w] = (input3 - gridPar3[index3]) * slope + dataField[int_index];
        }

        // Interpolazione sul parametro 2 (redshift)
        slope = (value3[1] - value3[0]) /
                std::log((1.0 + gridPar2[index2 + 1]) / (1.0 + gridPar2[index2]));

        value2[q] = std::log((1.0 + input2) / (1.0 + gridPar2[index2])) * slope + value3[0];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}

void interpolate_2Df3D_g(
    double input1, double input3,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    int64_t index2,
    const std::vector<double>& gridPar3, double dgridPar3,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1, index3, int_index;
    int q;
    double slope, value3[2];

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2,
              std::max(int64_t(0),
              int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index3 = std::min(gridDim[2] - 2,
              std::max(int64_t(0),
              int64_t((input3 - gridPar3[0]) / dgridPar3)));

    for (q = 0; q < 2; ++q)
    {
        // Interpolazione sul parametro 3
        int_index = ((index1 + q) * gridDim[1] + index2) * gridDim[2] + index3;

        slope = (dataField[int_index + 1] - dataField[int_index]) /
                (gridPar3[index3 + 1] - gridPar3[index3]);

        value3[q] = (input3 - gridPar3[index3]) * slope + dataField[int_index];
    }

    // Interpolazione sul parametro 1
    slope = (value3[1] - value3[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value3[0];
}

void interpolate_4D_g(
    double input1, double input2, double input3, double input4,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, double dgridPar2,
    const std::vector<double>& gridPar3, double dgridPar3,
    const std::vector<double>& gridPar4, double dgridPar4,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1, index2, index3, index4, int_index;
    int q, w, e;
    double slope, value4[2], value3[2], value2[2];

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2, std::max(int64_t(0),
             int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index2 = std::min(gridDim[1] - 2, std::max(int64_t(0),
             int64_t((input2 - gridPar2[0]) / dgridPar2)));
    index3 = std::min(gridDim[2] - 2, std::max(int64_t(0),
             int64_t((input3 - gridPar3[0]) / dgridPar3)));
    index4 = std::min(gridDim[3] - 2, std::max(int64_t(0),
             int64_t((input4 - gridPar4[0]) / dgridPar4)));

    for (q = 0; q < 2; ++q)
    {
        for (w = 0; w < 2; ++w)
        {
            for (e = 0; e < 2; ++e)
            {
                // Interpolazione sul parametro 4
                int_index = (((index1 + q) * gridDim[1] + (index2 + w)) * gridDim[2] + (index3 + e)) * gridDim[3] + index4;

                slope = (dataField[int_index + 1] - dataField[int_index]) /
                        (gridPar4[index4 + 1] - gridPar4[index4]);

                value4[e] = (input4 - gridPar4[index4]) * slope + dataField[int_index];
            }

            // Interpolazione sul parametro 3
            slope = (value4[1] - value4[0]) /
                    (gridPar3[index3 + 1] - gridPar3[index3]);

            value3[w] = (input3 - gridPar3[index3]) * slope + value4[0];
        }

        // Interpolazione sul parametro 2
        slope = (value3[1] - value3[0]) /
                (gridPar2[index2 + 1] - gridPar2[index2]);

        value2[q] = (input2 - gridPar2[index2]) * slope + value3[0];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}

void interpolate_5D_g(
    double input1, double input2, double input3, double input4, double input5,
    const std::vector<int64_t>& gridDim,
    const std::vector<double>& gridPar1, double dgridPar1,
    const std::vector<double>& gridPar2, double dgridPar2,
    const std::vector<double>& gridPar3, double dgridPar3,
    const std::vector<double>& gridPar4, double dgridPar4,
    const std::vector<double>& gridPar5, double dgridPar5,
    int64_t dataSize,
    const std::vector<double>& dataField,
    double& value)
{
    // Variabili locali
    int64_t index1, index2, index3, index4, index5, int_index;
    int q, w, e, r;
    double slope, value5[2], value4[2], value3[2], value2[2];

    // Calcolo degli indici di interpolazione
    index1 = std::min(gridDim[0] - 2, std::max(int64_t(0),
             int64_t((input1 - gridPar1[0]) / dgridPar1)));
    index2 = std::min(gridDim[1] - 2, std::max(int64_t(0),
             int64_t((input2 - gridPar2[0]) / dgridPar2)));
    index3 = std::min(gridDim[2] - 2, std::max(int64_t(0),
             int64_t((input3 - gridPar3[0]) / dgridPar3)));
    // index4 may require special handling (e.g., bisection method), omitted for brevity
    index4 = std::min(gridDim[3] - 2, std::max(int64_t(0),
             int64_t((input4 - gridPar4[0]) / dgridPar4)));
    index5 = std::min(gridDim[4] - 2, std::max(int64_t(0),
             int64_t((input5 - gridPar5[0]) / dgridPar5)));

    for (q = 0; q < 2; ++q)
    {
        for (w = 0; w < 2; ++w)
        {
            for (e = 0; e < 2; ++e)
            {
                for (r = 0; r < 2; ++r)
                {
                    // Interpolazione sul parametro 5
                    int_index = ((((index1 + q) * gridDim[1] + (index2 + w)) * gridDim[2] + (index3 + e)) * gridDim[3] + (index4 + r)) * gridDim[4] + index5;

                    slope = (dataField[int_index + 1] - dataField[int_index]) /
                            (gridPar5[index5 + 1] - gridPar5[index5]);

                    value5[r] = (input5 - gridPar5[index5]) * slope + dataField[int_index];
                }

                // Interpolazione sul parametro 4
                slope = (value5[1] - value5[0]) /
                        (gridPar4[index4 + 1] - gridPar4[index4]);

                value4[e] = (input4 - gridPar4[index4]) * slope + value5[0];
            }

            // Interpolazione sul parametro 3
            slope = (value4[1] - value4[0]) /
                    (gridPar3[index3 + 1] - gridPar3[index3]);

            value3[w] = (input3 - gridPar3[index3]) * slope + value4[0];
        }

        // Interpolazione sul parametro 2
        slope = (value3[1] - value3[0]) /
                (gridPar2[index2 + 1] - gridPar2[index2]);

        value2[q] = (input2 - gridPar2[index2]) * slope + value3[0];
    }

    // Interpolazione sul parametro 1
    slope = (value2[1] - value2[0]) /
            (gridPar1[index1 + 1] - gridPar1[index1]);

    value = (input1 - gridPar1[index1]) * slope + value2[0];
}
