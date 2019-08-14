#include "pzreal.h"

void ProdT(REAL *v1, REAL *v2, REAL *mat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            mat[i * 3 + j] = v1[i] * v2[j];
        }
    }
}

void MatrixByScalar(REAL val, REAL *mat1, REAL *matres) {
    for(int i = 0; i < 3 * 3; i++) {
        matres[i] = val * mat1[i];
    }
}

void FromMatToVoight(REAL *mat, REAL *voi) {
    int k = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            voi[k++] = mat[i * 3 + j];
        }
    }
}

void TangentOperator(REAL *gradient, REAL *stress_eigenvalues, REAL *strain_eigenvalues, REAL *strain_eigenvectors, REAL *Tangent, REAL G, REAL lambda){

    //Montando a matriz tangente
    unsigned int kival[] = {0, 0, 0, 1, 1, 2};
    unsigned int kjval[] = {0, 1, 2, 1, 2, 2};

    // Coluna da matriz tangente
    for (unsigned int k = 0; k < 6; ++k) {
        const unsigned int ki = kival[k];
        const unsigned int kj = kjval[k];
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                REAL temp;
                if (ki == kj) {
                    temp = 2 * G * strain_eigenvectors[j * 3 + kj] * strain_eigenvectors[j * 3 + ki];
                    temp += lambda;
                } else {
                    temp = 2 * G  * strain_eigenvectors[j * 3 + kj] * strain_eigenvectors[j * 3 + ki];
                    temp *= 1.;
                }
                for (int l = 0; l < 6; ++l) {
                    const unsigned int li = kival[l];
                    const unsigned int lj = kjval[l];
                    Tangent[l * 6 + k] += temp * gradient[i * 3 + j] * strain_eigenvectors[i * 3 + li] * strain_eigenvectors[i * 3 + lj];
                }/// l
            }///j
        }///i
    }///k

//    for(int i = 0; i < 6; i++) {
//        for(int j = 0; j < 6; j++) {
//        std::cout << Tangent[i * 6 + j] << "\t";
//        }
//        std::cout << std::endl;
//    }


    REAL deigensig = 0., deigeneps = 0.;
    REAL tempMat[3 * 3];
    REAL temp_mat[3 * 3];
    REAL ColCorrV[6];

    for(int i = 0; i < 3 * 3; i++) tempMat[i] = 0;

    // Correction of the eigenvectors variation
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = i + 1; j < 3; ++j) {
            deigeneps = strain_eigenvalues[i] - strain_eigenvalues[j];
            deigensig = stress_eigenvalues[i] - stress_eigenvalues[j];

            REAL factor = 0.;
            if (!((fabs(deigeneps) < 1.e-12) || (deigeneps < 0.0))) {
                factor = deigensig / deigeneps;
            } else {
                factor = 0.5 * G * (gradient[i * 3 + i] - gradient[i * 3 + j] - gradient[j * 3 + i] + gradient[j * 3 + j]); // expression C.20
            }

            ProdT(&strain_eigenvectors[3 * i], &strain_eigenvectors[3 * j], temp_mat);
            for (unsigned int it = 0; it < 3; ++it) {
                for (unsigned int jt = 0; jt < 3; ++jt) {
                    tempMat[it * 3 + jt] += temp_mat[it * 3 + jt];
                }
            }

            ProdT(&strain_eigenvectors[3 * j], &strain_eigenvectors[3 * i], temp_mat);
            for (unsigned int it = 0; it < 3; ++it) {
                for (unsigned int jt = 0; jt < 3; ++jt) {
                    tempMat[it * 3 + jt] += temp_mat[it * 3 + jt];
                }
            }

            // expression C.14
            for (unsigned int k = 0; k < 6; ++k) {
                const unsigned int ki = kival[k];
                const unsigned int kj = kjval[k];
                if (ki == kj) {
                    REAL val = strain_eigenvectors[j * 3 + ki] * strain_eigenvectors[i * 3 + kj] * factor;
                    MatrixByScalar(val, tempMat, temp_mat);
                } else {
                    REAL val = (strain_eigenvectors[j * 3 + ki] * strain_eigenvectors[i * 3 + kj] + strain_eigenvectors[j * 3 + kj] * strain_eigenvectors[i * 3 + ki]) * factor;
                    MatrixByScalar(val, tempMat, temp_mat);
                }
                FromMatToVoight(temp_mat, ColCorrV);
                for (int l = 0; l < 6; l++) {
                    Tangent[l * 6 + k] += ColCorrV[l];
                }
            }
        } // j
    } // i
//    for(int i = 0; i < 6; i++) {
//        for(int j = 0; j < 6; j++) {
//            std::cout << Tangent[i * 6 + j] << "\t";
//        }
//        std::cout << std::endl;
//    }


}