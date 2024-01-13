#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

// Function Prototypes
double calculateDistance(double Cn, double Cn0, double MS_min, double MS_max, double MS_Cn_min, double MS_Cn_max, double CritCnmin, double CritCnmax, double CritMsmin, double CritMsmax, double CritMsCnmin, double CritMsCnmax);
double calculateLineDistance(double value1, double value2, double center1, double center2, double centerRangeFactor);
double distanceToCenter(double value, double center);
int valueInRange(double value, double rangeStart, double rangeEnd);

// Function Prototypes for the calculations of values "Portance", "MS", and "Couple"
double calculateCn(const char* Type_masquage, double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail, double Q_can, double E_can, double D_can, double f_can, double m_can, double n_can, double D_og, double D2j, double D1j, double D2r, double D1r, double Q_int, double E_int, double D_int, double f_int, double m_int, double n_int);
double calculateCn0(const char* Type_masquage, double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail, double Q_can, double E_can, double D_can, double f_can, double m_can, double n_can, double D_og, double D2j, double D1j, double D2r, double D1r);
double calculateMS_min(double XCp, double XcgPlein, double D_ref);
double calculateMS_max(double XCp0, double XcgVide, double D_ref);
double calculateMS_Cn_min(double MS_min, double Cn);
double calculateMS_Cn_max(double MS_max, double Cn0);

// Function Prototypes for the calculations of the critical values for "Portance", "MS", and "Couple"
double calculateCritCnmin(const char* Type_fusee);
double calculateCritCnmax(const char* Type_fusee);
double calculateCritMsmin(const char* Type_fusee);
double calculateCritMsmax(const char* Type_fusee);
double calculateCritMsCnmin(const char* Type_fusee);
double calculateCritMsCnmax(const char* Type_fusee);

// This code only works for two-stage rockets with a bi-empennage
const char* Type_fusee = "Fusée expérimentale.";  // Rocket type
const char* Type_masquage_bi = "bi-empennage"; // Masking type for the entire rocket
const char* Type_masquage_ship = "Mono-empennage"; // Masking type for the ship
const char* withOrWithoutPropu = "sans propu"; // With or without propellant

//This code takes x arguments in this order: 
int main(int argc, char *argv[]) {
    // Input values
    if (argc < 5) {
        fprintf(stderr, "Not enough arguments.\n");
        return 1;
    }

    int currentStep = atoi(argv[1]);
    int Long_tot = atoi(argv[2]);
    int D_ref = atoi(argv[3]);
    char* Forme_ogive = argv[4];
    int Long_ogive = atoi(argv[5]);
    int D_og = atoi(argv[6]);
    char* Propu = argv[7];
    int XpropuRef = atoi(argv[8]);
    int Q_ail = atoi(argv[9]);
    int Q_can = atoi(argv[10]);
    int X_ail = atoi(argv[11]);
    int X_can = atoi(argv[12]);
    int D2j = atoi(argv[13]);

    // Formulas
    int D_ail = D_ref;
    int D_can = D_ref;
    int D1j = D_og;
    int D2r = D_og;
    int D1r = D2j;

    // Variable declarations that will stock the best values found
    double best_m_ail = 0, best_n_ail = 0, best_p_ail = 0, best_E_ail = 0;
    double best_m_can = 0, best_n_can = 0, best_p_can = 0, best_E_can = 0;
    double minDistanceGlobal = DBL_MAX;
    double Q_int = (Q_ail == Q_can) ? Q_ail : 0;
    double D_int = D_ail;

    // Critical values
    double CritCnmin = calculateCritCnmin(Type_fusee);
    double CritCnmax = calculateCritCnmax(Type_fusee);
    double CritMsmin = calculateCritMsmin(Type_fusee);
    double CritMsmax = calculateCritMsmax(Type_fusee);
    double CritMsCnmin = calculateCritMsCnmin(Type_fusee);
    double CritMsCnmax = calculateCritMsCnmax(Type_fusee);

    // Iterate over the range of values
    for (double E_ail = 100; E_ail <= 400; E_ail += currentStep) {
        for (double m_ail = 50; m_ail <= 300; m_ail += currentStep) {
            for (double n_ail = 50; n_ail <= m_ail; n_ail += currentStep) {
                for (double p_ail = 0; p_ail <= n_ail; p_ail += currentStep) {
                    for (double E_can = 100; E_can <= 400; E_can += currentStep) {
                        for (double m_can = 50; m_can <= 300; m_can += currentStep) {
                            for (double n_can = 50; n_can <= m_can; n_can += currentStep) {
                                for (double p_can = 0; p_can <= n_can; p_can += currentStep) {
                                    // Calculations
                                    double E_int;
                                    if ((D_can / 2.0 + E_can) <= (D_ail / 2.0)) {
                                        E_int = 0;
                                    } else if ((D_can / 2.0 + E_can) >= (D_ail / 2.0 + E_ail)) {
                                        E_int = E_ail;
                                    } else {
                                        E_int = (D_can / 2.0 + E_can) - (D_ail / 2.0);
                                    }

                                    double m_int = m_ail;
                                    double n_int = n_ail + (m_ail-n_ail) * (1 - E_int/E_ail);
                                    double p_int = p_ail * E_int/E_ail;

                                    double f_int = sqrt(pow(p_int + n_int / 2.0 - m_int / 2.0, 2) + pow(E_int, 2));

                                    // Calculations of values "Portance", "MS", and "Couple" for the entire rocket
                                    double Cn = calculateCn(Type_masquage_bi, Q_ail, E_ail, D_ref, D_ail, X_ail, m_ail, n_ail, Q_can, E_can, D_can, X_can, m_can, n_can, D_og, D2j, D1j, D2r, D1r, Q_int, E_int, D_int, f_int, m_int, n_int);
                                    double Cn0 = calculateCn0(Type_masquage_bi, Q_ail, E_ail, D_ref, D_ail, X_ail, m_ail, n_ail, Q_can, E_can, D_can, X_can, m_can, n_can, D_og, D2j, D1j, D2r, D1r);
                                    double MS_min = calculateMS_min(X_ail, XpropuRef, D_ref);
                                    double MS_max = calculateMS_max(X_can, XpropuRef, D_ref);
                                    double MS_Cn_min = calculateMS_Cn_min(MS_min, Cn);
                                    double MS_Cn_max = calculateMS_Cn_max(MS_max, Cn0);

                                    // Calculations of values "Portance", "MS", and "Couple" for the ship
                                    double CnShip = calculateCn(Type_masquage_ship, Q_ail, E_can, D_ref, D_ail, X_ail, m_can, n_can, Q_can, E_can, D_can, X_can, m_can, n_can, D_og, D2j, D1j, D2r, D1r, Q_int, E_int, D_int, f_int, m_int, n_int);
                                    double Cn0Ship = calculateCn0(Type_masquage_ship, Q_ail, E_can, D_ref, D_ail, X_ail, m_can, n_can, Q_can, E_can, D_can, X_can, m_can, n_can, D_og, D2j, D1j, D2r, D1r);
                                    double MS_minShip = calculateMS_min(X_ail, XpropuRef, D_ref);
                                    double MS_maxShip = calculateMS_max(X_can, XpropuRef, D_ref);
                                    double MS_Cn_minShip = calculateMS_Cn_min(MS_minShip, CnShip);
                                    double MS_Cn_maxShip = calculateMS_Cn_max(MS_maxShip, Cn0Ship);

                                    // Calculate the distance
                                    double distance = calculateDistance(Cn, Cn0, MS_min, MS_max, MS_Cn_min, MS_Cn_max, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);
                                    double distanceShip = calculateDistance(CnShip, Cn0Ship, MS_minShip, MS_maxShip, MS_Cn_minShip, MS_Cn_maxShip, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

                                    // Update the minimal distance and corresponding values
                                    if (distance + distanceShip < minDistanceGlobal) {
                                        minDistanceGlobal = distance + distanceShip;
                                        best_m_ail = m_ail;
                                        best_n_ail = n_ail;
                                        best_p_ail = p_ail;
                                        best_E_ail = E_ail;
                                        best_m_can = m_can;
                                        best_n_can = n_can;
                                        best_p_can = p_can;
                                        best_E_can = E_can;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Output the best values found
    printf("Best values found:\n");
    printf("m_ail: %f, n_ail: %f, p_ail: %f, E_ail: %f\n", best_m_ail, best_n_ail, best_p_ail, best_E_ail);
    printf("m_can: %f, n_can: %f, p_can: %f, E_can: %f\n", best_m_can, best_n_can, best_p_can, best_E_can);
    printf("Minimum Distance: %f\n", minDistanceGlobal);

    return 0;
}

//Calculation of values "Portance", "MS", and "Couple"
double calculateCn(const char* Type_masquage, double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail, double Q_can, double E_can, double D_can, double f_can, double m_can, double n_can, double D_og, double D2j, double D1j, double D2r, double D1r, double Q_int, double E_int, double D_int, double f_int, double m_int, double n_int) {
    // Cni calculation
    double Cni = 0;
    if (Type_masquage[0] == 'B') {
        Cni = 4 * Q_int * pow(E_int / D_ref, 2) * (1 + D_int / (2 * E_int + D_int)) / (1 + sqrt(1 + pow(2 * f_int / (m_int + n_int), 2)));
    }

    // Cnail calculation
    double Cnail = 4 * Q_ail * pow(E_ail / D_ref, 2) * (1 + D_ail / (2 * E_ail + D_ail)) / (1 + sqrt(1 + pow(2 * f_ail / (m_ail + n_ail), 2)));

    // Cnai calculation
    double Cnai = Cnail - Cni / 2;

    // Cnc calculation
    double Cnc = 0;
    if (Type_masquage[0] != 'M') {
        Cnc = 4 * Q_can * pow(E_can / D_ref, 2) * (1 + D_can / (2 * E_can + D_can)) / (1 + sqrt(1 + pow(2 * f_can / (m_can + n_can), 2)));
    }

    // Cno calculation
    double Cno = 2 * pow(D_og / D_ref, 2);

    // Cnj calculation
    double Cnj = (D2j == 0) ? 0 : 2 * (pow(D2j / D_ref, 2) - pow(D1j / D_ref, 2));

    // Cnr calculation
    double Cnr = (D2r == 0) ? 0 : 2 * (pow(D2r / D_ref, 2) - pow(D1r / D_ref, 2));

    // Calculate Cn
    return Cnai + Cnc + Cno + Cnj + Cnr;
}

double calculateCn0(const char* Type_masquage, double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail, double Q_can, double E_can, double D_can, double f_can, double m_can, double n_can, double D_og, double D2j, double D1j, double D2r, double D1r) {
    // Cnail calculation
    double Cnail = 4 * Q_ail * pow(E_ail / D_ref, 2) * (1 + D_ail / (2 * E_ail + D_ail)) / (1 + sqrt(1 + pow(2 * f_ail / (m_ail + n_ail), 2)));

    // Cnc calculation
    double Cnc = 0;
    if (Type_masquage[0] != 'M') {
        Cnc = 4 * Q_can * pow(E_can / D_ref, 2) * (1 + D_can / (2 * E_can + D_can)) / (1 + sqrt(1 + pow(2 * f_can / (m_can + n_can), 2)));
    }

    // Cno calculation
    double Cno = 2 * pow(D_og / D_ref, 2);

    // Cnj calculation
    double Cnj = (D2j == 0) ? 0 : 2 * (pow(D2j / D_ref, 2) - pow(D1j / D_ref, 2));

    // Cnr calculation
    double Cnr = (D2r == 0) ? 0 : 2 * (pow(D2r / D_ref, 2) - pow(D1r / D_ref, 2));

    // Calculate Cn0
    return Cnail + Cnc + Cno + Cnj + Cnr;
}

double calculateMS_min(double XCp, double XcgPlein, double D_ref) {
    double MS_min = 0;

    MS_min = (XCp - XcgPlein) / D_ref;

    return MS_min;
}

double calculateMS_max(double XCp0, double XcgVide, double D_ref) {
    double MS_max = 0;

    MS_max = (XCp0 - XcgVide) / D_ref;

    return MS_max;
}

double calculateMS_Cn_min(double MS_min, double Cn) {
    double MS_Cn_min = 0;

    MS_Cn_min = MS_min * Cn;

    return MS_Cn_min;
}

double calculateMS_Cn_max(double MS_max, double Cn0) {
    double MS_Cn_max = 0;

    MS_Cn_max = MS_max * Cn0;

    return MS_Cn_max;
}

//Calculation of the distance
double calculateDistance(double Cn, double Cn0, double MS_min, double MS_max, double MS_Cn_min, double MS_Cn_max, double CritCnmin, double CritCnmax, double CritMsmin, double CritMsmax, double CritMsCnmin, double CritMsCnmax) {
    double totalDistance = 0;

    totalDistance += calculateLineDistance(Cn, Cn0, CritCnmin, CritCnmax, 50);
    totalDistance += calculateLineDistance(MS_min, MS_max, CritMsmin, CritMsmax, 50);
    totalDistance += calculateLineDistance(MS_Cn_min, MS_Cn_max, CritMsCnmin, CritMsCnmax, 50);

    return totalDistance;
}

double calculateLineDistance(double value1, double value2, double center1, double center2, double centerRangeFactor) {
    double distance = 0;

    if (valueInRange(value1, center1, center2) && valueInRange(value2, center1, center2)) {
        distance = distanceToCenter(value1, (center1 + center2) / 2) + distanceToCenter(value2, (center1 + center2) / 2);
    } else {
        distance = centerRangeFactor * fabs(value1 - (center1 + center2) / 2) + centerRangeFactor * fabs(value2 - (center1 + center2) / 2);
    }

    return distance;
}

double distanceToCenter(double value, double center) {
    return fabs(value - center);
}

int valueInRange(double value, double rangeStart, double rangeEnd) {
    return (value >= rangeStart && value <= rangeEnd);
}


// Function definitions for the calculations of the critical values
double calculateCritCnmin(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 15.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0 || strncmp(Type_fusee, "Micro", 5) == 0) {
        return 15.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 15.0;
    }
    return 15.0; // Default case, if none of the above conditions are met
}

double calculateCritCnmax(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 40.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0) {
        return 30.0;
    } else if (strncmp(Type_fusee, "Micro", 5) == 0) {
        return 30.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 30.0;
    }
    return 30.0; // Default case
}

double calculateCritMsmin(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 2.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0) {
        return 1.5;
    } else if (strncmp(Type_fusee, "Micro", 5) == 0) {
        return 1.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 1.0;
    }
    return 1.0; // Default case
}

double calculateCritMsmax(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 6.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0) {
        return 6.0;
    } else if (strncmp(Type_fusee, "Micro", 5) == 0) {
        return 3.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 3.0;
    }
    return 3.0; // Default case
}

double calculateCritMsCnmin(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 40.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0) {
        return 30.0;
    } else if (strncmp(Type_fusee, "Micro", 5) == 0) {
        return 15.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 15.0;
    }
    return 15.0; // Default case
}

double calculateCritMsCnmax(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 100.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0) {
        return 100.0;
    } else if (strncmp(Type_fusee, "Micro", 5) == 0) {
        return 100.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 90.0;
    }
    return 90.0; // Default case
}

