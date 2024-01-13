#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int Long_propu;
    double MpropuVide;
    double MpropuPlein;
    int XpropuVide;
    int XpropuPlein;
} PropellerChar; // Propeller characteristics. For the moment, there's only data for Pro54-5G C, Pro75-3G C and Pro24-6G BS propellers. You can add data for other propellers if you want.

//Propeller characteristics prototype
PropellerChar propellerCaracteristics(char* Propu);

// Function Prototypes
int calculateDistance(int Cn, int Cn0, int MS_min, int MS_max, int MS_Cn_min, int MS_Cn_max, int CritCnmin, int CritCnmax, int CritMsmin, int CritMsmax, int CritMsCnmin, int CritMsCnmax);
int calculateLineDistance(int value1, int value2, int center1, int center2, int centerRangeFactor);
int distanceToCenter(int value, int center);
int valueInRange(int value, int rangeStart, int rangeEnd);

// Function Prototypes for the calculations of values "Portance", "MS", and "Couple"
int calculateCni(const char* Type_masquage, int Q_int, int E_int, int D_ref, int D_int, int f_int, int m_int, int n_int);
int calculateCnail(int Q_ail, int E_ail, int D_ref, int D_ail, int f_ail, int m_ail, int n_ail);
int calculateCnai(int Cnail, int Cni);
int calculateCnc(const char* Type_masquage, int Q_can, int E_can, int D_ref, int D_can, int f_can, int m_can, int n_can);
int calculateCno(int D_og, int D_ref);
int calculateCnj(int D2j, int D1j, int D_ref);
int calculateCnr(int D2r, int D1r, int D_ref);

int calculateCn(int Cni, int Cnail, int Cnai, int Cnc, int Cno, int Cnj, int Cnr);
int calculateCn0(int Cnail, int Cnc, int Cno, int Cnj, int Cnr);
int calculateMS_min(int XCp, int XcgPlein, int D_ref);
int calculateMS_max(int XCp0, int XcgVide, int D_ref);
int calculateMS_Cn_min(int MS_min, int Cn);
int calculateMS_Cn_max(int MS_max, int Cn0);

// Function Prototypes for the calculations of Xcp etc
int calculateXCpi(const char* Type_masquage, int X_int, int m_int, int p_int, int n_int);
int calculateXCpai(int XCpa, int Cnail, int XCpi, int Cni);
int calculateXCpc(const char* Type_masquage, int X_can, int m_can, int p_can, int n_can);
int calculateXCpj(int X_j, int l_j, int D1j, int D2j);
int calculateXCpr(int X_r, int l_r, int D1r, int D2r);
int calculateXCpo(const char* Forme_ogive, int Long_ogive);
int calculateXCpa(int X_ail, int m_ail, int p_ail, int n_ail);
int calculateXcgSans(const char* D12, int C12, int MasseVide, int XpropuRef, int Long_propu, int XpropuVide, int MasseSans, int MpropuPlein, int MpropuVide, int XpropuPlein);
int calculateMasseSans(const char* D11, int C11, int MpropuVide, int MpropuPlein);

int calculateXCp(int Cnai, int XCpai, int Cnc, int XCpc, int Cnj, int XCpj, int Cnr, int XCpr, int Cno, int XCpo);
int calculateXCp0(int Cnail, int XCpa, int Cnc, int XCpc, int Cnj, int XCpj, int Cnr, int XCpr, int Cno, int XCpo);
int calculateXcgPlein(int XcgSans, int MasseSans, int XpropuRef, int Long_propu, int XpropuPlein, int MpropuPlein);
int calculateXcgVide(int XcgSans, int MasseSans, int XpropuRef, int Long_propu, int XpropuVide, int MpropuVide);

// Function Prototypes for the calculations of the critical values for "Portance", "MS", and "Couple"
int calculateCritCnmin(const char* Type_fusee);
int calculateCritCnmax(const char* Type_fusee);
int calculateCritMsmin(const char* Type_fusee);
int calculateCritMsmax(const char* Type_fusee);
int calculateCritMsCnmin(const char* Type_fusee);
int calculateCritMsCnmax(const char* Type_fusee);

// This code only works for two-stage rockets with a bi-empennage
const char* Type_fusee = "Fusée expérimentale.";  // Rocket type
const char* Type_masquage_bi = "bi-empennage"; // Masking type for the entire rocket
const char* Type_masquage_ship = "Mono-empennage"; // Masking type for the ship
const char* withOrWithoutPropu = "sans propu"; // With or without propellant for the entire booster (corresponds to D11 and D12)
const char* withOrWithoutPropuShip = "sans propu"; // With or without propellant for the ship (corresponds to D11 and D12)

// This code takes 19 arguments in this order: 
int main(int argc, char *argv[]) {
    // Input values
    if (argc < 20) {
        fprintf(stderr, "Not enough arguments.\n");
        return 1;
    }

    int currentStep = atoi(argv[1]);

    int Long_tot = atoi(argv[2]); // Length of the booster + ship
    int D_ref = atoi(argv[3]); // Diameter of the booster + ship
    int Masse = atoi(argv[4]); // Mass of the booster + ship
    int CG = atoi(argv[5]); // Center of gravity of the booster + ship
    char* Forme_ogive = argv[6]; // Shape of the ogive
    int Long_ogive = atoi(argv[7]); // Length of the ogive
    int D_og = atoi(argv[8]); // Diameter of the ogive
    char* Propu = argv[9]; // Propeller used for the booster
    int XpropuRef = atoi(argv[10]); // Position of the propeller for the booster
    int Q_ail = atoi(argv[11]); // Number of fins for the booster
    int X_ail = atoi(argv[12]); // Position of the fins for the booster
    
    int Long_totShip = atoi(argv[13]); // Length of the ship
    int MasseShip = atoi(argv[14]); // Mass of the ship
    int CGShip = atoi(argv[15]); // Center of gravity of the ship
    char* PropuShip = argv[16]; // Propeller used for the ship
    int XpropuRefShip = atoi(argv[17]); // Position of the propeller for the ship
    int Q_can = atoi(argv[18]); // Number of fins for the ship
    int X_can = atoi(argv[19]); // Position of the fins for the ship

    // To VERIFY on the Excel file
    int X_j = 300;
    int X_r = 500;
    int l_j = 50;
    int l_r = 50;
    int D2j = 80;

    // Formulas
    int D_ail = D_ref;
    int D_can = D_ref;
    int D1j = D_og;
    int D2r = D_og;
    int D1r = D2j;
    int X_int = X_ail;
    int X_intShip = X_can;

    int XpropuRef = Long_tot;
    int XpropuRefShip = Long_totShip;

    // Variable declarations that will stock the best values found
    int best_m_ail = 0, best_n_ail = 0, best_p_ail = 0, best_E_ail = 0;
    int best_m_can = 0, best_n_can = 0, best_p_can = 0, best_E_can = 0;
    int minDistanceGlobal = DBL_MAX;
    int Q_int = (Q_ail == Q_can) ? Q_ail : 0;
    int D_int = D_ail;

    PropellerChar propeller = propellerCaracteristics(Propu);
    PropellerChar propellerShip = propellerCaracteristics(PropuShip);

    // Critical values
    int CritCnmin = calculateCritCnmin(Type_fusee);
    int CritCnmax = calculateCritCnmax(Type_fusee);
    int CritMsmin = calculateCritMsmin(Type_fusee);
    int CritMsmax = calculateCritMsmax(Type_fusee);
    int CritMsCnmin = calculateCritMsCnmin(Type_fusee);
    int CritMsCnmax = calculateCritMsCnmax(Type_fusee);

    // Iterate over the range of values
    for (int E_ail = 100; E_ail <= 400; E_ail += currentStep) {
        for (int m_ail = 50; m_ail <= 300; m_ail += currentStep) {
            for (int n_ail = 50; n_ail <= m_ail; n_ail += currentStep) {
                for (int p_ail = 0; p_ail <= n_ail; p_ail += currentStep) {
                    for (int E_can = 100; E_can <= 400; E_can += currentStep) {
                        for (int m_can = 50; m_can <= 300; m_can += currentStep) {
                            for (int n_can = 50; n_can <= m_can; n_can += currentStep) {
                                for (int p_can = 0; p_can <= n_can; p_can += currentStep) {
                                    
                                    // Calculations of values to calculate the values to calculate "Portance", "MS", and "Couple"
                                    int E_int;
                                    if ((D_can / 2.0 + E_can) <= (D_ail / 2.0)) {
                                        E_int = 0;
                                    } else if ((D_can / 2.0 + E_can) >= (D_ail / 2.0 + E_ail)) {
                                        E_int = E_ail;
                                    } else {
                                        E_int = (D_can / 2.0 + E_can) - (D_ail / 2.0);
                                    }

                                    int m_int = m_ail;
                                    int n_int = n_ail + (m_ail-n_ail) * (1 - E_int/E_ail);
                                    int p_int = p_ail * E_int/E_ail;

                                    int f_int = sqrt(pow(p_int + n_int / 2.0 - m_int / 2.0, 2) + pow(E_int, 2));

                                    // Calculations of values to calculate "Portance", "MS", and "Couple" for the entire rocket
                                    int Cni = calculateCni(Type_masquage_bi, Q_int, E_int, D_ref, D_int, f_int, m_int, n_int);
                                    int Cnail = calculateCnail(Q_ail, E_ail, D_ref, D_ail, f_int, m_ail, n_ail);
                                    int Cnai = calculateCnai(Cnail, Cni);
                                    int Cnc = calculateCnc(Type_masquage_bi, Q_can, E_can, D_ref, D_can, f_int, m_can, n_can);
                                    int Cno = calculateCno(D_og, D_ref);
                                    int Cnj = calculateCnj(D2j, D1j, D_ref);
                                    int Cnr = calculateCnr(D2r, D1r, D_ref);

                                    // Calculations of values to calculate "Portance", "MS", and "Couple" for the ship
                                    int CniShip = calculateCni(Type_masquage_ship, Q_int, E_int, D_ref, D_int, f_int, m_int, n_int);
                                    int CnailShip = calculateCnail(Q_can, E_can, D_ref, D_can, f_int, m_can, n_can);
                                    int CnaiShip = calculateCnai(CnailShip, CniShip);
                                    int CncShip = calculateCnc(Type_masquage_ship, Q_can, E_can, D_ref, D_can, f_int, m_can, n_can);
                                    int CnoShip = calculateCno(D_og, D_ref);
                                    int CnjShip = calculateCnj(D2j, D1j, D_ref);
                                    int CnrShip = calculateCnr(D2r, D1r, D_ref);

                                    // Calculations of Xcp for the entire rocket etc
                                    int XCpi = calculateXCpi(Type_masquage_bi, X_int, m_int, p_int, n_int);
                                    int XCpa = calculateXCpa(X_ail, m_ail, p_ail, n_ail);
                                    int XCpai = calculateXCpai(XCpa, Cnail, XCpi, Cni);
                                    int XCpc = calculateXCpc(Type_masquage_bi, X_can, m_can, p_can, n_can);
                                    int XCpj = calculateXCpj(X_j, l_j, D1j, D2j);
                                    int XCpr = calculateXCpr(X_r, l_r, D1r, D2r);
                                    int XCpo = calculateXCpo(Forme_ogive, Long_ogive);
                                    int MasseSans = calculateMasseSans(withOrWithoutPropu, Masse, propeller.MpropuVide, propeller.MpropuPlein);
                                    int MasseVide = MasseSans + propeller.MpropuVide;
                                    int XcgSans = calculateXcgSans(withOrWithoutPropu, CG, MasseVide, XpropuRef, propeller.Long_propu, propeller.XpropuVide, MasseSans, propeller.MpropuPlein, propeller.MpropuVide, propeller.XpropuPlein);

                                    int Xcp = calculateXCp(Cnai, XCpai, Cnc, XCpc, Cnj, XCpj, Cnr, XCpr, Cno, XCpo);
                                    int Xcp0 = calculateXCp0(Cnail, XCpa, Cnc, XCpc, Cnj, XCpj, Cnr, XCpr, Cno, XCpo);
                                    int XcgPlein = calculateXcgPlein(XcgSans, MasseSans, XpropuRef, propeller.Long_propu, propeller.XpropuPlein, propeller.MpropuPlein);
                                    int XcgVide = calculateXcgVide(XcgSans, MasseSans, XpropuRef, propeller.Long_propu, propeller.XpropuVide, propeller.MpropuVide);

                                    // Calculations of Xcp for the ship etc
                                    int XCpiShip = calculateXCpi(Type_masquage_ship, X_intShip, m_int, p_int, n_int);
                                    int XCpaShip = calculateXCpa(X_can, m_can, p_can, n_can);
                                    int XCpaiShip = calculateXCpai(XCpaShip, CnailShip, XCpiShip, CniShip);
                                    int XCpcShip = calculateXCpc(Type_masquage_ship, X_can, m_can, p_can, n_can);
                                    int XCpjShip = calculateXCpj(X_j, l_j, D1j, D2j);
                                    int XCprShip = calculateXCpr(X_r, l_r, D1r, D2r);
                                    int XCpoShip = calculateXCpo(Forme_ogive, Long_ogive);
                                    int MasseSansShip = calculateMasseSans(withOrWithoutPropuShip, MasseShip, propellerShip.MpropuVide, propellerShip.MpropuPlein);
                                    int MasseVideShip = MasseSansShip + propellerShip.MpropuVide;
                                    int XcgSansShip = calculateXcgSans(withOrWithoutPropuShip, CGShip, MasseVideShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuVide, MasseSansShip, propellerShip.MpropuPlein, propellerShip.MpropuVide, propellerShip.XpropuPlein);

                                    int XcpShip = calculateXCp(CnaiShip, XCpaiShip, CncShip, XCpcShip, CnjShip, XCpjShip, CnrShip, XCprShip, CnoShip, XCpoShip);
                                    int Xcp0Ship = calculateXCp0(CnailShip, XCpaShip, CncShip, XCpcShip, CnjShip, XCpjShip, CnrShip, XCprShip, CnoShip, XCpoShip);
                                    int XcgPleinShip = calculateXcgPlein(XcgSansShip, MasseSansShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuPlein, propellerShip.MpropuPlein);
                                    int XcgVideShip = calculateXcgVide(XcgSansShip, MasseSansShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuVide, propellerShip.MpropuVide);

                                    // Calculations of values "Portance", "MS", and "Couple" for the entire rocket
                                    int Cn = calculateCn(Cni, Cnail, Cnai, Cnc, Cno, Cnj, Cnr);
                                    int Cn0 = calculateCn0(Cnail, Cnc, Cno, Cnj, Cnr);
                                    int MS_min = calculateMS_min(Xcp, XcgPlein, D_ref);
                                    int MS_max = calculateMS_max(Xcp0, XcgVide, D_ref);
                                    int MS_Cn_min = calculateMS_Cn_min(MS_min, Cn);
                                    int MS_Cn_max = calculateMS_Cn_max(MS_max, Cn0);

                                    // Calculations of values "Portance", "MS", and "Couple" for the ship
                                    int CnShip = calculateCn(CniShip, CnailShip, CnaiShip, CncShip, CnoShip, CnjShip, CnrShip);
                                    int Cn0Ship = calculateCn0(CnailShip, CncShip, CnoShip, CnjShip, CnrShip);
                                    int MS_minShip = calculateMS_min(XcpShip, XcgPleinShip, D_ref);
                                    int MS_maxShip = calculateMS_max(Xcp0Ship, XcgVideShip, D_ref);
                                    int MS_Cn_minShip = calculateMS_Cn_min(MS_minShip, CnShip);
                                    int MS_Cn_maxShip = calculateMS_Cn_max(MS_maxShip, Cn0Ship);

                                    // Calculate the distance
                                    int distance = calculateDistance(Cn, Cn0, MS_min, MS_max, MS_Cn_min, MS_Cn_max, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);
                                    int distanceShip = calculateDistance(CnShip, Cn0Ship, MS_minShip, MS_maxShip, MS_Cn_minShip, MS_Cn_maxShip, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

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

// Function to get propeller characteristics based on its name
PropellerChar propellerCaracteristics(char* Propu) {
    PropellerChar propeller;

    // Original (Pro75-3G C)
    if (strcmp(Propu, "Pro75-3G C") == 0) {
        propeller.Long_propu = 486;
        propeller.MpropuVide = 1.638;
        propeller.MpropuPlein = 3.511;
        propeller.XpropuVide = 243;
        propeller.XpropuPlein = 243;
    }
    // Pandora (Pro24-6G BS)
    else if (strcmp(Propu, "Pro24-6G BS") == 0) {
        propeller.Long_propu = 228;
        propeller.MpropuVide = 0.0843;
        propeller.MpropuPlein = 0.1599;
        propeller.XpropuVide = 114;
        propeller.XpropuPlein = 114;
    }
    // Barasinga (Pro54-5G C)
    else if (strcmp(Propu, "Pro54-5G C") == 0) {
        propeller.Long_propu = 488;
        propeller.MpropuVide = 0.652;
        propeller.MpropuPlein = 1.685;
        propeller.XpropuVide = 240;
        propeller.XpropuPlein = 250;
    }
    // Default case to handle unknown propellers
    else {
        printf("Unknown propeller type.\n");
        propeller.Long_propu = -1;
        propeller.MpropuVide = -1;
        propeller.MpropuPlein = -1;
        propeller.XpropuVide = -1;
        propeller.XpropuPlein = -1;
    }

    return propeller;
}

//Calculation of values "Portance", "MS", and "Couple"
int calculateCni(const char* Type_masquage, int Q_int, int E_int, int D_ref, int D_int, int f_int, int m_int, int n_int) {
    if (Type_masquage[0] == 'B') {
        return 4 * Q_int * pow(E_int / D_ref, 2) * (1 + D_int / (2 * E_int + D_int)) / (1 + sqrt(1 + pow(2 * f_int / (m_int + n_int), 2)));
    }
    return 0;
}

int calculateCnail(int Q_ail, int E_ail, int D_ref, int D_ail, int f_ail, int m_ail, int n_ail) {
    return 4 * Q_ail * pow(E_ail / D_ref, 2) * (1 + D_ail / (2 * E_ail + D_ail)) / (1 + sqrt(1 + pow(2 * f_ail / (m_ail + n_ail), 2)));
}

int calculateCnai(int Cnail, int Cni) {
    return Cnail - Cni / 2;
}

int calculateCnc(const char* Type_masquage, int Q_can, int E_can, int D_ref, int D_can, int f_can, int m_can, int n_can) {
    if (Type_masquage[0] != 'M') {
        return 4 * Q_can * pow(E_can / D_ref, 2) * (1 + D_can / (2 * E_can + D_can)) / (1 + sqrt(1 + pow(2 * f_can / (m_can + n_can), 2)));
    }
    return 0;
}

int calculateCno(int D_og, int D_ref) {
    return 2 * pow(D_og / D_ref, 2);
}

int calculateCnj(int D2j, int D1j, int D_ref) {
    int Cnj = (D2j == 0) ? 0 : 2 * (pow(D2j / D_ref, 2) - pow(D1j / D_ref, 2));
    return Cnj;
}

int calculateCnr(int D2r, int D1r, int D_ref) {
    int Cnr = (D2r == 0) ? 0 : 2 * (pow(D2r / D_ref, 2) - pow(D1r / D_ref, 2));
    return Cnr;
}


int calculateCn(int Cni, int Cnail, int Cnai, int Cnc, int Cno, int Cnj, int Cnr) {
    return Cnai + Cnc + Cno + Cnj + Cnr;
}

int calculateCn0(int Cnail, int Cnc, int Cno, int Cnj, int Cnr) {
    return Cnail + Cnc + Cno + Cnj + Cnr;
}

int calculateMS_min(int XCp, int XcgPlein, int D_ref) {
    int MS_min = 0;

    MS_min = (XCp - XcgPlein) / D_ref;

    return MS_min;
}

int calculateMS_max(int XCp0, int XcgVide, int D_ref) {
    int MS_max = 0;

    MS_max = (XCp0 - XcgVide) / D_ref;

    return MS_max;
}

int calculateMS_Cn_min(int MS_min, int Cn) {
    int MS_Cn_min = 0;

    MS_Cn_min = MS_min * Cn;

    return MS_Cn_min;
}

int calculateMS_Cn_max(int MS_max, int Cn0) {
    int MS_Cn_max = 0;

    MS_Cn_max = MS_max * Cn0;

    return MS_Cn_max;
}

//Calculation of Xcp etc
int calculateXCpi(const char* Type_masquage, int X_int, int m_int, int p_int, int n_int) {
    if (Type_masquage[0] == 'B') {
        return X_int - m_int + p_int * (m_int + 2 * n_int) / (3 * (m_int + n_int)) + (m_int + n_int - m_int * n_int / (m_int + n_int)) / 6;
    }
    return 0;
}

int calculateXCpai(int XCpa, int Cnail, int XCpi, int Cni) {
    return (XCpa * Cnail - 0.5 * XCpi * Cni) / (Cnail - Cni / 2);
}

int calculateXCpc(const char* Type_masquage, int X_can, int m_can, int p_can, int n_can) {
    if (Type_masquage[0] != 'M') {
        return X_can - m_can + p_can * (m_can + 2 * n_can) / (3 * (m_can + n_can)) + (m_can + n_can - m_can * n_can / (m_can + n_can)) / 6;
    }
    return 0;
}

int calculateXCpj(int X_j, int l_j, int D1j, int D2j) {
    if (D2j != 0) {
        return X_j + l_j / 3 * (1 + 1 / (1 + D1j / D2j));
    }
    return 0;
}

int calculateXCpr(int X_r, int l_r, int D1r, int D2r) {
    if (D2r != 0) {
        return X_r + l_r / 3 * (1 + 1 / (1 + D1r / D2r));
    }
    return 0;
}

int calculateXCpo(const char* Forme_ogive, int Long_ogive) {
    if (strncmp(Forme_ogive, "Parab", 5) == 0) {
        return 0.5 * Long_ogive;
    } else if (strncmp(Forme_ogive, "Ogiv", 4) == 0) {
        return 7.0 / 15.0 * Long_ogive;
    } else if (strncmp(Forme_ogive, "Con", 3) == 0) {
        return 2.0 / 3.0 * Long_ogive;
    }
    return 0;
}

int calculateXCpa(int X_ail, int m_ail, int p_ail, int n_ail) {
    return X_ail - m_ail + p_ail * (m_ail + 2 * n_ail) / (3 * (m_ail + n_ail)) + (m_ail + n_ail - m_ail * n_ail / (m_ail + n_ail)) / 6;
}

int calculateXcgSans(const char* D12, int C12, int MasseVide, int XpropuRef, int Long_propu, int XpropuVide, int MasseSans, int MpropuPlein, int MpropuVide, int XpropuPlein) {
    if (strcmp(D12, "sans propu") == 0 || strcmp(D12, "without motor") == 0) {
        return C12;
    } else if (strcmp(D12, "avec propu vide") == 0 || strcmp(D12, "with empty motor") == 0) {
        return (C12 * MasseVide - (XpropuRef - Long_propu + XpropuVide) * MpropuVide) / MasseSans;
    } else if (strcmp(D12, "avec propu plein") == 0 || strcmp(D12, "with loaded motor") == 0) {
        return (C12 * MasseSans - (XpropuRef - Long_propu + XpropuPlein) * MpropuPlein) / MasseSans;
    }
    return -1; // Error case
}

int calculateMasseSans(const char* D11, int C11, int MpropuVide, int MpropuPlein) {
    if (strcmp(D11, "sans propu") == 0 || strcmp(D11, "without motor") == 0) {
        return C11 / 1000;
    } else if (strcmp(D11, "avec propu vide") == 0 || strcmp(D11, "with empty motor") == 0) {
        return C11 / 1000 - MpropuVide;
    } else if (strcmp(D11, "avec propu plein") == 0 || strcmp(D11, "with loaded motor") == 0) {
        return C11 / 1000 - MpropuPlein;
    }
    return -1; // Error case
}

int calculateXCp(int Cnai, int XCpai, int Cnc, int XCpc, int Cnj, int XCpj, int Cnr, int XCpr, int Cno, int XCpo) {
    int sumCnX = Cnai * XCpai + Cnc * XCpc + Cnj * XCpj + Cnr * XCpr + Cno * XCpo;
    int sumCn = Cnai + Cnc + Cnj + Cnr + Cno;
    return sumCn == 0 ? 0 : sumCnX / sumCn;
}

int calculateXCp0(int Cnail, int XCpa, int Cnc, int XCpc, int Cnj, int XCpj, int Cnr, int XCpr, int Cno, int XCpo) {
    int sumCnX = Cnail * XCpa + Cnc * XCpc + Cnj * XCpj + Cnr * XCpr + Cno * XCpo;
    int sumCn = Cnail + Cnc + Cnj + Cnr + Cno;
    return sumCn == 0 ? 0 : sumCnX / sumCn;
}

int calculateXcgPlein(int XcgSans, int MasseSans, int XpropuRef, int Long_propu, int XpropuPlein, int MpropuPlein) {
    return (XcgSans * MasseSans + (XpropuRef - Long_propu + XpropuPlein) * MpropuPlein) / (MasseSans + MpropuPlein);
}

int calculateXcgVide(int XcgSans, int MasseSans, int XpropuRef, int Long_propu, int XpropuVide, int MpropuVide) {
    return (XcgSans * MasseSans + (XpropuRef - Long_propu + XpropuVide) * MpropuVide) / (MasseSans + MpropuVide);
}

//Calculation of the distance
int calculateDistance(int Cn, int Cn0, int MS_min, int MS_max, int MS_Cn_min, int MS_Cn_max, int CritCnmin, int CritCnmax, int CritMsmin, int CritMsmax, int CritMsCnmin, int CritMsCnmax) {
    int totalDistance = 0;

    totalDistance += calculateLineDistance(Cn, Cn0, CritCnmin, CritCnmax, 50);
    totalDistance += calculateLineDistance(MS_min, MS_max, CritMsmin, CritMsmax, 50);
    totalDistance += calculateLineDistance(MS_Cn_min, MS_Cn_max, CritMsCnmin, CritMsCnmax, 50);

    return totalDistance;
}

int calculateLineDistance(int value1, int value2, int center1, int center2, int centerRangeFactor) {
    int distance = 0;

    if (valueInRange(value1, center1, center2) && valueInRange(value2, center1, center2)) {
        distance = distanceToCenter(value1, (center1 + center2) / 2) + distanceToCenter(value2, (center1 + center2) / 2);
    } else {
        distance = centerRangeFactor * fabs(value1 - (center1 + center2) / 2) + centerRangeFactor * fabs(value2 - (center1 + center2) / 2);
    }

    return distance;
}

int distanceToCenter(int value, int center) {
    return fabs(value - center);
}

int valueInRange(int value, int rangeStart, int rangeEnd) {
    return (value >= rangeStart && value <= rangeEnd);
}


// Function definitions for the calculations of the critical values
int calculateCritCnmin(const char* Type_fusee) {
    if (Type_fusee[strlen(Type_fusee) - 1] == '.') {
        return 15.0;
    } else if (Type_fusee[0] == 'R' || Type_fusee[0] == ',' || strncmp(Type_fusee, "Mini", 4) == 0 || strncmp(Type_fusee, "Micro", 5) == 0) {
        return 15.0;
    } else if (Type_fusee[strlen(Type_fusee) - 1] == ' ') {
        return 15.0;
    }
    return 15.0; // Default case, if none of the above conditions are met
}

int calculateCritCnmax(const char* Type_fusee) {
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

int calculateCritMsmin(const char* Type_fusee) {
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

int calculateCritMsmax(const char* Type_fusee) {
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

int calculateCritMsCnmin(const char* Type_fusee) {
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

int calculateCritMsCnmax(const char* Type_fusee) {
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

