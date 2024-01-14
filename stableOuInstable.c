#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    double Long_propu;
    double MpropuVide;
    double MpropuPlein;
    double XpropuVide;
    double XpropuPlein;
} PropellerChar; // Propeller characteristics. For the moment, there's only data for Pro54-5G C, Pro75-3G C and Pro24-6G BS propellers. You can add data for other propellers if you want.

//Propeller characteristics prototype
PropellerChar propellerCaracteristics(char* Propu);

// Function Prototypes
double calculateDistance(double Cn, double Cn0, double MS_min, double MS_max, double MS_Cn_min, double MS_Cn_max, double CritCnmin, double CritCnmax, double CritMsmin, double CritMsmax, double CritMsCnmin, double CritMsCnmax);
double calculateLineDistance(double value1, double value2, double center1, double center2, double centerRangeFactor);
double distanceToCenter(double value, double center);
double valueInRange(double value, double rangeStart, double rangeEnd);

// Function Prototypes for the calculations of values "Portance", "MS", and "Couple"
double calculateCni(const char* Type_masquage, double Q_int, double E_int, double D_ref, double D_int, double f_int, double m_int, double n_int);
double calculateCnail(double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail);
double calculateCnai(double Cnail, double Cni);
double calculateCnc(const char* Type_masquage, double Q_can, double E_can, double D_ref, double D_can, double f_can, double m_can, double n_can);
double calculateCno(double D_og, double D_ref);
double calculateCnj(double D2j, double D1j, double D_ref);
double calculateCnr(double D2r, double D1r, double D_ref);

double calculateCn(double Cni, double Cnail, double Cnai, double Cnc, double Cno, double Cnj, double Cnr);
double calculateCn0(double Cnail, double Cnc, double Cno, double Cnj, double Cnr);
double calculateMS_min(double XCp, double XcgPlein, double D_ref);
double calculateMS_max(double XCp0, double XcgVide, double D_ref);
double calculateMS_Cn_min(double MS_min, double Cn);
double calculateMS_Cn_max(double MS_max, double Cn0);

// Function Prototypes for the calculations of Xcp etc
double calculateXCpi(const char* Type_masquage, double X_int, double m_int, double p_int, double n_int);
double calculateXCpai(double XCpa, double Cnail, double XCpi, double Cni);
double calculateXCpc(const char* Type_masquage, double X_can, double m_can, double p_can, double n_can);
double calculateXCpj(double X_j, double l_j, double D1j, double D2j);
double calculateXCpr(double X_r, double l_r, double D1r, double D2r);
double calculateXCpo(const char* Forme_ogive, double Long_ogive);
double calculateXCpa(double X_ail, double m_ail, double p_ail, double n_ail);
double calculateXcgSans(const char* D12, double C12, double MasseVide, double XpropuRef, double Long_propu, double XpropuVide, double MasseSans, double MpropuPlein, double MpropuVide, double XpropuPlein);
double calculateMasseSans(const char* D11, double C11, double MpropuVide, double MpropuPlein);

double calculateXCp(double Cnai, double XCpai, double Cnc, double XCpc, double Cnj, double XCpj, double Cnr, double XCpr, double Cno, double XCpo);
double calculateXCp0(double Cnail, double XCpa, double Cnc, double XCpc, double Cnj, double XCpj, double Cnr, double XCpr, double Cno, double XCpo);
double calculateXcgPlein(double XcgSans, double MasseSans, double XpropuRef, double Long_propu, double XpropuPlein, double MpropuPlein);
double calculateXcgVide(double XcgSans, double MasseSans, double XpropuRef, double Long_propu, double XpropuVide, double MpropuVide);

// Function Prototypes for the calculations of the critical values for "Portance", "MS", and "Couple"
double calculateCritCnmin(const char* Type_fusee);
double calculateCritCnmax(const char* Type_fusee);
double calculateCritMsmin(const char* Type_fusee);
double calculateCritMsmax(const char* Type_fusee);
double calculateCritMsCnmin(const char* Type_fusee);
double calculateCritMsCnmax(const char* Type_fusee);

// Stable or not?
char* checkStability(double Cn, double Cn0, double MS_min, double MS_max, double MS_Cn_min, double MS_Cn_max, double CritCnmin, double CritCnmax, double CritMsmin, double CritMsmax, double CritMsCnmin, double CritMsCnmax);

// ------------------------------------------------------------------------------------------------------------------------------------------------------

// This code only works for two-stage rockets with a bi-empennage
const char* Type_fusee = "Fusée expérimentale.";  // Rocket type
const char* Type_masquage_bi = "Bi-empennage"; // Masking type for the entire rocket
const char* Type_masquage_ship = "Mono-empennage"; // Masking type for the ship
const char* withOrWithoutPropu = "sans propu"; // With or without propellant for the entire booster (corresponds to D11 and D12)
const char* withOrWithoutPropuShip = "sans propu"; // With or without propellant for the ship (corresponds to D11 and D12)

// This code takes 19 arguments in this order: 
int main(int argc, char *argv[]) {
    // Input values
    if (argc < 19) {
        fprintf(stderr, "Not enough arguments.\n");
        return 1;
    }

    int currentStep = atoi(argv[1]);

    double Long_tot = atof(argv[2]); // Length of the booster + ship
    double D_ref = atof(argv[3]); // Diameter of the booster + ship
    double Masse = atof(argv[4]); // Mass of the booster + ship
    double CG = atof(argv[5]); // Center of gravity of the booster + ship
    char* Forme_ogive = argv[6]; // Shape of the ogive
    double Long_ogive = atof(argv[7]); // Length of the ogive
    double D_og = atof(argv[8]); // Diameter of the ogive
    char* Propu = argv[9]; // Propeller used for the booster
    double Q_ail = atof(argv[10]); // Number of fins for the booster
    double X_ail = atof(argv[11]); // Position of the fins for the booster
    
    double Long_totShip = atof(argv[12]); // Length of the ship
    double MasseShip = atof(argv[13]); // Mass of the ship
    double CGShip = atof(argv[14]); // Center of gravity of the ship
    char* PropuShip = argv[15]; // Propeller used for the ship
    double Q_can = atof(argv[16]); // Number of fins for the ship
    double X_can = atof(argv[17]); // Position of the fins for the ship

    double alpha = atof(argv[18]); // Importance of the stability of the entire rocket over the ship alone

    // To VERIFY on the Excel file
    double X_j = 300;
    double X_r = 500;
    double l_j = 50;
    double l_r = 50;
    double D2j = 80;

    // Formulas
    double D_ail = D_ref;
    double D_can = D_ref;
    double D1j = D_og;
    double D2r = D_og;
    double D1r = D2j;
    double X_int = X_ail;
    double X_intShip = X_can;

    double XpropuRef = Long_tot;
    double XpropuRefShip = Long_totShip;

    // Variable declarations that will stock the best values found
    int best_m_ail = 0, best_n_ail = 0, best_p_ail = 0, best_E_ail = 0;
    int best_m_can = 0, best_n_can = 0, best_p_can = 0, best_E_can = 0;
    double minDistanceGlobal = 1000000;
    char best_stability[20];
    char best_stabilityShip[20];

    double best_cn = 0, best_cn0 = 0, best_ms_min = 0, best_ms_max = 0, best_ms_cn_min = 0, best_ms_cn_max = 0;
    double best_cn_ship = 0, best_cn0_ship = 0, best_ms_min_ship = 0, best_ms_max_ship = 0, best_ms_cn_min_ship = 0, best_ms_cn_max_ship = 0;
    
    double Q_int = (Q_ail == Q_can) ? Q_ail : 0;
    double D_int = D_ail;

    // Get propeller characteristics
    PropellerChar propeller = propellerCaracteristics(Propu);
    PropellerChar propellerShip = propellerCaracteristics(PropuShip);

    // Critical values
    double CritCnmin = calculateCritCnmin(Type_fusee);
    double CritCnmax = calculateCritCnmax(Type_fusee);
    double CritMsmin = calculateCritMsmin(Type_fusee);
    double CritMsmax = calculateCritMsmax(Type_fusee);
    double CritMsCnmin = calculateCritMsCnmin(Type_fusee);
    double CritMsCnmax = calculateCritMsCnmax(Type_fusee);

    // Initialize start time
    clock_t start_time = clock();
    double elapsed_seconds = 0;

    // File opening
    FILE *pf;

    pf=fopen("valeursOptimales.out","w");
    if (!pf) {
        perror("error.out");
        pf=stderr;
    }
    else
        setvbuf(pf, NULL, _IOLBF, 0);
    
    long iterationsCompleted = 0;

    // Iterate over the range of values
    for (int E_ail = 100; E_ail <= 400; E_ail += currentStep) {
        for (int m_ail = 50; m_ail <= 300; m_ail += currentStep) {
            for (int n_ail = 50; n_ail <= m_ail; n_ail += currentStep) {
                for (int p_ail = 0; p_ail <= n_ail; p_ail += currentStep) {
                    for (int E_can = 100; E_can <= 400; E_can += currentStep) {
                        for (int m_can = 50; m_can <= 300; m_can += currentStep) {
                            for (int n_can = 50; n_can <= m_can; n_can += currentStep) {
                                for (int p_can = 0; p_can <= n_can; p_can += currentStep) {

                                    // Update the number of iterations completed
                                    iterationsCompleted++;

                                    // Calculations of values to calculate the values to calculate "Portance", "MS", and "Couple"
                                    double E_int;
                                    if ((D_can / 2.0 + E_can) <= (D_ail / 2.0)) {
                                        E_int = 0;
                                    } else if ((D_can / 2.0 + E_can) >= (D_ail / 2.0 + E_ail)) {
                                        E_int = E_ail;
                                    } else {
                                        E_int = (D_can / 2.0 + E_can) - (D_ail / 2.0);
                                    }

                                    double f_ail = sqrt(pow((p_ail + n_ail / 2.0 - m_ail / 2.0), 2) + pow(E_ail, 2));
                                    double f_can = sqrt(pow((p_can + n_can / 2.0 - m_can / 2.0), 2) + pow(E_can, 2));

                                    // Calculations for the entire rocket
                                    double m_int = m_ail;
                                    double n_int = n_ail + (m_ail-n_ail) * (1 - E_int/E_ail);
                                    double p_int = p_ail * E_int/E_ail;
                                    double f_int = sqrt(pow(p_int + n_int / 2.0 - m_int / 2.0, 2) + pow(E_int, 2));

                                    // Calculations for the ship
                                    double m_intShip = m_can;
                                    double n_intShip = n_can + (m_can-n_can) * (1 - E_int/E_can);
                                    double p_intShip = p_can * E_int/E_can;
                                    double f_intShip = sqrt(pow(p_intShip + n_intShip / 2.0 - m_intShip / 2.0, 2) + pow(E_int, 2));

                                    // Calculations of values to calculate "Portance", "MS", and "Couple" for the entire rocket
                                    double Cni = calculateCni(Type_masquage_bi, Q_int, E_int, D_ref, D_int, f_int, m_int, n_int);
                                    double Cnail = calculateCnail(Q_ail, E_ail, D_ref, D_ail, f_ail, m_ail, n_ail);
                                    double Cnai = calculateCnai(Cnail, Cni);
                                    double Cnc = calculateCnc(Type_masquage_bi, Q_can, E_can, D_ref, D_can, f_can, m_can, n_can);
                                    double Cno = calculateCno(D_og, D_ref);
                                    double Cnj = calculateCnj(D2j, D1j, D_ref);
                                    double Cnr = calculateCnr(D2r, D1r, D_ref);

                                    // Calculations of values to calculate "Portance", "MS", and "Couple" for the ship
                                    double CniShip = calculateCni(Type_masquage_ship, Q_int, E_int, D_ref, D_int, f_intShip, m_intShip, n_intShip);
                                    double CnailShip = calculateCnail(Q_can, E_can, D_ref, D_can, f_can, m_can, n_can);
                                    double CnaiShip = calculateCnai(CnailShip, CniShip);
                                    double CncShip = calculateCnc(Type_masquage_ship, Q_can, E_can, D_ref, D_can, f_can, m_can, n_can);
                                    double CnoShip = calculateCno(D_og, D_ref);
                                    double CnjShip = calculateCnj(D2j, D1j, D_ref);
                                    double CnrShip = calculateCnr(D2r, D1r, D_ref);

                                    // Calculations of Xcp for the entire rocket etc
                                    double XCpi = calculateXCpi(Type_masquage_bi, X_int, m_int, p_int, n_int);
                                    double XCpa = calculateXCpa(X_ail, m_ail, p_ail, n_ail);
                                    double XCpai = calculateXCpai(XCpa, Cnail, XCpi, Cni);
                                    double XCpc = calculateXCpc(Type_masquage_bi, X_can, m_can, p_can, n_can);
                                    double XCpj = calculateXCpj(X_j, l_j, D1j, D2j);
                                    double XCpr = calculateXCpr(X_r, l_r, D1r, D2r);
                                    double XCpo = calculateXCpo(Forme_ogive, Long_ogive);
                                    double MasseSans = calculateMasseSans(withOrWithoutPropu, Masse, propeller.MpropuVide, propeller.MpropuPlein);
                                    double MasseVide = MasseSans + propeller.MpropuVide;
                                    double XcgSans = calculateXcgSans(withOrWithoutPropu, CG, MasseVide, XpropuRef, propeller.Long_propu, propeller.XpropuVide, MasseSans, propeller.MpropuPlein, propeller.MpropuVide, propeller.XpropuPlein);

                                    double Xcp = calculateXCp(Cnai, XCpai, Cnc, XCpc, Cnj, XCpj, Cnr, XCpr, Cno, XCpo);
                                    double Xcp0 = calculateXCp0(Cnail, XCpa, Cnc, XCpc, Cnj, XCpj, Cnr, XCpr, Cno, XCpo);
                                    double XcgPlein = calculateXcgPlein(XcgSans, MasseSans, XpropuRef, propeller.Long_propu, propeller.XpropuPlein, propeller.MpropuPlein);
                                    double XcgVide = calculateXcgVide(XcgSans, MasseSans, XpropuRef, propeller.Long_propu, propeller.XpropuVide, propeller.MpropuVide);

                                    // Calculations of Xcp for the ship etc
                                    double XCpiShip = calculateXCpi(Type_masquage_ship, X_intShip, m_intShip, p_intShip, n_intShip);
                                    double XCpaShip = calculateXCpa(X_can, m_can, p_can, n_can);
                                    double XCpaiShip = calculateXCpai(XCpaShip, CnailShip, XCpiShip, CniShip);
                                    double XCpcShip = calculateXCpc(Type_masquage_ship, X_can, m_can, p_can, n_can);
                                    double XCpjShip = calculateXCpj(X_j, l_j, D1j, D2j);
                                    double XCprShip = calculateXCpr(X_r, l_r, D1r, D2r);
                                    double XCpoShip = calculateXCpo(Forme_ogive, Long_ogive);
                                    double MasseSansShip = calculateMasseSans(withOrWithoutPropuShip, MasseShip, propellerShip.MpropuVide, propellerShip.MpropuPlein);
                                    double MasseVideShip = MasseSansShip + propellerShip.MpropuVide;
                                    double XcgSansShip = calculateXcgSans(withOrWithoutPropuShip, CGShip, MasseVideShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuVide, MasseSansShip, propellerShip.MpropuPlein, propellerShip.MpropuVide, propellerShip.XpropuPlein);

                                    double XcpShip = calculateXCp(CnaiShip, XCpaiShip, CncShip, XCpcShip, CnjShip, XCpjShip, CnrShip, XCprShip, CnoShip, XCpoShip);
                                    double Xcp0Ship = calculateXCp0(CnailShip, XCpaShip, CncShip, XCpcShip, CnjShip, XCpjShip, CnrShip, XCprShip, CnoShip, XCpoShip);
                                    double XcgPleinShip = calculateXcgPlein(XcgSansShip, MasseSansShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuPlein, propellerShip.MpropuPlein);
                                    double XcgVideShip = calculateXcgVide(XcgSansShip, MasseSansShip, XpropuRefShip, propellerShip.Long_propu, propellerShip.XpropuVide, propellerShip.MpropuVide);

                                    // Calculations of values "Portance", "MS", and "Couple" for the entire rocket
                                    double Cn = calculateCn(Cni, Cnail, Cnai, Cnc, Cno, Cnj, Cnr);
                                    double Cn0 = calculateCn0(Cnail, Cnc, Cno, Cnj, Cnr);
                                    double MS_min = calculateMS_min(Xcp, XcgPlein, D_ref);
                                    double MS_max = calculateMS_max(Xcp0, XcgVide, D_ref);
                                    double MS_Cn_min = calculateMS_Cn_min(MS_min, Cn);
                                    double MS_Cn_max = calculateMS_Cn_max(MS_max, Cn0);

                                    // Calculations of values "Portance", "MS", and "Couple" for the ship
                                    double CnShip = calculateCn(CniShip, CnailShip, CnaiShip, CncShip, CnoShip, CnjShip, CnrShip);
                                    double Cn0Ship = calculateCn0(CnailShip, CncShip, CnoShip, CnjShip, CnrShip);
                                    double MS_minShip = calculateMS_min(XcpShip, XcgPleinShip, D_ref);
                                    double MS_maxShip = calculateMS_max(Xcp0Ship, XcgVideShip, D_ref);
                                    double MS_Cn_minShip = calculateMS_Cn_min(MS_minShip, CnShip);
                                    double MS_Cn_maxShip = calculateMS_Cn_max(MS_maxShip, Cn0Ship);

                                    // Calculate the distance
                                    double distance = alpha*calculateDistance(Cn, Cn0, MS_min, MS_max, MS_Cn_min, MS_Cn_max, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);
                                    double distanceShip = calculateDistance(CnShip, Cn0Ship, MS_minShip, MS_maxShip, MS_Cn_minShip, MS_Cn_maxShip, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

                                    // Check if the rocket is stable
                                    char* stability = checkStability(Cn, Cn0, MS_min, MS_max, MS_Cn_min, MS_Cn_max, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

                                    // Check if the ship is stable
                                    char* stabilityShip = checkStability(CnShip, Cn0Ship, MS_minShip, MS_maxShip, MS_Cn_minShip, MS_Cn_maxShip, CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

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

                                        best_cn = Cn;
                                        best_cn0 = Cn0;
                                        best_ms_min = MS_min;
                                        best_ms_max = MS_max;
                                        best_ms_cn_min = MS_Cn_min;
                                        best_ms_cn_max = MS_Cn_max;

                                        best_cn_ship = CnShip;
                                        best_cn0_ship = Cn0Ship;
                                        best_ms_min_ship = MS_minShip;
                                        best_ms_max_ship = MS_maxShip;
                                        best_ms_cn_min_ship = MS_Cn_minShip;
                                        best_ms_cn_max_ship = MS_Cn_maxShip;

                                        strcpy(best_stability, stability);
                                        strcpy(best_stabilityShip, stabilityShip);
                                    }

                                    // Write the best values found for the moment to the file and print them in the console every 10 000 000 iterations
                                    if (iterationsCompleted % 10000000 == 0) {
                                        fprintf(pf, "Temps passé depuis l'exécution : %d s\n", (int)elapsed_seconds);
                                        fprintf(pf, "Nombre d'itérations : %ld\n\n", iterationsCompleted);

                                        fprintf(pf, "Valeures optimales trouvées :\n");
                                        fprintf(pf, "Pour la Fusée entière : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n", best_m_ail, best_n_ail, best_p_ail, best_E_ail);
                                        fprintf(pf, "Pour le Ship : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n\n", best_m_can, best_n_can, best_p_can, best_E_can);

                                        fprintf(pf, "Résultats (à vérifier avec StabTraj) :\n");
                                        fprintf(pf, "Pour la Fusée entière : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n", best_cn, best_cn0, best_ms_min, best_ms_max, best_ms_cn_min, best_ms_cn_max);
                                        fprintf(pf, "Pour le Ship : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n\n", best_cn_ship, best_cn0_ship, best_ms_min_ship, best_ms_max_ship, best_ms_cn_min_ship, best_ms_cn_max_ship);

                                        fprintf(pf, "Limites critiques : Cnmin : %f, Cnmax : %f, MSmin : %f, MSmax : %f, MSCnmin : %f, MSCnmax : %f\n\n", CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

                                        fprintf(pf, "La Fusée entière est : %s\n", best_stability);
                                        fprintf(pf, "Le Ship est : %s\n\n", best_stabilityShip);

                                        fprintf(pf, "----------------------------------------------------------------------\n\n");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // Calculate elapsed time
        clock_t current_time = clock();
        elapsed_seconds = (double)(current_time - start_time) / CLOCKS_PER_SEC;

        if (E_ail > 100) {
            double progress_fraction = (double)(E_ail - 100) / 300;
            double estimated_total_time = elapsed_seconds / progress_fraction;
            double remaining_time_seconds = estimated_total_time - elapsed_seconds;

            // Display the progress
            printf("Progrès : %d%%\n", (int) (100.0 * progress_fraction));
            printf("Temps passé depuis l'exécution : %d s\n", (int)elapsed_seconds);
            printf("Estimation du temps restant : %d s\n\n", (int)remaining_time_seconds);
        }
    }

    // Output the best values found
    printf("\nValeures optimales trouvées :\n");
    printf("Pour la Fusée entière : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n", best_m_ail, best_n_ail, best_p_ail, best_E_ail);
    printf("Pour le Ship : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n\n", best_m_can, best_n_can, best_p_can, best_E_can);

    printf("Résultats (à vérifier avec StabTraj) :\n");
    printf("Pour la Fusée entière : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n", best_cn, best_cn0, best_ms_min, best_ms_max, best_ms_cn_min, best_ms_cn_max);
    printf("Pour le Ship : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n\n", best_cn_ship, best_cn0_ship, best_ms_min_ship, best_ms_max_ship, best_ms_cn_min_ship, best_ms_cn_max_ship);

    printf("Limites critiques : Cnmin : %f, Cnmax : %f, MSmin : %f, MSmax : %f, MSCnmin : %f, MSCnmax : %f\n\n", CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

    printf("La Fusée entière est : %s\n", best_stability);
    printf("Le Ship est : %s\n", best_stabilityShip);

    // Output the best values found to the file
    fprintf(pf, "Temps passé depuis l'exécution : %d s\n", (int)elapsed_seconds);
    fprintf(pf, "Nombre d'itérations : %ld\n\n", iterationsCompleted);

    fprintf(pf, "Valeures optimales trouvées :\n");
    fprintf(pf, "Pour la Fusée entière : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n", best_m_ail, best_n_ail, best_p_ail, best_E_ail);
    fprintf(pf, "Pour le Ship : Emplanture : %d, Saumon : %d, Flèche : %d, Envergure : %d\n\n", best_m_can, best_n_can, best_p_can, best_E_can);

    fprintf(pf, "Résultats (à vérifier avec StabTraj) :\n");
    fprintf(pf, "Pour la Fusée entière : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n", best_cn, best_cn0, best_ms_min, best_ms_max, best_ms_cn_min, best_ms_cn_max);
    fprintf(pf, "Pour le Ship : Cn : %f, Cn0 : %f, MS_min : %f, MS_max : %f, MS_Cn_min : %f, MS_Cn_max : %f\n\n", best_cn_ship, best_cn0_ship, best_ms_min_ship, best_ms_max_ship, best_ms_cn_min_ship, best_ms_cn_max_ship);

    fprintf(pf, "Limites critiques : Cnmin : %f, Cnmax : %f, MSmin : %f, MSmax : %f, MSCnmin : %f, MSCnmax : %f\n\n", CritCnmin, CritCnmax, CritMsmin, CritMsmax, CritMsCnmin, CritMsCnmax);

    fprintf(pf, "La Fusée entière est : %s\n", best_stability);
    fprintf(pf, "Le Ship est : %s\n\n", best_stabilityShip);

    fprintf(pf, "-----------------------------END-----------------------------------------");

    fclose(pf);

    return 0;
}

// Function to get propeller characteristics based on its name
PropellerChar propellerCaracteristics(char* Propu) {
    PropellerChar propeller;

    // Original (Pro75-3G C)
    if (strcmp(Propu, "Pro75-3G") == 0) {
        propeller.Long_propu = 486;
        propeller.MpropuVide = 1.638;
        propeller.MpropuPlein = 3.511;
        propeller.XpropuVide = 243;
        propeller.XpropuPlein = 243;
    }
    // Pandora (Pro24-6G BS)
    else if (strcmp(Propu, "Pro24-6G") == 0) {
        propeller.Long_propu = 228;
        propeller.MpropuVide = 0.0843;
        propeller.MpropuPlein = 0.1599;
        propeller.XpropuVide = 114;
        propeller.XpropuPlein = 114;
    }
    // Barasinga (Pro54-5G C)
    else if (strcmp(Propu, "Pro54-5G") == 0) {
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
double calculateCni(const char* Type_masquage, double Q_int, double E_int, double D_ref, double D_int, double f_int, double m_int, double n_int) {
    if (Type_masquage[0] == 'B') {
        return 4 * Q_int * pow(E_int / D_ref, 2) * (1 + D_int / (2 * E_int + D_int)) / (1 + sqrt(1 + pow(2 * f_int / (m_int + n_int), 2)));
    }
    return 0;
}

double calculateCnail(double Q_ail, double E_ail, double D_ref, double D_ail, double f_ail, double m_ail, double n_ail) {
    return 4 * Q_ail * pow(E_ail / D_ref, 2) * (1 + D_ail / (2 * E_ail + D_ail)) / (1 + sqrt(1 + pow(2 * f_ail / (m_ail + n_ail), 2)));
}

double calculateCnai(double Cnail, double Cni) {
    return Cnail - Cni / 2;
}

double calculateCnc(const char* Type_masquage, double Q_can, double E_can, double D_ref, double D_can, double f_can, double m_can, double n_can) {
    if (Type_masquage[0] != 'M') {
        return 4 * Q_can * pow(E_can / D_ref, 2) * (1 + D_can / (2 * E_can + D_can)) / (1 + sqrt(1 + pow(2 * f_can / (m_can + n_can), 2)));
    }
    return 0;
}

double calculateCno(double D_og, double D_ref) {
    return 2 * pow(D_og / D_ref, 2);
}

double calculateCnj(double D2j, double D1j, double D_ref) {
    double Cnj = (D2j == 0) ? 0 : 2 * (pow(D2j / D_ref, 2) - pow(D1j / D_ref, 2));
    return Cnj;
}

double calculateCnr(double D2r, double D1r, double D_ref) {
    double Cnr = (D2r == 0) ? 0 : 2 * (pow(D2r / D_ref, 2) - pow(D1r / D_ref, 2));
    return Cnr;
}


double calculateCn(double Cni, double Cnail, double Cnai, double Cnc, double Cno, double Cnj, double Cnr) {
    return Cnai + Cnc + Cno + Cnj + Cnr;
}

double calculateCn0(double Cnail, double Cnc, double Cno, double Cnj, double Cnr) {
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

//Calculation of Xcp etc
double calculateXCpi(const char* Type_masquage, double X_int, double m_int, double p_int, double n_int) {
    if (Type_masquage[0] == 'B') {
        return X_int - m_int + p_int * (m_int + 2 * n_int) / (3 * (m_int + n_int)) + (m_int + n_int - m_int * n_int / (m_int + n_int)) / 6;
    }
    return 0;
}

double calculateXCpai(double XCpa, double Cnail, double XCpi, double Cni) {
    return (XCpa * Cnail - 0.5 * XCpi * Cni) / (Cnail - Cni / 2);
}

double calculateXCpc(const char* Type_masquage, double X_can, double m_can, double p_can, double n_can) {
    if (Type_masquage[0] != 'M') {
        return X_can - m_can + p_can * (m_can + 2 * n_can) / (3 * (m_can + n_can)) + (m_can + n_can - m_can * n_can / (m_can + n_can)) / 6;
    }
    return 0;
}

double calculateXCpj(double X_j, double l_j, double D1j, double D2j) {
    if (D2j != 0) {
        return X_j + l_j / 3 * (1 + 1 / (1 + D1j / D2j));
    }
    return 0;
}

double calculateXCpr(double X_r, double l_r, double D1r, double D2r) {
    if (D2r != 0) {
        return X_r + l_r / 3 * (1 + 1 / (1 + D1r / D2r));
    }
    return 0;
}

double calculateXCpo(const char* Forme_ogive, double Long_ogive) {
    if (strncmp(Forme_ogive, "Parab", 5) == 0) {
        return 0.5 * Long_ogive;
    } else if (strncmp(Forme_ogive, "Ogiv", 4) == 0) {
        return 7.0 / 15.0 * Long_ogive;
    } else if (strncmp(Forme_ogive, "Con", 3) == 0) {
        return 2.0 / 3.0 * Long_ogive;
    }
    return 0;
}

double calculateXCpa(double X_ail, double m_ail, double p_ail, double n_ail) {
    return X_ail - m_ail + p_ail * (m_ail + 2 * n_ail) / (3 * (m_ail + n_ail)) + (m_ail + n_ail - m_ail * n_ail / (m_ail + n_ail)) / 6;
}

double calculateXcgSans(const char* D12, double C12, double MasseVide, double XpropuRef, double Long_propu, double XpropuVide, double MasseSans, double MpropuPlein, double MpropuVide, double XpropuPlein) {
    if (strcmp(D12, "sans propu") == 0 || strcmp(D12, "without motor") == 0) {
        return C12;
    } else if (strcmp(D12, "avec propu vide") == 0 || strcmp(D12, "with empty motor") == 0) {
        return (C12 * MasseVide - (XpropuRef - Long_propu + XpropuVide) * MpropuVide) / MasseSans;
    } else if (strcmp(D12, "avec propu plein") == 0 || strcmp(D12, "with loaded motor") == 0) {
        return (C12 * MasseSans - (XpropuRef - Long_propu + XpropuPlein) * MpropuPlein) / MasseSans;
    }
    return -1; // Error case
}

double calculateMasseSans(const char* D11, double C11, double MpropuVide, double MpropuPlein) {
    if (strcmp(D11, "sans propu") == 0 || strcmp(D11, "without motor") == 0) {
        return C11 / 1000;
    } else if (strcmp(D11, "avec propu vide") == 0 || strcmp(D11, "with empty motor") == 0) {
        return C11 / 1000 - MpropuVide;
    } else if (strcmp(D11, "avec propu plein") == 0 || strcmp(D11, "with loaded motor") == 0) {
        return C11 / 1000 - MpropuPlein;
    }
    
    return -1; // Error case
}

double calculateXCp(double Cnai, double XCpai, double Cnc, double XCpc, double Cnj, double XCpj, double Cnr, double XCpr, double Cno, double XCpo) {
    double sumCnX = Cnai * XCpai + Cnc * XCpc + Cnj * XCpj + Cnr * XCpr + Cno * XCpo;
    double sumCn = Cnai + Cnc + Cnj + Cnr + Cno;
    return sumCn == 0 ? 0 : sumCnX / sumCn;
}

double calculateXCp0(double Cnail, double XCpa, double Cnc, double XCpc, double Cnj, double XCpj, double Cnr, double XCpr, double Cno, double XCpo) {
    double sumCnX = Cnail * XCpa + Cnc * XCpc + Cnj * XCpj + Cnr * XCpr + Cno * XCpo;
    double sumCn = Cnail + Cnc + Cnj + Cnr + Cno;
    return sumCn == 0 ? 0 : sumCnX / sumCn;
}

double calculateXcgPlein(double XcgSans, double MasseSans, double XpropuRef, double Long_propu, double XpropuPlein, double MpropuPlein) {
    return (XcgSans * MasseSans + (XpropuRef - Long_propu + XpropuPlein) * MpropuPlein) / (MasseSans + MpropuPlein);
}

double calculateXcgVide(double XcgSans, double MasseSans, double XpropuRef, double Long_propu, double XpropuVide, double MpropuVide) {
    return (XcgSans * MasseSans + (XpropuRef - Long_propu + XpropuVide) * MpropuVide) / (MasseSans + MpropuVide);
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

double valueInRange(double value, double rangeStart, double rangeEnd) {
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

// Know if the rocket is stable
char* checkStability(double Cn, double Cn0, double MS_min, double MS_max, double MS_Cn_min, double MS_Cn_max, double CritCnmin, double CritCnmax, double CritMsmin, double CritMsmax, double CritMsCnmin, double CritMsCnmax) {
    // Check if within stable range
    double isStable = CritCnmin < Cn && Cn0 < CritCnmax && CritMsmin < MS_min && MS_max < CritMsmax && CritMsCnmin < MS_Cn_min && MS_Cn_max < CritMsCnmax;

    // Check if outside stable range
    double isUnstable = Cn < CritCnmin || MS_min < CritMsmin || MS_Cn_min < CritMsCnmin;

    if (isStable) {
        return "STABLE";
    } else if (isUnstable) {
        return "INSTABLE";
    } else {
        return "SURSTABLE";
    }
}

