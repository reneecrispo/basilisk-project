#include <stdio.h>
#include <math.h>

#define MAX_LINES 10000

// Funzione per caricare i dati (tempo e velocità) da un file
int load_data(const char *filename, double *time, double *velocity, int max_lines) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return 0;
    }

    int n_lines = 0;
    while (fscanf(file, "%lf %lf", &time[n_lines], &velocity[n_lines]) == 2) {
        n_lines++;
        if (n_lines >= max_lines) break;
    }
    fclose(file);
    return n_lines;
}

// Funzione per calcolare l'errore in norma L2 tra due set di dati di velocità
double compute_L2_error(double *time1, double *vel1, int n1, double *time2, double *vel2, int n2) {
    double error = 0.0;
    int i1 = 0, i2 = 0;

    // Calcoliamo l'errore solo per i punti comuni, ma gestiamo i punti extra
    while (i1 < n1 || i2 < n2) {
        if (i1 < n1 && i2 < n2) {
            // Se i tempi sono abbastanza vicini, calcola l'errore
            if (fabs(time1[i1] - time2[i2]) < 1e-6) {
                error += pow(vel1[i1] - vel2[i2], 2);  // Calcola il quadrato della differenza
                i1++;
                i2++;
            } else if (time1[i1] < time2[i2]) {
                // Se il tempo nel primo file è minore, considera solo il primo file
                i1++;
            } else {
                // Se il tempo nel secondo file è minore, considera solo il secondo file
                i2++;
            }
        } else if (i1 < n1) {
            // Se ci sono più righe nel primo file, calcoliamo l'errore per il primo file
            i1++;
        } else if (i2 < n2) {
            // Se ci sono più righe nel secondo file, calcoliamo l'errore per il secondo file
            i2++;
        }
    }

    return sqrt(error);  // Ritorna la radice quadrata della somma dei quadrati
}

int main() {
    // Variabili per memorizzare i dati (tempo e velocità)
    double time1[MAX_LINES], vel1[MAX_LINES];
    double time2[MAX_LINES], vel2[MAX_LINES];
    int n1, n2;

    // Carica i dati dai due file
    n1 = load_data("velocity_log_1", time1, vel1, MAX_LINES);
    n2 = load_data("velocity_log_2", time2, vel2, MAX_LINES);

    if (n1 == 0 || n2 == 0) {
        fprintf(stderr, "Errore nel caricare i dati dai file.\n");
        return 1;
    }

    // Calcola l'errore in norma L2 tra le due velocità
    double error = compute_L2_error(time1, vel1, n1, time2, vel2, n2);
    printf("Errore in norma L2 tra le due velocità: %g\n", error);

    return 0;
}
