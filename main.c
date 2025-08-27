#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Função para simulação em malha aberta
void malha_aberta(double y0, double r, double T_int, int n_int, double K, double tau, double *y_aberto)
{
    double y = y0;
    for (int k = 0; k < n_int; k++)
    {
        double dy = (-y + K * r) / tau;
        y += dy * T_int;
        y_aberto[k] = y;
    }
}

// Função para simulação em malha fechada sem controlador
void malha_fechada_sem_controlador(double y0, double r, double T_int, int n_int, double K, double tau, double *y_hist)
{
    double y = y0;
    for (int k = 0; k < n_int; k++)
    {
        double u = r - y;
        double dy = (-y + K * u) / tau;
        y += dy * T_int;
        y_hist[k] = y;
    }
}

// Função para simulação em malha fechada com PI
void malha_fechada(double y0, double r, double T_int, double integral_error0, int n_int, double T_s, double Kp, double Ki, double tau, double u0, double *y_hist, double *u_hist)
{
    double y = y0;
    double u = u0;
    double integral_error = integral_error0;
    int amostras_por_Ts = (int)(T_s / T_int);

    for (int k = 0; k < n_int; k++)
    {
        if (k % amostras_por_Ts == 0)
        {
            double e = r - y;
            double i_candidate = integral_error + Ki * T_s * e;
            double u_unsat = Kp * e + i_candidate;
            u = u_unsat; // Sem saturação
            if (u == u_unsat)
            {
                integral_error = i_candidate;
            }
        }
        double dy = (-y + 2.0 * u) / tau;
        y += dy * T_int;
        y_hist[k] = y;
        u_hist[k] = u;
    }
}

int main()
{
    // Parâmetros da planta
    double K = 2.0;
    double tau = 10.0;
    double tau_min = tau;

    // Configurações de simulação
    double T_int = tau_min / 100.0;
    double T_s = tau_min / 10.0;
    double T_final = 100.0;
    int n_int = (int)(T_final / T_int);

    // Parâmetros do controlador PI
    double Kp = 3.0;
    double Ki = 0.4;

    // Variáveis de estado
    double y = 0.0;
    double u = 0.0;
    double integral_error = 0.0;
    double r = 20.0;

    // Alocação dos vetores de histórico
    double *time = (double *)malloc(n_int * sizeof(double));
    double *y_hist = (double *)malloc(n_int * sizeof(double));
    double *u_hist = (double *)malloc(n_int * sizeof(double));
    double *y_aberto = (double *)malloc(n_int * sizeof(double));
    double *y_fechado_sc = (double *)malloc(n_int * sizeof(double));

    for (int i = 0; i < n_int; i++)
    {
        time[i] = i * T_int;
    }

    // Simulação do sistema em malha aberta
    malha_aberta(y, r, T_int, n_int, K, tau, y_aberto);

    // Simulação do sistema em malha fechada sem controlador
    malha_fechada_sem_controlador(0.0, r, T_int, n_int, K, tau, y_fechado_sc);

    // Simulação do sistema em malha fechada com PI
    malha_fechada(0.0, r, T_int, 0.0, n_int, T_s, Kp, Ki, tau, 0.0, y_hist, u_hist);

    // Exemplo de saída simples (substitua por gráficos em Python, se desejar)
    printf("Tempo\tAberta\tFechada_SC\tFechada_PI\tControle_PI\n");
    for (int i = 0; i < n_int; i += n_int / 20)
    { // Mostra 20 pontos
        printf("%.1f\t%.2f\t%.2f\t\t%.2f\t\t%.2f\n", time[i], y_aberto[i], y_fechado_sc[i], y_hist[i], u_hist[i]);
    }

    // Liberação de memória
    free(time);
    free(y_hist);
    free(u_hist);
    free(y_aberto);
    free(y_fechado_sc);

    return 0;
}