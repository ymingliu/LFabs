
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

double myround(double a, int n)
{
    double s;
    if (a<0) {
        s=-a;
    }else{
        s=a;
    }
    s = s * pow(10.0, 1.0*n);
    s = s + 0.5;
    s = (int)s;
    s = s / pow(10.0, 1.0*n);
    if (a<0) {
        return -s;
    }else{
        return s;
    }
}

void lmL(double *Loss, double *LLoss, double *y, double *xbeta, int *param, double tau, double difflaplace)
{
    int i;
    double sum=0.0, a;
    int n = param[0];
    
    for (i = 0; i < n; i++) {
        a = y[i] - xbeta[i];
        sum += a*a;
    }
    
    *Loss = (1.0/n)*sum;
    *LLoss = *Loss + tau * difflaplace;
}

void dlmL(double *der, double *y, double *xbeta, double *x, int *param,
          double tau, double *gammabeta, double *Lbeta)
{
    int i, u;
    
    int n = param[0];
    int p = param[1];
    
    for (u = 0; u < p; u++) {
        for (i = 0; i < n; i++) {
            der[u] += 1.0*x[i*p + u] * (y[i] - xbeta[i]);
        }
        der[u] /= -n;
    }
    
    if (tau!=0) {
        for (i = 0; i < p; i++) der[i] += tau * (gammabeta[i] - Lbeta[i]);
    }
}

void coxL(double *Loss, double *LLoss, double *xbeta, int *param, int *delta, 
    double tau, double difflaplace)
{
    int i, j;
    double sum=0.0, d;
    int n = param[0];

    for (i = 0; i < n; i++){
        if(delta[i]){
            d = 0.0;
            for (j=0; j<=i; j++){
                d += exp(xbeta[j]);
            }
            sum += 1.0*(xbeta[i]-log(d));
        }
    }

    *Loss = -(1.0/n)*sum;
    *LLoss = *Loss + tau * difflaplace;
}
    
void dcoxL(double *der, double *xbeta, double *x, int *param, 
    int *delta, double tau, double *gammabeta, double *Lbeta)
{
    int i, j, u;
    double a, d;

    int n = param[0];
    int p = param[1];

    for (u = 0; u < p; u++) {
        for (i = 0; i < n; i++) {
            if(delta[i]){
                a = d = 0.0;
                for (j = 0; j <= i; j++){
                    a += x[j*p + u]*exp(xbeta[j]);
                    d += exp(xbeta[j]);
                }
                der[u] += 1.0*(x[i*p + u] - a/d);
            }
        }
        der[u] /= -n;
    }

    if (tau!=0) {
        for (i = 0; i < p; i++) der[i] += tau * (gammabeta[i] - Lbeta[i]);
    }
}
    
void sprL(double *Loss, double *LLoss, double *xbeta, int *param, int *delta, 
    double sigma, double tau, double difflaplace)
{
    int i, j;
    double sum=0, u=0;
    int n = param[0];
    for (i=0; i<n; i++) {
        for (j = i+1; j < n; j++) {
            if (delta[j]) {
                u = (xbeta[i] - xbeta[j])/sigma;
                sum += 1.0/(1.0 + exp(-1.0*u));
            }
        }
    }
    *Loss = -sum/(n*(n-1));
    *LLoss = *Loss + tau * difflaplace;
}
    
void dsprL(double *der, double *xbeta, double *x, int *param, int *delta, 
    double sigma, double tau, double *gammabeta, double *Lbeta)
{
    int i, j, k;
    int n = param[0];
    int p = param[1];
    double u, *xi, *xj;

    for (k = 0; k < p; k++) der[k] = 0.0;
    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            if (delta[j]) {
                u = (xbeta[i] - xbeta[j])/sigma;
                if (u >= -500) {
                    u = exp(-1.0*u);
                    u = u/( (1.0 + u)*(1.0 + u) );
                    xi = x + i*p;
                    xj = x + j*p;
                    for (k = 0; k < p; k++) {
                        der[k] += u * (xi[k] - xj[k]);
                    }
                }
            }
        }
    }
    for (k = 0; k < p; k++) der[k] /= -n*(n-1)*sigma;
    if (tau!=0) {
        for (k = 0; k < p; k++) der[k] = der[k] + tau * (gammabeta[k] - Lbeta[k]);
    }
}
    
void loglikhd(double *corr, double *corr0, int is_loss, double *y, double *xbeta, int *param, int *delta,
    double sigma, double *tau, double difflaplace)
{
    if (is_loss == 0) {
        sprL(corr, corr0, xbeta, param, delta, sigma, tau[0], difflaplace);
    } else if (is_loss == 1) {
        coxL(corr, corr0, xbeta, param, delta, tau[0], difflaplace);
    } else {
        lmL(corr, corr0, y, xbeta, param, tau[0], difflaplace);
    }
}
    
void dloglikhd(double *der, int is_loss, double *y, double *x, double *xbeta, int *param,
    int *delta, double sigma, double *tau, double *gammabeta, double *Lbeta)
{
    if (is_loss == 0) {
        dsprL(der, xbeta, x, param, delta, sigma, tau[0], gammabeta, Lbeta);
    } else if (is_loss == 1) {
        dcoxL(der, xbeta, x, param, delta, tau[0], gammabeta, Lbeta);
    } else {
        dlmL(der, y, xbeta, x, param, tau[0], gammabeta, Lbeta);
    }
}

double calculate_bic(double *Loss, int *param, int is_loss)
{
    int n = param[0];
    int p = param[1];
    int df = param[3];

    if (is_loss == 0) {
        return -2*n*log(-(*Loss)) + df*log(n);
    } else if (is_loss == 1) {
        return 2*n*(*Loss) + df*log(n);
    } else {
        return 2*n*(*Loss) + df*log(n) + 2*lchoose(p, df);
    }
}

// update sparse index
void usi(int *s_i, int *s_j, int *tt_b, int *tt_a, int *act, int ns, int iter)
{
    int i;
    s_i += *tt_a;
    s_j += *tt_a;
    for (i = 0; i < ns; ++i)
    {
        s_i[i] = act[i];
        s_j[i] = iter;
    }
    *tt_b = *tt_a;
    *tt_a += ns;
}

void update_Lbate(double *edgeweight, int *edgeweightLocRow, 
    int *edgeweightLocCol, double *Lbeta, double value, int k, int *param)
{
    int i;

    for (i = 0; i < param[2]; ++i)
    {
        if ( edgeweightLocCol[i] == k ) {
            Lbeta[ edgeweightLocRow[i] ] += edgeweight[i] * value;
        }
        if (edgeweightLocRow[i] == k) {
            Lbeta[ edgeweightLocCol[i] ] += edgeweight[i] * value;
        }
    }
}

// take an initial step
void Initial(double *der, int is_loss, double *y, double *x, double *xbeta, int *param,
    int *delta, double sigma, double *tau, double *gammabeta, double *Lbeta, 
    double *weight, double eps, double *beta, double *Loss, double *LLoss, 
    double *gamma, double *edge, int *edgeRow, int *edgeCol, double *lambda, 
    int *direction, int *active, double *bic, int *nactive, double *difflaplace)
{
    // der doesn't need initial value
    int i, k=0;
    int n = param[0], p = param[1];
    double temp, temp2, value;

    // calculate the derivative
    dloglikhd(der, is_loss, y, x, xbeta, param, delta, sigma, tau, gammabeta, Lbeta);
    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i){
        value = fabs(der[i])/weight[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }
    // calculate increment, update xb, lambda, beta and loss_after
    value = eps/weight[k];
    if (der[k] > 0.0) value *= -1.0;
    *beta = value;
    // calculate Loss and LLoss for beta = 0.
    loglikhd(&temp, &temp2, is_loss, y, xbeta, param, delta, sigma, tau, *difflaplace);

    // update xbeta, Loss, LLoss, difflaplace, gammabeta, Lbeta, lambda
    for (i = 0; i < n; ++i)
        xbeta[i] = x[i*p+k]*value;
    gammabeta[k] = gamma[k]*value;
    update_Lbate(edge, edgeRow, edgeCol, Lbeta, value, k, param);
    *difflaplace = gamma[k]*value*value/2.0;
    loglikhd(Loss, LLoss, is_loss, y, xbeta, param, delta, sigma, tau, *difflaplace);
    *lambda = (temp2 - *LLoss)/eps;

    *direction = 1;
    param[3]   = 1;
    *active    = k;
    *bic       = calculate_bic(Loss, param, is_loss);
    *nactive   = 1;
}

int Backward(double *der, int is_loss, double *y, double *x, double *xbeta, int *param,
    int *delta, double sigma, double *tau, double *gammabeta, double *Lbeta, 
    double *weight, double eps, double *beta, double *Loss, double *LLoss, 
    double *gamma, double *edge, int *edgeRow, int *edgeCol, 
    double *lambda, int *direction, int *active,double *bic, int *nactive, 
    double *beta0, double *xbtemp, double xi, double *difflaplace)
{
    int i, c, k;
    int n = param[0], p = param[1], ns = param[3];
    double temp, value;
    
    // update der
    dloglikhd(der, is_loss, y, x, xbeta, param, delta, sigma, tau, gammabeta, Lbeta);

    // find a backward direction
    k = 0;
    c = active[0];
    temp = der[c] / weight[c];
    if (beta[0] < 0) temp *= -1.0;
    beta0[0]  = beta[0];
    for (i = 1; i < ns; ++i){
        c = active[i];
        value = der[c] / weight[c];
        if (beta[i] < 0) value *= -1.0;
        if (value > temp){
            temp = value;
            k = i;
        }
        beta0[i]  = beta[i];
    }
    // try a backward step.
    c = active[k];
    value = eps/weight[c];
    if (beta[k] > 0) value *= -1.0;
    beta0[k] += value;
    for (i = 0; i < n; ++i)
        xbtemp[i] = xbeta[i] + x[i*p+c]*value;
    temp = *difflaplace + gamma[c]*value*value/2.0 + value*(gammabeta[c] - Lbeta[c]);
    loglikhd(Loss+1, LLoss+1, is_loss, y, xbtemp, param, delta, sigma, tau, temp);
    if (myround(LLoss[1] - LLoss[0] - (*lambda)*eps/weight[c], 6) < (-1.0*xi))
    {
        // adopt a backward step.
        // update direction, bic, nactive, lambda, difflaplace, gammabeta, Lbeta, 
        // active, param[3], beta, xbeta, Loss, LLoss
        direction[1] = -1;
        bic[1]       = calculate_bic(Loss+1, param, is_loss);
        nactive[1]   = nactive[0];
        lambda[1]    = lambda[0];
        *difflaplace = temp;
        for (i = 0; i < n; ++i)
            xbeta[i] = xbtemp[i];
        update_Lbate(edge, edgeRow, edgeCol, Lbeta, value, c, param);
        gammabeta[c] += gamma[c]*value;
        // test if vanish
        if (fabs(beta0[k]) < 1.0*xi) {
            for (i = k; i < param[3]; ++i)
            {
                active[i] = active[i+1];
                beta0[i] = beta0[i+1];
            }
            param[3]--;
            nactive[1]--;
        }
        return 0;
    }
    beta0[k] -= value;
    return 1;
}

void Forward(double *der, int is_loss, double *y, double *x, double *xbeta, int *param,
    int *delta, double sigma, double *tau, double *gammabeta, double *Lbeta, 
    double *weight, double eps, double *beta, double *Loss, double *LLoss, 
    double *gamma, double *edge, int *edgeRow, int *edgeCol, double *lambda, 
    int *direction, int *active, double *bic, int *nactive, double *difflaplace, double xi)
{
    // der doesn't need initial value
    int i, j, k=0;
    int n = param[0], p = param[1], ns = param[3];
    double temp, value;
    
    // d has calculated in Backward.
    // update der
    // find a forward direction
    temp = 0.0;
    for (i = 0; i < p; ++i){
        value = fabs(der[i])/weight[i];
        if (value > temp){
            temp = value;
            k = i;
        }
    }
    // calculate increment
    value = eps/weight[k];
    if (der[k] > 0.0) value *= -1.0;
    // beta has been assigned in backward
    // update beta, active, nactive, param[3]
    if (k > active[ns-1])
    {
        active[ns] = k;
        beta[ns] = value;
        nactive[1] = nactive[0]+1;
        param[3]++;
    } else {
        for (i = 0; i < ns; ++i)
        {
            if (active[i] < k) continue;
            if (active[i] == k)
            {
                beta[i] += value;
                nactive[1] = nactive[0];
            } else {
                for (j = ns; j > i; --j)
                {
                    active[j] = active[j-1];
                    beta[j] = beta[j-1];
                }
                active[i] = k;
                beta[i] = value;
                param[3]++;
                nactive[1] = nactive[0]+1;
            }
            break;
        }
    }
    // update xbeta, Loss, LLoss, difflaplace, gammabeta, Lbeta, lambda, direction, bic
    for (i = 0; i < n; ++i)
        xbeta[i] += x[i*p+k]*value;
    *difflaplace += gamma[k]*value*value/2 + value*(gammabeta[k] - Lbeta[k]);
    gammabeta[k] += gamma[k]*value;
    update_Lbate(edge, edgeRow, edgeCol, Lbeta, value, k, param);
    loglikhd(Loss+1, LLoss+1, is_loss, y, xbeta, param, delta, sigma, tau, difflaplace[0]);
    temp = (LLoss[0] - LLoss[1] - xi)/eps;
    if (temp < *lambda) lambda[1] = temp;
    else lambda[1] = lambda[0];
    // test if vanish
    direction[1] = 1;
    bic[1]       = calculate_bic(Loss, param, is_loss);
}


void LFabs_single_tau(double *der, int is_loss, double *y, double *x, double *xbeta, int *param, int *delta, double sigma, double *tau, double *gammabeta, double *Lbeta,
                      double *weight, double eps, double *beta, double *Loss, double *LLoss,
                      double *gamma, double *edge, int *edgeRow, int *edgeCol, double *lambda,
                      int *direction, int *active, double *bic, int*nactive, int iter, double *xbtemp,
                      double xi, int max_s, int *sparse_i, int *sparse_j, double lam_m, int stoping)
{
    int i, k;
    double difflaplace = 0.0;
    
    // step 1: initial step (forward)
    Initial(der, is_loss, y, x, xbeta, param, delta, sigma, tau, gammabeta, Lbeta,
            weight, eps, beta, Loss, LLoss, gamma, edge, edgeRow, edgeCol,
            lambda, direction, active, bic, nactive, &difflaplace);
    int tt_act_b = 0;
    int tt_act_a = 0;
    usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], 0);
    
    // step 2: forward and backward
    for (i = 0; i < iter-1; ++i)
    {
        k = Backward(der, is_loss, y, x, xbeta, param, delta, sigma, tau, gammabeta,
                     Lbeta, weight, eps, beta+tt_act_b, Loss+i, LLoss+i, gamma, edge, edgeRow,
                     edgeCol, lambda+i, direction+i, active, bic+i, nactive+i, beta+tt_act_a, xbtemp, xi, &difflaplace);
        if (k) {
            Forward(der, is_loss, y, x, xbeta, param, delta, sigma, tau, gammabeta, Lbeta,
                    weight, eps, beta+tt_act_a, Loss+i, LLoss+i, gamma, edge, edgeRow, edgeCol,
                    lambda+i, direction+i, active, bic+i, nactive+i, &difflaplace, xi);
        }
        usi(sparse_i, sparse_j, &tt_act_b, &tt_act_a, active, param[3], i+1);
        
        if ( stoping && (lambda[i+1] <= lambda[0] * lam_m) ) {
            iter = i+2;
            if (lambda[i+1] < 0)
            {
                iter--;
                tt_act_a -= param[3];
            }
            break;
        }
        if (param[3] > max_s) {
//            Rprintf("Warning! Max nonzero number is larger than predetermined threshold. Program ended early.\n");
            iter = i+2;
            break;
        }
        if (i == iter-2) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }
    
    param[4] = tt_act_a;
    param[5] = iter;
}

SEXP LFabs(SEXP Y, SEXP X, SEXP Weight, SEXP Delta, SEXP Epsilon, SEXP Lam_min,
           SEXP Sigma, SEXP Xi, SEXP Stoping, SEXP Iter, SEXP Param, SEXP Is_loss,
           SEXP Edge, SEXP EdgeRow, SEXP EdgeCol, SEXP Gamma, SEXP Max_S, SEXP Tau)
{
    int i, n, p, ltau, ntau, iter, max_s;
    double *y, *x, *weight, eps, lam_m, xi, sigma, *gamma, *edge, *tau;
    int *param, stoping, is_loss, *edgeRow, *edgeCol, *delta;
    
    param     = INTEGER(Param);
    n         = param[0];
    p         = param[1];
    ntau      = param[6];
    y         = REAL(Y);
    x         = REAL(X);
    weight    = REAL(Weight);
    eps       = REAL(Epsilon)[0];
    lam_m     = REAL(Lam_min)[0];
    xi        = REAL(Xi)[0];
    stoping   = INTEGER(Stoping)[0];
    iter      = INTEGER(Iter)[0];
    max_s     = INTEGER(Max_S)[0];
    is_loss   = INTEGER(Is_loss)[0];
    delta     = INTEGER(Delta);
    sigma     = REAL(Sigma)[0];
    gamma     = REAL(Gamma);
    edge      = REAL(Edge);
    edgeRow   = INTEGER(EdgeRow);
    edgeCol   = INTEGER(EdgeCol);
    tau       = REAL(Tau);
    
    double *der, *beta, *lambda, *xbeta, *bic, *loss;
    double *lloss, *xbtemp, *gammabeta, *Lbeta;
    int *direction, *nactive, *active;
    int *sparse_i, *sparse_j;

    beta      = (double*)malloc(sizeof(double)*iter*max_s);
    sparse_i  =    (int*)calloc(iter*max_s, sizeof(int));
    sparse_j  =    (int*)calloc(iter*max_s, sizeof(int));
    lambda    = (double*)calloc(iter, sizeof(double));
    direction =    (int*)malloc(sizeof(int)   *iter);
    bic       = (double*)malloc(sizeof(double)*iter);
    loss      = (double*)malloc(sizeof(double)*iter);
    lloss     = (double*)malloc(sizeof(double)*iter);
    nactive   =    (int*)malloc(sizeof(int)   *iter);
    xbtemp    = (double*)malloc(sizeof(double)  *n);
    xbeta     = (double*)calloc(n, sizeof(double)); // x^T beta.
    der       = (double*)calloc(p, sizeof(double)); // 1st order Taylor Formula.
    active    =    (int*)calloc(max_s+1, sizeof(int));
    gammabeta = (double*)calloc(p, sizeof(double));
    Lbeta     = (double*)calloc(p, sizeof(double));
    
    SEXP AllResult;
    PROTECT(AllResult = allocVector(VECSXP, ntau));
    
    for (ltau = 0; ltau < ntau; ++ltau)
    {
        LFabs_single_tau(der, is_loss, y, x, xbeta, param, delta, sigma, tau+ltau,
                         gammabeta, Lbeta, weight, eps, beta, loss, lloss, gamma, edge,
                         edgeRow, edgeCol, lambda, direction, active, bic, nactive, iter,
                         xbtemp, xi, max_s, sparse_i, sparse_j, lam_m, stoping);
        
        int tt_act_a = param[4];
        int niter = param[5];
        
        SEXP Beta, Lam, Drct, Loops, LLoss, BIC, Loss, Nactive, Idi, Idj, Result, R_names;
        char *names[10] = {"beta", "lambda", "direction", "iter", "bic", "loss",
            "losslaplace", "n_active", "index_i", "index_j"};
        PROTECT(Beta    = allocVector(REALSXP, tt_act_a));
        PROTECT(Idi     = allocVector(INTSXP,  tt_act_a));
        PROTECT(Idj     = allocVector(INTSXP,  tt_act_a));
        PROTECT(Lam     = allocVector(REALSXP, niter));
        PROTECT(BIC     = allocVector(REALSXP, niter));
        PROTECT(Loss    = allocVector(REALSXP, niter));
        PROTECT(LLoss   = allocVector(REALSXP, niter));
        PROTECT(Drct    = allocVector(INTSXP,  niter));
        PROTECT(Nactive = allocVector(INTSXP,  niter));
        PROTECT(Loops   = allocVector(INTSXP,  1));
        PROTECT(Result  = allocVector(VECSXP,  10));
        PROTECT(R_names = allocVector(STRSXP,  10));
        
        for(i = 0; i < 10; ++i) SET_STRING_ELT(R_names, i,  mkChar(names[i]));
        INTEGER(Loops)[0] = niter;
        for (i = 0; i < tt_act_a; ++i)
        {
            REAL(Beta)[i]   = beta[i];
            INTEGER(Idi)[i] = sparse_i[i];
            INTEGER(Idj)[i] = sparse_j[i];
        }
        for (i = 0; i < niter; ++i)
        {
            REAL(BIC)[i]        = bic[i];
            REAL(Loss)[i]       = loss[i];
            REAL(LLoss)[i]      = lloss[i];
            REAL(Lam)[i]        = lambda[i];
            INTEGER(Drct)[i]    = direction[i];
            INTEGER(Nactive)[i] = nactive[i];
        }
        
        SET_VECTOR_ELT(Result, 0, Beta);
        SET_VECTOR_ELT(Result, 1, Lam);
        SET_VECTOR_ELT(Result, 2, Drct);
        SET_VECTOR_ELT(Result, 3, Loops);
        SET_VECTOR_ELT(Result, 4, BIC);
        SET_VECTOR_ELT(Result, 5, Loss);
        SET_VECTOR_ELT(Result, 6, LLoss);
        SET_VECTOR_ELT(Result, 7, Nactive);
        SET_VECTOR_ELT(Result, 8, Idi);
        SET_VECTOR_ELT(Result, 9, Idj);
        setAttrib(Result, R_NamesSymbol, R_names);
        
        SET_VECTOR_ELT(AllResult, ltau, Result);
        
        if (ltau < ntau-1)
        {
            param[3] = 0;
            for (i = 0; i < tt_act_a; ++i)
            {
                beta[i]     = 0.0;
                sparse_i[i] = 0;
                sparse_j[i] = 0;
            }
            
            for (i = 0; i < niter; ++i)
            {
                bic[i]       = 0.0;
                loss[i]      = 0.0;
                lloss[i]     = 0.0;
                lambda[i]    = 0.0;
                direction[i] = 0;
                nactive[i]   = 0;
            }
            
            for (i = 0; i < n; ++i) xbeta[i] = 0.0; // Need to be initialized.
            for (i = 0; i < max_s+1; ++i) active[i] = 0;
            for (i = 0; i < p; ++i)
            {
                gammabeta[i] = 0.0;  // Need to be initialized.
                Lbeta[i] = 0.0;
            }
        }
    }
    
    free(beta      );
    free(sparse_i  );
    free(sparse_j  );
    free(lambda    );
    free(direction );
    free(bic       );
    free(loss      );
    free(lloss     );
    free(nactive   );
    free(xbtemp    );
    free(xbeta     );
    free(der       );
    free(active    );
    free(gammabeta );
    free(Lbeta     );
    
    UNPROTECT(12*ntau+1);
    return AllResult;
}




    
    
 


























