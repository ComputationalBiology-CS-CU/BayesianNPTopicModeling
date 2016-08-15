//
//  main.cpp
//  TopicModelforGenomics
//
//  Created by Jun Ho Yoon on 5/26/16.
//  Copyright © 2016 PeerLab. All rights reserved.
//

#include <string>
#include <iostream>
#include <sstream>
#include <random>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>

#include "utils.hpp"
#include "md_vector.h"

using namespace std;

int main() {
    //cout << "main\n";
    //  plate sizes
    const unsigned M = 10000;
    const unsigned N = 300;
    const unsigned T = 10;
    const unsigned S = 15;
    const unsigned K = 50;
    const unsigned J = 10;
    
    //  hyperparameters
    double alpha = 2;
    double omega = 8;
    double gamma = 2;
    double eta [3] = { 1, 1, 1 };
    
    //  declaration: variational parameters
    md_vector<double> v_alpha (S,J,2);
    md_vector<double> v_omega (K,2);
    md_vector<double> v_gamma (K,N,T,2);
    md_vector<double> v_eta (K,M,3);
    md_vector<double> phiE (S,N,J);
    md_vector<double> phiO (S,M,J);
    md_vector<double> v_nu (S,J,K);
    
    //  declaration: observations
    md_vector<int> E (S,N,T);
    md_vector<int> O_onehot (S,M,3);
    
    //  initialization: variational parameters
    string folder_path = "../data/data4/";
    
    double err = 0.1;
    unsigned seed1 = 150;
    default_random_engine generator(seed1);
    uniform_real_distribution<double> err_unif(-err, err);
    
    for (int k = 0; k < K; k++) {
        //if ((k+1)%10 == 0) cout << k << ' ';
        
        v_omega (k,0) = 1 + err_unif(generator);
        v_omega (k,1) = omega + err_unif(generator);
        
        for (int m = 0; m < M; m++) {
            v_eta (k,m,0) = eta[0] + err_unif(generator);
            v_eta (k,m,1) = eta[1] + err_unif(generator);
            v_eta (k,m,2) = eta[2] + err_unif(generator);
        }
        for (int n = 0; n < N; n++) {
            for (int t = 0; t < T; t++) {
                v_gamma (k,n,t,0) = gamma + err_unif(generator);
                v_gamma (k,n,t,1) = 1 + err_unif(generator);
            }
        }
    }
    for (int s = 0; s < S; s++) {
        for (int j = 0; j < J; j++) {
            v_alpha (s,j,0) = 1 + err_unif(generator);
            v_alpha (s,j,1) = alpha + err_unif(generator);
        }
    }
    cout << "read" << '\n';
    

    
    ifstream f_pi ( folder_path + "pi.csv" );
    int check_S = 0;
    while (f_pi) {
        //cout << "file" << '\n';
        string line;
        if (!getline (f_pi, line)) break;
        istringstream ss (line);
        
        int check_J = 0;
        while (ss) {
            //cout << "ss" << '\n';
            string sss;
            if (!getline (ss, sss, '\t')) break;
            
            for (int m = 0; m < M; m++) {
                phiO (check_S,m,check_J) = stod(sss) + abs( err_unif(generator) );
            }
            for (int n = 0; n < N; n++) {
                phiE (check_S,n,check_J) = stod(sss) + abs( err_unif(generator) );
            }
            check_J++;
        }
        check_S++;
    }
    //  normalize
    for (int s = 0; s < S; s++) {
        for (int n = 0; n < N; n++) {
            double j_sum = 0;
            for (int j = 0; j < J; j++) {
                j_sum += phiE (s,n,j);
            }
            for (int j = 0; j < J; j++) {
                phiE (s,n,j) = phiE (s,n,j) / j_sum;
            }
        }
    }
    for (int s = 0; s < S; s++) {
        for (int m = 0; m < M; m++) {
            double j_sum = 0;
            for (int j = 0; j < J; j++) {
                j_sum += phiO (s,m,j);
            }
            for (int j = 0; j < J; j++) {
                phiO (s,m,j) = phiO (s,m,j) / j_sum;
            }
        }
    }
    
    ifstream f_v ( folder_path+"v.csv" );
    int check_K = 0;
    while (f_v) {
        //cout << "file" << '\n';
        string line;
        if (!getline (f_v, line)) break;
        
        for (int s = 0; s < S; s++) {
            for (int j = 0; j < J; j++) {
                v_nu (s,j,check_K) = stod(line) + abs( err_unif(generator) );
            }
        }
        check_K++;
    }
    //  normalize
    for (int s = 0; s < S; s++) {
        for (int j = 0; j < J; j++) {
            double k_sum = 0;
            for (int k = 0; k < K; k++) {
                k_sum += v_nu (s,j,k);
            }
            for (int k = 0; k < K; k++) {
                v_nu (s,j,k) = v_nu (s,j,k) / k_sum;
            }
        }
    }
    
    //  read data E and O (and one-hot encoding)
    //  E
    string fname_E1 = folder_path + "E";
    string fname_E2 = ".csv";
    ostringstream oss;
    for (int t = 0; t < T; t++) {
        oss << t;
        ifstream f_E (fname_E1 + oss.str() + fname_E2);
        oss.str("");
        
        check_S = 0;
        while (f_E) {
            string line;
            if (!getline (f_E, line)) break;
            istringstream ss (line);
            
            int check_N = 0;
            while (ss) {
                string sss;
                if (!getline (ss, sss, '\t')) break;
                
                E (check_S,check_N,t) = stoi(sss);
                
                check_N++;
            }
            check_S++;
        }
    }
    //  O and one-hot encoding
    for (int s = 0; s < S; s++) {
        for (int m = 0; m < M; m++) {
            for (int l = 0; l < 3; l++) {
                O_onehot (s,m,l) = 0;
            }
        }
    }
    ifstream f_O ( folder_path + "O.csv" );
    check_S = 0;
    while (f_O) {
        string line;
        if (!getline (f_O, line)) break;
        istringstream ss (line);
        
        int check_M = 0;
        while (ss) {
            string sss;
            if (!getline (ss, sss, '\t')) break;
            
            O_onehot (check_S,check_M,stoi(sss)) = 1;
            
            check_M++;
        }
        check_S++;
    }
    
    //  variational updates
    double conv_crit = 0.0001;
    unsigned step = 1;
    double prev_ELBO = -10000000;
    vector<double> ELBOs;
    
    while (true) {
        cout << "phiE" << '\n';
        //  phiE update
        for (int s = 0; s < S; s++) {
            for (int n = 0; n < N; n++) {
                double dialpha_sum = 0;
                for (int j = 0; j < J; j++) {
                    phiE (s,n,j) = 0;
                    
                    if (j > 0) {
                        dialpha_sum += digamma(v_alpha (s,j-1,1)) - digamma(v_alpha (s,j-1,0) + v_alpha (s,j-1,1) );
                    }
                    for (int t = 0; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            phiE (s,n,j) += v_nu (s,j,k) * (E (s,n,t) * ( digamma(v_gamma (k,n,t,0)) - log(v_gamma (k,n,t,1)) ) - v_gamma (k,n,t,0)/v_gamma (k,n,t,1) - stirling_approx(E (s,n,t)));
                        }
                    }
                    phiE (s,n,j) += digamma(v_alpha (s,j,0)) - digamma(v_alpha (s,j,0) + v_alpha (s,j,1)) + dialpha_sum;
                }
            }
        }
        //  softmax normalize
        for (int s = 0; s < S; s++) {
            for (int n = 0; n < N; n++) {
                double jmax = -100000000;
                for (int j = 0; j < J; j++) {
                    if (phiE (s,n,j) > jmax) {
                        jmax = phiE (s,n,j);
                    }
                }
                double j_sum = 0;
                for (int j = 0; j < J; j++) {
                    j_sum += exp( phiE (s,n,j) - jmax );
                }
                for (int j = 0; j < J; j++) {
                    phiE (s,n,j) = exp( phiE (s,n,j) - jmax ) / j_sum;
                }
            }
        }
        
        //  phiO update
        cout << "phiO" << '\n';
        for (int s = 0; s < S; s++) {
            for (int m = 0; m < M; m++) {
                double dialpha_sum = 0;
                for (int j = 0; j < J; j++) {
                    phiO (s,m,j) = 0;
                    if (j > 0) {
                        dialpha_sum += digamma(v_alpha (s,j-1,1)) - digamma(v_alpha (s,j-1,0) + v_alpha (s,j-1,1) );
                    }
                    for (int l = 0; l < 3; l++) {
                        for (int k = 0; k < K; k++) {
                            phiO (s,m,j) += O_onehot (s,m,l) * v_nu (s,j,k) * ( digamma(v_eta (k,m,l)) - digamma(v_eta (k,m,0) + v_eta (k,m,1) + v_eta (k,m,2)) );
                        }
                    }
                    phiO (s,m,j) += digamma(v_alpha (s,j,0)) - digamma(v_alpha (s,j,0) + v_alpha (s,j,1)) + dialpha_sum;
                }
            }
        }
        //  softmax normalize
        for (int s = 0; s < S; s++) {
            for (int m = 0; m < M; m++) {
                double jmax = -100000000;
                for (int j = 0; j < J; j++) {
                    if (phiO (s,m,j) > jmax) {
                        jmax = phiO (s,m,j);
                    }
                }
                double j_sum = 0;
                for (int j = 0; j < J; j++) {
                    j_sum += exp( phiO (s,m,j) - jmax );
                }
                for (int j = 0; j < J; j++) {
                    phiO (s,m,j) = exp( phiO (s,m,j) - jmax ) / j_sum;
                }
            }
        }
        
        //  v_alpha update
        cout << "v_alpha" << '\n';
        for (int s = 0; s < S; s++) {
            for (int j = 0; j < J; j++) {
                v_alpha (s,j,0) = 1;
                v_alpha (s,j,1) = alpha;
                for (int m = 0; m < M; m++) {
                    v_alpha (s,j,0) += phiO (s,m,j);
                    
                    for (int i = j+1; i < J; i++) {
                        v_alpha (s,j,1) += phiO (s,m,i);
                    }
                }
                for (int n = 0; n < N; n++) {
                    v_alpha (s,j,0) += phiE (s,n,j);
                    
                    for (int i = j+1; i < J; i++) {
                        v_alpha (s,j,1) += phiE (s,n,i);
                    }
                }
            }
        }
        
        //  v_nu update
        cout << "v_nu" << '\n';
        for (int s = 0; s < S; s++) {
            for (int j = 0; j < J; j++) {
                double diomega_sum = 0;
                for (int k = 0; k < K; k++) {
                    v_nu (s,j,k) = 0;
                    if (k > 0) {
                        diomega_sum += digamma(v_omega (k-1,1)) - digamma(v_omega (k-1,0) + v_omega (k-1,1) );
                    }
                    for (int m = 0; m < M; m++) {
                        for (int l = 0; l < 3; l++) {
                            v_nu (s,j,k) += O_onehot (s,m,l) * phiO (s,m,j) * ( digamma(v_eta (k,m,l)) - digamma(v_eta (k,m,0) + v_eta (k,m,1) + v_eta (k,m,2)) );
                        }
                    }
                    for (int n = 0; n < N; n++) {
                        for (int t = 0; t < T; t++) {
                            v_nu (s,j,k) += phiE (s,n,j) * (E (s,n,t) * ( digamma(v_gamma (k,n,t,0)) - log(v_gamma (k,n,t,1)) ) - v_gamma (k,n,t,0)/v_gamma (k,n,t,1) - stirling_approx(E (s,n,t)));
                        }
                    }
                    v_nu (s,j,k) += digamma(v_omega (k,0)) - digamma(v_omega (k,0) + v_omega (k,1)) + diomega_sum;
                }
            }
        }
        //  softmax normalize
        for (int s = 0; s < S; s++) {
            for (int j = 0; j < J; j++) {
                double kmax = -100000000;
                for (int k = 0; k < K; k++) {
                    if (v_nu (s,j,k) > kmax) {
                        kmax = v_nu (s,j,k);
                    }
                }
                double k_sum = 0;
                for (int k = 0; k < K; k++) {
                    k_sum += exp( v_nu (s,j,k) - kmax );
                }
                for (int k = 0; k < K; k++) {
                    v_nu (s,j,k) = exp( v_nu (s,j,k) - kmax ) / k_sum;
                }
            }
        }
        
        //  v_omega update
        cout << "v_omega" << '\n';
        for (int k = 0; k < K; k++) {
            v_omega (k,0) = 1;
            v_omega (k,1) = omega;
            for (int s = 0; s < S; s++) {
                for (int j = 0; j < J; j++) {
                    v_omega (k,0) += v_nu (s,j,k);
                    for (int i = k+1; i < K; i++) {
                        v_omega (k,1) += v_nu (s,j,i);
                    }
                }
            }
        }
        
        //  v_gamma update
        cout << "v_gamma" << '\n';
        for (int k = 0; k < K; k++) {
            for (int n = 0; n < N; n++) {
                for (int t = 0; t < T; t++) {
                    v_gamma (k,n,t,0) = gamma;
                    v_gamma (k,n,t,1) = 1;
                    for (int s = 0; s < S; s++) {
                        for (int j = 0; j < J; j++) {
                            v_gamma (k,n,t,0) += E (s,n,t) * v_nu (s,j,k) * phiE (s,n,j);
                            v_gamma (k,n,t,1) += v_nu (s,j,k) * phiE (s,n,j);
                        }
                    }
                }
            }
        }
        
        //  v_eta update
        cout << "v_eta" << '\n';
        for (int k = 0; k < K; k++) {
            for (int m = 0; m < M; m++) {
                for (int l = 0; l < 3; l++) {
                    v_eta (k,m,l) = eta [l];
                    for (int s = 0; s < S; s++) {
                        for (int j = 0; j < J; j++) {
                            v_eta (k,m,l) += O_onehot (s,m,l) * v_nu (s,j,k) * phiO (s,m,j);
                        }
                    }
                }
            }
        }
        
        //ELBO calculation
        cout << "ELBOk" << '\n';
        double ELBO_kk = 0;
        double ELBO_km = 0;
        double ELBO_knt = 0;
        
        for (int k = 0; k < K; k++) {
            ELBO_kk += (omega - v_omega (k,1)) * (digamma(v_omega (k,1)) - digamma(v_omega (k,0) + v_omega (k,1))) + (1-v_omega (k,0)) * (digamma(v_omega (k,0)) - digamma(v_omega (k,0) + v_omega (k,1))) - log(omega) + log_beta(v_omega (k,0), v_omega (k,1));
            
            for (int m = 0; m < M; m++) {
                double di_eta = digamma(v_eta (k,m,0) + v_eta (k,m,1) + v_eta (k,m,2));
                for (int l = 0; l < 3; l++) {
                    ELBO_km += (eta [l] - v_eta (k,m,l)) * (digamma(v_eta (k,m,l)) - di_eta);
                }
                ELBO_km += log_beta(v_eta (k,m,0), v_eta (k,m,1), v_eta (k,m,2));
            }
            
            for (int n = 0; n < N; n++) {
                for (int t = 0; t < T; t++) {
                    ELBO_knt += (gamma - v_gamma (k,n,t,0)) * (digamma(v_gamma (k,n,t,0)) - log(v_gamma (k,n,t,1))) - v_gamma (k,n,t,0)/v_gamma (k,n,t,1) + v_gamma (k,n,t,0) - v_gamma (k,n,t,0) * log(v_gamma (k,n,t,1)) + stirling_approx(v_gamma (k,n,t,0)-1) - stirling_approx(gamma-1);
                }
            }
        }
        
        cout << "ELBOs" << '\n';
        double ELBO_sj = 0;
        double ELBO_sm = 0;
        double ELBO_sml = 0;
        double ELBO_sn = 0;
        double ELBO_snt = 0;
        
        double diomega_i = 0;
        for (int k = 0; k < K; k++) {
            if (k > 0) {
                diomega_i += digamma(v_omega (k-1,1)) - digamma(v_omega (k-1,0) + v_omega (k-1,1));
            }
            double diomega_k = digamma(v_omega (k,0)) - digamma(v_omega (k,0) + v_omega (k,1)) + diomega_i;
            for (int s = 0; s < S; s++) {
                for (int j = 0; j < J; j++) {
                    if (v_nu (s,j,k) > 0.0001) {
                        ELBO_sj += v_nu (s,j,k) * ( diomega_k - log(v_nu (s,j,k)) );
                    } else {
                        ELBO_sj += v_nu (s,j,k) * ( diomega_k );
                    }
                }
            }
        }
        for (int s = 0; s < S; s++) {
            for (int j = 0; j < J; j++) {
                ELBO_sj += (alpha - v_alpha (s,j,1)) * ( digamma(v_alpha (s,j,1)) - digamma(v_alpha (s,j,0) + v_alpha (s,j,1)) ) - (v_alpha (s,j,0) - 1) * (digamma(v_alpha (s,j,0)) - digamma(v_alpha (s,j,0) + v_alpha (s,j,1))) + log_beta(v_alpha (s,j,0), v_alpha (s,j,1));
            }
        }
        
        for (int s = 0; s < S; s++) {
            double dialpha_i = 0;
            for (int j = 0; j < J; j++) {
                if (j > 0) {
                    dialpha_i += digamma(v_alpha (s,j-1,1)) - digamma(v_alpha (s,j-1,0) + v_alpha (s,j-1,1));
                }
                double dialpha_sj = digamma(v_alpha (s,j,0)) - digamma(v_alpha (s,j,0) + v_alpha (s,j,1)) + dialpha_i;
                for (int m = 0; m < M; m++) {
                    if (phiO (s,m,j) > 0.0001) {
                        ELBO_sm += phiO (s,m,j) * ( dialpha_sj - log(phiO (s,m,j)) );
                    } else {
                        ELBO_sm += phiO (s,m,j) * ( dialpha_sj );
                    }
                    for (int k = 0; k < K; k++) {
                        double di_eta = digamma(v_eta (k,m,0) + v_eta (k,m,1) + v_eta (k,m,2));
                        for (int l = 0; l < 3; l++) {
                            ELBO_sml += v_nu (s,j,k) * phiO (s,m,j) * O_onehot (s,m,l) * ( digamma(v_eta (k,m,l)) - di_eta );
                        }
                    }
                }
                for (int n = 0; n < N; n++) {
                    if (phiE (s,n,j) > 0.0001) {
                        ELBO_sn += phiE (s,n,j) * ( dialpha_sj - log(phiE (s,n,j)) );
                    } else {
                        ELBO_sn += phiE (s,n,j) * ( dialpha_sj );
                    }
                    for (int k = 0; k < K; k++) {
                        for (int t = 0; t < T; t++) {
                            ELBO_snt += v_nu (s,j,k) * phiE (s,n,j) * ( E (s,n,t) * (digamma(v_gamma (k,n,t,0)) - log(v_gamma (k,n,t,1))) - v_gamma (k,n,t,0)/v_gamma (k,n,t,1) - stirling_approx(E (s,n,t)) );
                        }
                    }
                }
            }
        }
        cout << ELBO_kk << " " << ELBO_km << " " << ELBO_knt << " " << ELBO_sj << " " << ELBO_sm << " " << ELBO_sml << " " << ELBO_sn << " " << ELBO_snt << '\n';
        double ELBO = ELBO_kk + ELBO_km + ELBO_knt + ELBO_sj + ELBO_sm + ELBO_sml + ELBO_sn + ELBO_snt;
        
        //convergence test
        double conv = abs(prev_ELBO - ELBO) / abs(prev_ELBO);
        cout << step << " ELBO: " << ELBO << ", convergence: " << conv << '\n';
        if (conv < conv_crit) {
            cout << "converged." << '\n';
            break;
        }
        ELBOs.push_back(ELBO);
        prev_ELBO = ELBO;
        step++;
    }
    
    return 0;
}
