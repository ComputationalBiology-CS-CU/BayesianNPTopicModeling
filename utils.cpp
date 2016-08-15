//
//  utils.cpp
//  TopicModelforGenomics
//
//  Created by Jun Ho Yoon on 8/16/16.
//  Copyright © 2016 PeerLab. All rights reserved.
//

#include "utils.hpp"

double stirling_approx(double x) {
    if (x < 1.01) {
        return 0;
    } else {
        return x*log(x) - x + 0.5*log(6.283184*x);
    }
}

double log_beta(double x, double y) {
    return stirling_approx(x-1) + stirling_approx(y-1) - stirling_approx(x+y-1);
}

double log_beta(double x, double y, double z) {
    return stirling_approx(x-1) + stirling_approx(y-1) + stirling_approx(z-1) - stirling_approx(x+y+z-1);
}

double digamma(double x)
{
    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
        0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}