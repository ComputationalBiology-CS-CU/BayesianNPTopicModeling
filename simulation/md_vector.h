//
//  md_vector.h
//  TopicModelforGenomics
//
//  Created by Jun Ho Yoon on 8/15/16.
//  Copyright © 2016 PeerLab. All rights reserved.
//

#ifndef md_vector_h
#define md_vector_h

#include <vector>

template <typename T>
class md_vector {
    size_t d1, d2, d3, d4; //  dimension sizes
    std::vector<T> data;  // actual container
    
public:
    //  constructors
    md_vector(size_t d1, size_t d2, size_t d3, size_t d4) : d1(d1), d2(d2), d3(d3), d4(d4), data(d1 * d2 * d3 * d4) {} //  4D
    md_vector(size_t d1, size_t d2, size_t d3) : d1(d1), d2(d2), d3(d3), data(d1 * d2 * d3) {} //  3D
    md_vector(size_t d1, size_t d2) : d1(d1), d2(d2), data(d1 * d2) {} ; //  2D
    
    // operator ()
    T &operator()(size_t i, size_t j, size_t k, size_t l) {
        return data[ i*d2*d3*d4 + j*d3*d4 + k*d4 + l ];
    }
    
    T &operator()(size_t i, size_t j, size_t k) {
        return data [ i*d2*d3 + j*d3 + k ];
    }
    
    T &operator()(size_t i, size_t j) {
        return data [ i*d2 + j ];
    }
};

#endif /* md_vector_h */
