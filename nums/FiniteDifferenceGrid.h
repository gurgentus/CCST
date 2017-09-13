//
// Created by Gurgen Hayrapetyan on 9/9/17.
//

#ifndef NUMS_FINITEDIFFERENCEGRID_H
#define NUMS_FINITEDIFFERENCEGRID_H

#include "Eigen/Dense"
#include <vector>
#include "Node.h"

class FiniteDifferenceGrid {

public:
    std::vector<Node> nodes_;
    FiniteDifferenceGrid(unsigned long num_nodes, double t_min, double t_max);
};

#endif //NUMS_FINITEDIFFERENCEGRID_H
