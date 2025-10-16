#include <vector>
#include <memory>
#include <iostream>
#include "shape.h"

int main(){
    std::vector<std::shared_ptr<Shape>> shapes_vector;
    // Do not use std::shared_ptr....; use make_shared instead

    shapes_vector.push_back(std::make_shared<Circle>(1.0));
    shapes_vector.push_back(std::make_shared<Rectangle>(3., 1.0));

    // By creating a constant reference to this type, we avoid 
    // using more memory than creating useless copies.
    // In case we have to modify, do not insert const.
    for (const auto & s: shapes_vector){
        std::cout << s -> getArea() << ", " << s -> getName() << std::endl;
        // std::cerr flushed the buffer automatically, if needed
    }

    return 0;
}