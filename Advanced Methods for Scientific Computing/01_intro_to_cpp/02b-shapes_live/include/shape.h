// You should avoid entering multiple times the
// same libraries across main files and headers, as it throws a compile error.
// Thus, use ifdef instead

#ifndef SHAPES_H
#define SHAPES_H

#include <cmath>

class Shape{
    // The default of class is private
    public:
        // Constructor and destructors
        Shape () = default;
        virtual ~Shape() = default;     // recall that virtual means that it can be overridden (?? CHECK THIS)

        // Method declarations
        virtual double getArea() = 0;
        constexpr virtual const char *getName() = 0;

};

class Circle: public Shape{
    private:
        const double radius;
    public:
        Circle(const double r);
        virtual ~Circle() override = default;
        // Inline definition
        double getArea() override {return (M_PI*radius*radius);};
        constexpr virtual const char *getName() override {return "Circle";};
};

class Rectangle: public Shape{
    private:
        const double height;
        const double basis;

    public:
        Rectangle(const double, const double);
        virtual ~Rectangle() override = default;
        
        double getArea() override {return (height * basis);};
};

#endif