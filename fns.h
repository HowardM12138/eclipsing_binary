#ifndef FNS_H
#define FNS_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <tgmath.h>
#include <cmath>
#include <windows.h>
#include <fstream>

/*
 *mass is in meter
 *radius is in meter
 *temperature is in celcius
 *brightness is in lux
 *x is in meter
 *y is in meter
 *z is in meter
 *separation is in meter
*/

const double PI = 3.14159265358979323846;
const double STEFAN_BOLTZMANN_CONSTANT = 5.670367 / pow(10, 8); //Units kg * s^-3 * K^-4
const double SUN_MASS = 1.989 * pow(10, 30);
const double SUN_RADIUS = 695.51 * pow(10, 6);

struct Planet {
    std::string name;
    int index;
    double mass;
    double r;
    double area;
    double temperature;
    double luminosity;
    double x;
    double y;
    double z;
};

void init_planet (std::string name, double mass, double radius, double brightness);
void start(double s, int ppr, double i);
void run();
void calculate_xyz (Planet &planet);
double get_luminosity();
void print_planet(Planet &planet);
void print_coord (Planet &planet);

#endif
