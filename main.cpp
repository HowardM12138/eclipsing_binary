#include <iostream>
#include "fns.h"

int main() {

    //Parameter: name, mass ratio to sun's mass, radius ratio to sun, temperature in kelvin
    init_planet("SZ Her #2", 0.4072, 2.5, 6000);
    init_planet("SZ Her #1", 1, 2.5, 4000);
    //Parameter: Separation ratio to sun's radius, phase per revolution, inclination in degree
    start(5, 10000, 87.57);

    return 0;
}
