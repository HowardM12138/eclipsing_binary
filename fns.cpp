#include "fns.h"

double theta = PI; //in radian

double separation;
int phase_per_rev;
double inclination;

double q;

Planet planet_1;
Planet planet_2;

void start(double s, int ppr, double i) {
    if (planet_1.mass != 0 && planet_2.mass != 0) {
        separation = s * SUN_RADIUS;
        phase_per_rev = ppr;
        inclination = i * PI / 180;
        q = planet_2.mass / planet_1.mass;
        run();
    }
}

void init_planet (std::string name, double mass, double radius, double temperature) {
    if (planet_1.mass == 0) {
        planet_1.name = name;
        planet_1.index = 1;
        planet_1.mass = mass * SUN_MASS;
        planet_1.r = radius * SUN_RADIUS;
        planet_1.area = PI * pow(planet_1.r, 2);
        planet_1.temperature = temperature;
        planet_1.luminosity = STEFAN_BOLTZMANN_CONSTANT * planet_1.area * 4 * pow(temperature, 4);
        //std::cout << "init planet_1\n";
    } else if (planet_2.mass == 0) {
        planet_2.name = name;
        planet_2.index = 2;
        planet_2.mass = mass * SUN_MASS;
        planet_2.r = radius * SUN_RADIUS;
        planet_2.area = PI * pow(planet_2.r, 2);
        planet_2.temperature = temperature;
        planet_2.luminosity = STEFAN_BOLTZMANN_CONSTANT * planet_2.area * 4 * pow(temperature, 4);
        //std::cout << "init planet_2\n\n";
    } else {
        std::cout << "faled to init planet\n";
    }
}

void run() {
    int i = 0;
    std::ofstream data_file ("data.txt");
    if (!data_file.is_open()) std::cout << "Unable to create data file.\n";
    else {
        while (i <= phase_per_rev * 2) {
            calculate_xyz(planet_1);
            calculate_xyz(planet_2);
            /*
            std::cout << "Phase " << i+1 << ": \n";
            print_coord(planet_1);
            print_coord(planet_2);
            */
            double luminosity = get_luminosity();
            data_file << round((theta - PI) / (2*PI) * 100000) / 100000 << "  " << luminosity << "\n";
            /*
            std::cout << "luminosity: " << luminosity << "\n";
            std::cout << "\n";
            */
            //Sleep(100);
            theta = theta + (2*PI / phase_per_rev);
            i++;
        }
    }
    data_file.close();
}

void calculate_xyz (Planet &planet) {
    double x = separation * sin(theta);
    planet.x = (planet.index == 1) ? -1 * x / (1 + 1/q) : x / (1 + q);
    planet.x = round(planet.x * 100000) / 100000;
    if (planet.x == 0.0) planet.x = 0.0;
    double y = separation * cos(inclination) * cos(theta);
    planet.y = (planet.index == 1) ? -1 * y / (1 + 1/q) : y / (1 + q);
    planet.y = round(planet.y * 100000) / 100000;
    if (planet.y == 0.0) planet.y = 0.0;
    double z = separation * sin(inclination) * cos(theta);
    planet.z = (planet.index == 1) ? -1 * z / (1 + 1/q) : z / (1 + q);
    planet.z = round(planet.z * 100000) / 100000;
    if (planet.z == 0.0) planet.z = 0.0;
}

double get_luminosity() {
    double p = sqrt(pow(planet_2.x - planet_1.x, 2) + pow(planet_2.y - planet_1.y, 2));

    double r1 = planet_1.r;
    double r2 = planet_2.r;

    double theta_1 = acos((pow(r2, 2) - pow(r1, 2) - pow(p, 2)) / (-2 * r1 * p)) * 2;
    if (theta_1 > PI) theta_1 = 2*PI - theta_1;
    double theta_2 = acos((pow(r1, 2) - pow(r2, 2) - pow(p, 2)) / (-2 * r2 * p)) * 2;
    if (theta_2 > PI) theta_2 = 2*PI - theta_2;

    double arca_1 = 0.5 * pow(r1, 2) * (theta_1 - sin(theta_1));
    double arca_2 = 0.5 * pow(r2, 2) * (theta_2 - sin(theta_2));

    double a1 = planet_1.area;
    double a2 = planet_2.area;

    double luminosity_1 = 0.0, luminosity_2 = 0.0;

    if (p >= (r1 + r2)) {
        //no eclipse
        //std::cout << "no eclipsing\n";
        luminosity_1 = planet_1.luminosity;
        luminosity_2 = planet_2.luminosity;
    } else if (p < (r1+r2) && p >= sqrt(abs(pow(r1, 2) - pow(r2, 2)))) {
        //shallow eclipse
        //std::cout << "shallow eclipsing\n";
        if (planet_1.z > planet_2.z) {
            luminosity_1 = planet_1.luminosity;
            luminosity_2 = planet_2.luminosity * (a2 - arca_1 - arca_2) / a2;
        } else {
            luminosity_1 = planet_1.luminosity * (a1 - arca_1 - arca_2) / a1;
            luminosity_2 = planet_2.luminosity;
        }
    } else if (p < sqrt(abs(pow(r1 ,2) - pow(r2, 2))) && p > abs(r1 - r2)) {
        //deep eclipse
        //std::cout << "deep eclipsing\n";
        if (planet_1.z > planet_2.z) {
            luminosity_1 = planet_1.luminosity;
            if (r1 >= r2) luminosity_2 = planet_2.luminosity * (arca_2 - arca_1) / a2;
            else luminosity_2 = planet_2.luminosity * (a2 - a1 + arca_1 - arca_2) / a2;
        } else {
            if (r1 >= r2) luminosity_1 = planet_1.luminosity * (a1 - a2 + arca_2 - arca_1) / a1;
            else luminosity_1 = planet_1.luminosity * (arca_1 - arca_2) / a1;
            luminosity_2 = planet_2.luminosity;
        }
    } else if (p <= abs(r1 - r2)){
        //total eclipse
        //std::cout << "total eclipsing\n";
        if (planet_1.z > planet_2.z) {
            luminosity_1 = planet_1.luminosity;
            if (r1 >= r2) luminosity_2 = 0;
            else luminosity_2 = planet_2.luminosity * (a2 - a1) / a2;
        } else {
            if (r1 >= r2) luminosity_1 = planet_1.luminosity * (a1 - a2) / a1;
            else luminosity_1 = 0;
            luminosity_2 = planet_2.luminosity;
        }
    }

    return round((luminosity_1 + luminosity_2) * 10000) / 10000;
}

void print_planet (Planet &planet) {
    std::cout << "\nName: " << planet.name << "\n";
    std::cout << "Mass: " << planet.mass << "\n";
    std::cout << "Radius: " << planet.r << "\n";
    std::cout << "Area: " << planet.area << "\n";
    std::cout << "Temperature: " << planet.temperature << "\n";
    std::cout << "Luminosity: " << planet.luminosity << "\n\n";
}

void print_coord (Planet &planet) {
    std::cout << planet.name << " ";
    std::cout << "X: " << planet.x << " ";
    std::cout << "Y: " << planet.y << " ";
    std::cout << "Z: " << planet.z << " \n";
}
