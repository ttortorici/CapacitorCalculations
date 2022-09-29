#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
//#include <iostream>

double pi = 3.141592653589793115997963468544185161590576171875;
double ln2 = 0.69314718055994530941723212145817656807550013436;

namespace f
{

    double eps_0 = 0.00885418781762e-3;  // pF/mm
    double eps_silica = 3.9;


    /**
    * Calculate the modulus of an elliptic integral for interdigital capacitor calculation
    * All parameters must be in the same length units
    * @param w_g: gap distance
    * @param h: material thickness
    * @param w_u: unit cell size
    * @return: the modulus
    */
    double k_thin(double w_g, float h = 500., double w_u = 20.) {
        double invh = pi / (4 * h);
        double x = std::sinh((3 * w_u - w_g) * invh);
        double y = std::sinh((w_u + w_g) * invh);
        double z = std::sinh((w_u - w_g) * invh);
        return z / y * sqrt((x - y / x * y) / (x - z / x * z));
    };


    /**
    * Calculate the modulus of an elliptic integral for interdigital capacitor calculation for thick material like air
    * All parameters must be in the same length units
    * @param w_g: gap distance
    * @param w_u: unit cell size
    * @return: the modulus
    */
    double k_thick(double w_g, double w_u = 20.) {
        // approximates k for large h
        double diff = w_u - w_g;
        return diff / (w_u + w_g) * std::sqrt(2 * diff / (2 * w_u - w_g));
    };

    /**
    * Calculate the ratio K(k)/K'(k)
    * @param k: the modulus
    * @return: the ratio of the integral to its compliment
    */
    double ellint_ratio(double k) {
        // calculates elliptic integral divided by compliment
        double k_prime = std::sqrt(1 - k * k);
        return std::comp_ellint_1(k) / std::comp_ellint_1(k_prime);
    };


    /**
    * Calculate the ratio K(k)/K'(k) for k < 1e-7 (good for k < 0.1)
    * @param k: the modulus
    * @return: the ratio of the integral to its compliment
    */
    double ellint_ratio_approx(double k) {
        //double k_part = std::pow(1 - k * k, 0.25);
        double k2 = k * k;
        double ks = k2;
        double terms[3] = { 0.9375, 0.546875, 0.03759765625 };
        double k_part = 0.25 * ks;
        for (int ii = 0; ii < 3; ii++) {
            ks *= k2;
            k_part += terms[ii] * ks;
        }
        return pi / std::log(4 / k_part - 2);
    };

    /**
    * Calculate the ratio K'(k)/K(k) for k < 0.1
    * @param k: the modulus
    * @return: the ratio of the integral's compliment to itself
    */
    double ellint_inv_ratio_approx(double k) {
        double k2 = k * k;
        double ks = k2;
        double terms[3] = { 0.9375, 0.546875, 0.03759765625 };
        double k_part = 0.25 * ks;
        for (int ii = 0; ii < 3; ii++) {
            ks *= k2;
            k_part += terms[ii] * ks;
        }
        return std::log(4 / k_part - 2) / pi;
    };


    /**
    * Calculate the capacitance of a bare interdigital capacitor using its geometry
    * @param w_g: gap distance (must be same units as w_u and sub_thickness)
    * @param w_u: unit cell size (must be same units as w_g and sub_thickness)
    * @param sub_thickness: thickness of substrate (must be same units as w_u and w_g)
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @return: the capacitance in pF
    */
    double capacitance_bare(double w_g, double w_u = 20., double sub_thickness = 500., int N = 50, double l = 1., double eps_sub = eps_silica) {
        double k_air = k_thick(w_u, w_g);
        double k_sub = k_thin(w_u, w_g, sub_thickness);
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub));
    };


    /**
    * Calculate the capacitance of a bare interdigital capacitor using the modulus
    * @param k_air: modulus for air contribution
    * @param k_sub: modulus for substrate contribution
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @return: the capacitance in pF
    */
    double capacitance_bare_k(double k_air, double k_sub, int N = 50, double l = 1., double eps_sub = eps_silica) {
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub));
    };


    /**
    * Calculate the geometric capacitance of a capacitor with a given film thickness
    * @param w_g: gap distance (must be same units as w_u and thickness)
    * @param thickness: thickness of the film (must be same units as w_u and w_g)
    * @param w_u: unit cell size (must be same units as w_g and thickness)
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @return: the capacitance in pF
    */
    double capacitance_geometric(double w_g, double thickness, double w_u = 20., int N = 50, double l = 1.) {
        double k_film = k_thin(w_u, w_g, thickness);
        return eps_0 * (N - 1) * l * ellint_ratio_approx(k_film);
    }


    /**
    * Calculate the geometric capacitance of a capacitor with a given film modulus
    * @param k_film: modulus for film contribution
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @return: the capacitance in pF
    */
    double capacitance_geometric_k(double k_film, int N = 50, double l = 1.) {
        return eps_0 * (N - 1) * l * ellint_ratio_approx(k_film);
    }


    /**
    * Calculate the capacitance of an interdigital capacitor with a film using its geometry
    * @param w_g: gap distance (must be same units as w_u, sub_thickness, and film_thickness)
    * @param film_thickness: thickness of the film (must be same units as w_u, w_g, and sub_thickness)
    * @param w_u: unit cell size (must be same units as w_g, sub_thickness, and film_thickness)
    * @param sub_thickness: thickness of substrate (must be same units as w_u, w_g, and film_thickness)
    * @param eps_film: relative dielectric constant of film
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @return: the capacitance in pF
    */
    double capacitance_total(double w_g, double film_thickness, double eps_film, double w_u = 20.,
        double sub_thickness = 500., int N = 50, double l = 1., double eps_sub = eps_silica) {
        double k_air = k_thick(w_u, w_g);
        double k_sub = k_thin(w_u, w_g, sub_thickness);
        double k_film = k_thin(w_u, w_g, film_thickness);
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub) + (eps_film - 1) * ellint_ratio(k_film));
    }


    /**
    * Calculate the capacitance of an interdigital capacitor with a film using the modulii
    * @param k_air: modulus for air contribution
    * @param k_sub: modulus for substrate contribution
    * @param k_film: modulus for film contribution
    * @param eps_film: relative dielectric constant of film
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @return: the capacitance in pF
    */
    double capacitance_total_k(double k_air, double k_sub, double k_film, double eps_film, int N = 50, double l = 1., double eps_sub = eps_silica) {
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub) + (eps_film - 1) * ellint_ratio(k_film));
    }


    /**
    * Calculate the relative dielectric constant given the geometry
    * @param C_t: total capacitance in pF
    * @param C_0: bare capacitance in pF
    * @param w_g: gap distance (must be same units as w_u, sub_thickness, and film_thickness)
    * @param thickness: thickness of the film (must be same units as w_u, w_g, and sub_thickness)
    * @param w_u: unit cell size (must be same units as w_g, sub_thickness, and film_thickness)
    * @param N: total number of fingers
    * @param l: finger length in mm
    */
    double dielectric_constant_relative(double C_t, double C_0, double w_g, double thickness, double w_u = 20., int N = 50, double l = 1.) {
        double k_film = k_thin(w_u, w_g, thickness);
        double denominator = eps_0 * (N - 1) * l;
        return (C_t - C_0) * ellint_inv_ratio_approx(k_film) / denominator + 1;
    }


    /**
    * Calculate the relative dielectric constant given the modulus
    * @param C_t: total capacitance in pF
    * @param C_0: bare capacitance in pF
    * @param k_film: modulus for film contribution
    * @param N: total number of fingers
    * @param l: finger length in mm
    */
    double dielectric_constant_relative_k(double C_t, double C_0, double k_film, int N = 50, double l = 1.) {
        double denominator = eps_0 * (N - 1) * l;
        return (C_t - C_0) * ellint_inv_ratio_approx(k_film) / denominator + 1;
    }

}
#endif // FUNCTIONS_H