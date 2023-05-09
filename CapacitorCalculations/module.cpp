#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <Windows.h>
//#include <functions.h>
#include <cmath>

double pi = 3.141592653589793115997963468544185161590576171875;
double ln2 = 0.69314718055994530941723212145817656807550013436;

namespace f
{
    double eps_0 = 0.00885418781762;  // pF/mm

<<<<<<< HEAD
=======
    /**
    * Calculate the modulus of an elliptic integral for interdigital capacitor calculation
    * All parameters must be in the same length units
    * @param w_g: gap distance
    * @param w_u: unit cell size
    * @param h: material thickness
    * @return: the modulus
    */
    double k_thin(double w_g, double w_u, float h) {
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
    double k_thick(double w_g, double w_u) {
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
    * Calculate the ratio K'(k)/K(k)
    * @param k: the modulus
    * @return: the ratio of the integral to its compliment
    */
    double ellint_inv_ratio(double k) {
        double k_prime = std::sqrt(1 - k * k);
        return std::comp_ellint_1(k_prime) / std::comp_ellint_1(k);
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
    double capacitance_bare(double w_g, double w_u, double sub_thickness, double N, double l, double eps_sub) {
        double k_air = k_thick(w_u, w_g);
        double k_sub = k_thin(w_u, w_g, sub_thickness);
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub));
    };

    /**
    * Calculate the capacitance of a bare interdigital capacitor using the modulus
    * @param k_air: modulus for air contribution
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @param k_sub: modulus for substrate contribution
    * @return: the capacitance in pF
    */
    double capacitance_bare_k(double k_air, double N, double l, double eps_sub, double k_sub) {
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub));
    };

    /**
    * Calculate the geometric capacitance of a capacitor with a given film thickness
    * @param w_g: gap distance (must be same units as w_u and thickness)
    * @param w_u: unit cell size (must be same units as w_g and thickness)
    * @param thickness: thickness of the film (must be same units as w_u and w_g)
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @return: the capacitance in pF
    */
    double capacitance_geometric(double w_g, double w_u, double thickness, double N, double l) {
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
    double capacitance_geometric_k(double k_film, double N, double l) {
        return eps_0 * (N - 1) * l * ellint_ratio_approx(k_film);
    }

    /**
    * Calculate the capacitance of an interdigital capacitor with a film using its geometry
    * @param w_g: gap distance (must be same units as w_u, sub_thickness, and film_thickness)
    * @param w_u: unit cell size (must be same units as w_g, sub_thickness, and film_thickness)
    * @param sub_thickness: thickness of substrate (must be same units as w_u, w_g, and film_thickness)
    * @param film_thickness: thickness of the film (must be same units as w_u, w_g, and sub_thickness)
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @param eps_film: relative dielectric constant of film
    * @return: the capacitance in pF
    */
    double capacitance_total(double w_g, double w_u, double sub_thickness,
        double film_thickness, double N, double l, double eps_sub, double eps_film) {
        double k_air = k_thick(w_u, w_g);
        double k_sub = k_thin(w_u, w_g, sub_thickness);
        double k_film = k_thin(w_u, w_g, film_thickness);
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1)
               * ellint_ratio(k_sub) + (eps_film - 1) * ellint_ratio(k_film));
    }


    /**
    * Calculate the capacitance of an interdigital capacitor with a film using the modulii
    * @param k_air: modulus for air contribution
    * @param k_sub: modulus for substrate contribution
    * @param k_film: modulus for film contribution
    * @param N: total number of fingers
    * @param l: finger length in mm
    * @param eps_sub: relative dielectric constant of substrate
    * @param eps_film: relative dielectric constant of film
    * @return: the capacitance in pF
    */
    double capacitance_total_k(double k_air, double k_sub, double k_film, double N, double l, double eps_sub, double eps_film) {
        return eps_0 * (N - 1) * l * (2 * ellint_ratio(k_air) + (eps_sub - 1) * ellint_ratio(k_sub) + (eps_film - 1) * ellint_ratio(k_film));
    }


    /**
    * Calculate the relative dielectric constant given the geometry
    * @param C_t: total capacitance in pF
    * @param C_0: bare capacitance in pF
    * @param w_g: gap distance (must be same units as w_u, sub_thickness, and film_thickness)
    * @param w_u: unit cell size (must be same units as w_g, sub_thickness, and film_thickness)
    * @param thickness: thickness of the film (must be same units as w_u, w_g, and sub_thickness)
    * @param N: total number of fingers
    * @param l: finger length in mm
    */
    double dielectric_constant_relative(double C_t, double C_0, double w_g, double w_u, double thickness, double N, double l) {
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
    double dielectric_constant_relative_k(double C_t, double C_0, double k_film, double N, double l) {
        double denominator = eps_0 * (N - 1) * l;
        return (C_t - C_0) * ellint_inv_ratio_approx(k_film) / denominator + 1;
    }

}
>>>>>>> parent of 62ca22d (1.3)


PyObject* k_thin(PyObject*, PyObject* args) {
    double g, h, u;

    if (!PyArg_ParseTuple(args, "ddd", &g, &u, &h))
        return NULL;

    double k = f::k_thin(g, u, h);
    return PyFloat_FromDouble(k);
};


PyObject* k_thick(PyObject*, PyObject* args) {
    double g, u;

    if (!PyArg_ParseTuple(args, "dd", &g, &u))
        return NULL;

    double k = f::k_thick(g, u);
    return PyFloat_FromDouble(k);
};


PyObject* ellint_ratio(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = f::ellint_ratio(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* ellint_inv_ratio(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = f::ellint_inv_ratio(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* ellint_ratio_approx(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = f::ellint_ratio_approx(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* capacitance_bare(PyObject*, PyObject* args) {
    double g, u, h_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "dddddd", &g, &u, &h_s, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_bare(g, u, h_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_bare_k(PyObject*, PyObject* args) {
    double k_a, k_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "ddddd", &k_a, &k_s, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_bare_k(k_a, k_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_geometric(PyObject*, PyObject* args) {
    double g, h_f, u, l;
    int N;

    if (!PyArg_ParseTuple(args, "ddddd", &g, &u, &h_f, &N, &l))
        return NULL;

    double cap = f::capacitance_geometric(g, u, h_f, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_geometric_k(PyObject*, PyObject* args) {
    double k_f, l;
    int N;

    if (!PyArg_ParseTuple(args, "ddd", &k_f, &N, &l))
        return NULL;

    double cap = f::capacitance_geometric_k(k_f, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_total(PyObject*, PyObject* args) {
    double g, h_f, eps_f, u, h_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "dddddddd", &g, &u, &h_s, &h_f, &N, &l, &eps_s, &eps_f))
        return NULL;

    double cap = f::capacitance_total(g, u, h_s, h_f, N, l, eps_s, eps_f);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_total_k(PyObject*, PyObject* args) {
    double k_a, k_s, k_f, eps_f, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "ddddddd", &k_a, &k_s, &k_f, &N, &l, &eps_s, &eps_f))
        return NULL;

    double cap = f::capacitance_total_k(k_a, k_s, k_f, N, l, eps_s, eps_f);
    return PyFloat_FromDouble(cap);
}


PyObject* dielectric_constant_relative(PyObject*, PyObject* args) {
    double C_t, C_0, g, h_f, u, l;
    int N;

    if (!PyArg_ParseTuple(args, "ddddddd", &C_t, &C_0, &g, &u, &h_f, &N, &l))
        return NULL;

    double eps = f::dielectric_constant_relative(C_t, C_0, g, u, h_f, N, l);
    return PyFloat_FromDouble(eps);
}


PyObject* dielectric_constant_relative_k(PyObject*, PyObject* args) {
    double C_t, C_0, k_f, l;
    int N;

    if (!PyArg_ParseTuple(args, "ddddd", &C_t, &C_0, &k_f, &N, &l))
        return NULL;

    double eps = f::dielectric_constant_relative_k(C_t, C_0, k_f, N, l);
    return PyFloat_FromDouble(eps);
}


static PyMethodDef idcappy_methods[] = {
    { "k_thin", (PyCFunction)k_thin, METH_VARARGS, R"pbdoc(
        k_thin(g: float, u: float, h: float) -> float
        
        Calculate the modulus for an elliptic integral for small to moderate material thickness.
        Can use any units as long as u, g, & h are the same units.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :param h: material thickness.
        :return: elliptic integral modulus k.
    )pbdoc" },

    { "k_thick", (PyCFunction)k_thick, METH_VARARGS, R"pbdoc(
        k_thick(g: float, u: float) -> float
        
        Calculate the modulus for an elliptic integral for large material thickness (like air).
        Can use any units as long as u & g are the same units.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :return: elliptic integral modulus k.
    )pbdoc" },

    { "ellint_ratio", (PyCFunction)ellint_ratio, METH_O,  R"pbdoc(
        ellint_ratio(k: float) -> float

        Calculate the ratio K(k)/K'(k) (starts to fail for k < 1e-7).
        :param k: the modulus.
        :return: the ratio of a complete elliptic integral of the first kind to its complement.
    )pbdoc" },

    { "ellint_inv_ratio", (PyCFunction)ellint_inv_ratio, METH_O,  R"pbdoc(
        ellint_inv_ratio(k: float) -> float

        Calculate the ratio K'(k)/K(k).
        :param k: the modulus.
        :return: the ratio of the complement to the complete elliptic integral of the first kind to itself.
    )pbdoc" },

    { "ellint_ratio_approx", (PyCFunction)ellint_ratio_approx, METH_O,  R"pbdoc(
        ellint_ratio_approx(k: float) -> float

        Calculate the ratio K(k)/K'(k) for very small k (works well for k < 0.1; essential for k < 1e-7).
        :param k: the modulus.
        :return: the ratio of a complete elliptic integral of the first kind to its complement for very small k.
    )pbdoc" },

    { "capacitance_bare", (PyCFunction)capacitance_bare, METH_VARARGS, R"pbdoc(
        capacitance_bare(g: float, u: float, h_s: float, N: float, l: float, eps_s: float) -> float

        Calculate the bare capacitance of an interdigital capacitor in pF based on geometry.
        u, g, & h_s must have the same units, but l must be mm.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :param h_s: substrate thickness.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :param eps_s: relative dielectric constant of the substrate.
        :return: bare capacitance in pF.
    )pbdoc" },

    { "capacitance_bare_k", (PyCFunction)capacitance_bare_k, METH_VARARGS, R"pbdoc(
        capacitance_bare(k_a: float, k_s: float, N: float, l: float, eps_s: float) -> float

        Calculate the bare capacitance of an interdigital capacitor in pF based on moduli.
        :param k_a: modulus for air contribution.
        :param k_s: modulus for substrate contribution.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :param eps_s: relative dielectric constant of the substrate.
        :return: bare capacitance in pF.
    )pbdoc" },

    { "capacitance_geometric", (PyCFunction)capacitance_geometric, METH_VARARGS, R"pbdoc(
        capacitance_geometric(g: float, u: float, h_f: float, N: float, l: float) -> float

        Calculate the geometric capacitance for a film in pF based on geometry.
        u, g, & h_f must have the same units, but l must be mm.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :param h_f: film thickness.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :return: geometric capacitance in pF.
    )pbdoc" },

    { "capacitance_geometric_k", (PyCFunction)capacitance_geometric_k, METH_VARARGS, R"pbdoc(
        capacitance_geometric(k_f: float, N: float, l: float) -> float

        Calculate the geometric capacitance for a film in pF based on the modulus.
        :param k_f: modulus for film contribution.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :return: geometric capacitance in pF.
    )pbdoc" },

    { "capacitance_total", (PyCFunction)capacitance_total, METH_VARARGS, R"pbdoc(
        capacitance_total(g: float, u: float, h_s: float, h_f: float, N: float, l: float, eps_s: float, eps_f: float) -> float

        Calculate the bare capacitance of an interdigital capacitor in pF based on geometry.
        u, g, h_s, & h_f must have the same units, but l must be mm.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :param h_s: substrate thickness.
        :param h_f: film thickness.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :param eps_s: relative dielectric constant of substrate.
        :param eps_f: relative dielectric constant of film.
        :return: total capacitance in pF.
    )pbdoc" },

    { "capacitance_total_k", (PyCFunction)capacitance_total_k, METH_VARARGS, R"pbdoc(
        capacitance_total(k_a: float, k_s: float, k_f: float, N: float, L: float, eps_s: float, eps_f: float) -> float

        Calculate the bare capacitance of an interdigital capacitor in pF based on the moduli.
        :param k_a: modulus for air contribution.
        :param k_s: modulus for substrate contribution.
        :param k_f: modulus for film contribution.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :param eps_s: relative dielectric constant of substrate.
        :param eps_f: relative dielectric constant of film.
        :return: total capacitance in pF.
    )pbdoc" },

    { "dielectric_constant_relative", (PyCFunction)dielectric_constant_relative, METH_VARARGS, R"pbdoc(
        dielectric_constant_relative(C_t: float, C_0: float, g: float, u: float, h_f: float, N: float, l: float) -> float

        Calculate the relative dielectric constant given the geometry
        u, g, & h_f must have the same units, but l must be mm.
        :param C_t: measured capacitance.
        :param C_0: measured capacitance before film was grown.
        :param g: gap spacing between fingers.
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).
        :param h_f: film thickness.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :return: relative dielectric constant of substrate.
    )pbdoc" },

    { "dielectric_constant_relative_k", (PyCFunction)dielectric_constant_relative_k, METH_VARARGS, R"pbdoc(
        dielectric_constant_relative(C_t: float, C_0: float, k_f: float, N: float, l: float) -> float

        Calculate the relative dielectric constant given the modulus
        :param C_t: measured capacitance.
        :param C_0: measured capacitance before film was grown.
        :param k_f: modulus for film contribution.
        :param N: total number of fingers.
        :param l: finger length in mm.
        :retrun: relative dielectric constant of substrate.
    )pbdoc" },

    // Terminate the array with an object containing nulls.
    { nullptr, nullptr, 0, nullptr }
};


static PyModuleDef idcappy_module = {
    PyModuleDef_HEAD_INIT,
    "idcappy",                                  // Module name to use with Python import statements
    "Module Description",                       // Module description
    0,
    idcappy_methods                             // Structure that defines the methods of the module
};


PyMODINIT_FUNC PyInit_idcappy(void) {
    return PyModule_Create(&idcappy_module);
}