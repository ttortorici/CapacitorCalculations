#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <Windows.h>
//#include <functions.h>
#include <cmath>

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


PyObject* k_thin(PyObject*, PyObject* args) {
    double g, h, u;

    if (!PyArg_ParseTuple(args, "ddd", &g, &h, &u))
        return NULL;

    double k = f::k_thin(g, h, u);
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


PyObject* ellint_ratio_approx(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = f::ellint_ratio_approx(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* capacitance_bare(PyObject*, PyObject* args) {
    double g, u, h_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "dddidd", &g, &u, &h_s, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_bare(g, u, h_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_bare_k(PyObject*, PyObject* args) {
    double k_a, k_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "ddidd", &k_a, &k_s, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_bare_k(k_a, k_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_geometric(PyObject*, PyObject* args) {
    double g, h_f, u, l;
    int N;

    if (!PyArg_ParseTuple(args, "dddid", &g, &h_f, &u, &N, &l))
        return NULL;

    double cap = f::capacitance_geometric(g, h_f, u, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_geometric_k(PyObject*, PyObject* args) {
    double k_f, l;
    int N;

    if (!PyArg_ParseTuple(args, "did", &k_f, &N, &l))
        return NULL;

    double cap = f::capacitance_geometric_k(k_f, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_total(PyObject*, PyObject* args) {
    double g, h_f, eps_f, u, h_s, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "dddddidd", &g, &h_f, &eps_f, &u, &h_s, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_total(g, h_f, eps_f, u, h_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* capacitance_total_k(PyObject*, PyObject* args) {
    double k_a, k_s, k_f, eps_f, l, eps_s;
    int N;

    if (!PyArg_ParseTuple(args, "ddddidd", &k_a, &k_s, &k_f, &eps_f, &N, &l, &eps_s))
        return NULL;

    double cap = f::capacitance_total_k(k_a, k_s, k_f, eps_f, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* dielectric_constant_relative(PyObject*, PyObject* args) {
    double C_t, C_0, g, h_f, u, l;
    int N;

    if (!PyArg_ParseTuple(args, "dddddid", &C_t, &C_0, &g, &h_f, &u, &N, &l))
        return NULL;

    double eps = f::dielectric_constant_relative(C_t, C_0, g, h_f, u, N, l);
    return PyFloat_FromDouble(eps);
}


PyObject* dielectric_constant_relative_k(PyObject*, PyObject* args) {
    double C_t, C_0, k_f, l;
    int N;

    if (!PyArg_ParseTuple(args, "dddddid", &C_t, &C_0, &k_f, &N, &l))
        return NULL;

    double eps = f::dielectric_constant_relative(C_t, C_0, k_f, N, l);
    return PyFloat_FromDouble(eps);
}


static PyMethodDef idcappy_methods[] = {
    { "k_thin", (PyCFunction)k_thin, METH_VARARGS, R"pbdoc(
        k_thin(g: float, h: float = 500., u: float = 20.) -> float\n\n
        Calculate the modulus for an elliptic integral for small to moderate material thickness.\n
        Can use any units as long as u, g, & h are the same units.\n
        :param g: gap spacing between fingers.\n
        :param h: material thickness.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :return: elliptic integral modulus k.
    )pbdoc" },

    { "k_thick", (PyCFunction)k_thick, METH_VARARGS, R"pbdoc(
        k_thick(g: float, u: float = 20.) -> float\n\n
        Calculate the modulus for an elliptic integral for large material thickness (like air).\n
        Can use any units as long as u & g are the same units.\n
        :param g: gap spacing between fingers.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :return: elliptic integral modulus k.
    )pbdoc" },

    { "ellint_ratio", (PyCFunction)ellint_ratio, METH_O,  R"pbdoc(
        ellint_ratio(k: float) -> float\n\n
        Calculate the ratio K(k)/K'(k) (starts to fail for k < 1e-7).\n
        :param k: the modulus.\n
        :return: the ratio of a complete elliptic integral of the first kind to its complement.
    )pbdoc" },

    { "ellint_ratio_approx", (PyCFunction)ellint_ratio_approx, METH_O,  R"pbdoc(
        ellint_ratio_approx(k: float) -> float\n\n
        Calculate the ratio K(k)/K'(k) for very small k (works well for k < 0.1; essential for k < 1e-7).\n
        :param k: the modulus.\n
        :return: the ratio of a complete elliptic integral of the first kind to its complement for very small k.
    )pbdoc" },

    { "capacitance_bare", (PyCFunction)capacitance_bare, METH_VARARGS, R"pbdoc(
        capacitance_bare(g: float, u: float = 20., h_sub: float = 500., N: int = 50, l: float = 1., eps_s: float = [for silia]) -> float\n\n
        Calculate the bare capacitance of an interdigital capacitor in pF based on geometry.\n
        u, g, & h_s must have the same units, but l must be mm.\n
        :param g: gap spacing between fingers.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :param h_s: substrate thickness.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :param eps_s: relative dielectric constant of the substrate.\n
        :return: bare capacitance in pF.
    )pbdoc" },

    { "capacitance_bare_k", (PyCFunction)capacitance_bare_k, METH_VARARGS, R"pbdoc(
        capacitance_bare(k_a: float, k_s: float, N: int = 50, l: float = 1., eps_s: float = [for silica]) -> float\n\n
        Calculate the bare capacitance of an interdigital capacitor in pF based on moduli.\n
        :param k_a: modulus for air contribution.\n
        :param k_s: modulus for substrate contribution.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :param eps_s: relative dielectric constant of the substrate.\n
        :return: bare capacitance in pF.
    )pbdoc" },

    { "capacitance_geometric", (PyCFunction)capacitance_geometric, METH_VARARGS, R"pbdoc(
        capacitance_bare(g: float, h_f: float, u: float = 20., N: int = 50, l: float = 1.) -> float\n\n
        Calculate the geometric capacitance for a film in pF based on geometry.\n
        u, g, & h_f must have the same units, but l must be mm.\n
        :param g: gap spacing between fingers.\n
        :param h_f: film thickness.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :return: geometric capacitance in pF.
    )pbdoc" },

    { "capacitance_geometric_k", (PyCFunction)capacitance_geometric_k, METH_VARARGS, R"pbdoc(
        capacitance_bare(k_f: float, N: int = 50, l: float = 1.) -> float\n\n
        Calculate the geometric capacitance for a film in pF based on the modulus.\n
        :param k_f: modulus for film contribution.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :return: geometric capacitance in pF.
    )pbdoc" },

    { "capacitance_total", (PyCFunction)capacitance_total, METH_VARARGS, R"pbdoc(
        capacitance_bare(g: float, h_f: float, eps_f: float, u: float = 20., h_s: float = 500., N: int = 50, l: float = 1., eps_s: float = [for silica]) -> float\n\n
        Calculate the bare capacitance of an interdigital capacitor in pF based on geometry.\n
        u, g, h_s, & h_f must have the same units, but l must be mm.\n
        :param g: gap spacing between fingers.\n
        :param h_f: film thickness.\n
        :param eps_f: relative dielectric constant of film.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :param h_s: substrate thickness.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :param eps_s: relative dielectric constant of substrate.\n
        :return: total capacitance in pF.
    )pbdoc" },

    { "capacitance_total_k", (PyCFunction)capacitance_total_k, METH_VARARGS, R"pbdoc(
        capacitance_bare(k_a: float, k_s: float, k_f: float, eps_f, N: int = 50, L: float = 1., eps_s: float = [for silica]) -> float\n\n
        Calculate the bare capacitance of an interdigital capacitor in pF based on the moduli.\n
        :param k_a: modulus for air contribution.\n
        :param k_s: modulus for substrate contribution.\n
        :param k_f: modulus for film contribution.\n
        :param eps_f: relative dielectric constant of film.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :param eps_s: relative dielectric constant of substrate.\n
        :return: total capacitance in pF.
    )pbdoc" },

    { "dielectric_constant_relative", (PyCFunction)dielectric_constant_relative, METH_VARARGS, R"pbdoc(
        dielectric_constant_relative(delC: float, g: float, h_f: float, u: float = 20., N: int = 50, l: float = 1.) -> float\n\n
        Calculate the relative dielectric constant given the geometry\n
        u, g, & h_f must have the same units, but l must be mm.\n
        :param delC: total C - bare C,\n
        :param g: gap spacing between fingers.\n
        :param h_f: film thickness.\n
        :param u: unit cell of interdigital capacitor (finger width + gap spacing).\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :return: relative dielectric constant of substrate.
    )pbdoc" },

    { "dielectric_constant_relative_k", (PyCFunction)dielectric_constant_relative_k, METH_VARARGS, R"pbdoc(
        dielectric_constant_relative(delC: float, k_f, N: int = 50, l: float = 1.) -> float\n\n
        Calculate the relative dielectric constant given the modulus\n
        :param delC: total C - bare C
        :param k_f: modulus for film contribution.\n
        :param N: total number of fingers.\n
        :param l: finger length in mm.\n
        :retrun: relative dielectric constant of substrate.
    )pbdoc" },


    //{ "find_gap", (PyCFunction)find_gap, METH_VARARGS, "find_gap(capacitance: float, u: float, t_sub: float, N: int, L: float, eps_sub: float) -> float\n\n:param capacitance: in pF;\n:param u: unit cell size in um;\n:param N: number of fingers;\n:param L: finger length in um;\n:param eps_sub: dielectric constant of substrate;\n:param t_sub: thickness of substrate"},

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