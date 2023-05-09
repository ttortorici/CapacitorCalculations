#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <Windows.h>
#include <cmath>
//#include <iostream>

double pi = 3.141592653589793115997963468544185161590576171875;
double ln2 = 0.69314718055994530941723212145817656807550013436;
double eps_0 = 0.00885418781762;  // pF/mm





PyObject* Pyk_thin(PyObject*, PyObject* args) {
    double g, h, u;

    if (!PyArg_ParseTuple(args, "ddd", &g, &u, &h))
        return NULL;

    double k = k_thin(g, u, h);
    return PyFloat_FromDouble(k);
};


PyObject* Pyk_thick(PyObject*, PyObject* args) {
    double g, u;

    if (!PyArg_ParseTuple(args, "dd", &g, &u))
        return NULL;

    double k = k_thick(g, u);
    return PyFloat_FromDouble(k);
};


PyObject* Pyellint_ratio(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = ellint_ratio(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* Pyellint_inv_ratio(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = ellint_inv_ratio(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* Pyellint_ratio_approx(PyObject*, PyObject* modulus) {
    double k = PyFloat_AsDouble(modulus);
    double ratio = ellint_ratio_approx(k);
    return PyFloat_FromDouble(ratio);
};


PyObject* Pycapacitance_bare(PyObject*, PyObject* args) {
    double g, u, h_s, l, eps_s, N;

    if (!PyArg_ParseTuple(args, "dddddd", &g, &u, &h_s, &N, &l, &eps_s))
        return NULL;

    double cap = capacitance_bare(g, u, h_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* Pycapacitance_bare_k(PyObject*, PyObject* args) {
    double k_a, k_s, l, eps_s, N;

    if (!PyArg_ParseTuple(args, "ddddd", &k_a, &k_s, &N, &l, &eps_s))
        return NULL;

    double cap = capacitance_bare_k(k_a, k_s, N, l, eps_s);
    return PyFloat_FromDouble(cap);
}


PyObject* Pycapacitance_geometric(PyObject*, PyObject* args) {
    double g, h_f, u, l, N;

    if (!PyArg_ParseTuple(args, "ddddd", &g, &u, &h_f, &N, &l))
        return NULL;

    double cap = capacitance_geometric(g, u, h_f, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* Pycapacitance_geometric_k(PyObject*, PyObject* args) {
    double k_f, l, N;

    if (!PyArg_ParseTuple(args, "ddd", &k_f, &N, &l))
        return NULL;

    double cap = capacitance_geometric_k(k_f, N, l);
    return PyFloat_FromDouble(cap);
}


PyObject* Pycapacitance_total(PyObject*, PyObject* args) {
    double g, h_f, eps_f, u, h_s, l, eps_s, N;

    if (!PyArg_ParseTuple(args, "dddddddd", &g, &u, &h_s, &h_f, &N, &l, &eps_s, &eps_f))
        return NULL;

    double cap = capacitance_total(g, u, h_s, h_f, N, l, eps_s, eps_f);
    return PyFloat_FromDouble(cap);
}


PyObject* Pycapacitance_total_k(PyObject*, PyObject* args) {
    double k_a, k_s, k_f, eps_f, l, eps_s, N;

    if (!PyArg_ParseTuple(args, "ddddddd", &k_a, &k_s, &k_f, &N, &l, &eps_s, &eps_f))
        return NULL;

    double cap = capacitance_total_k(k_a, k_s, k_f, N, l, eps_s, eps_f);
    return PyFloat_FromDouble(cap);
}


PyObject* Pydielectric_constant_relative(PyObject*, PyObject* args) {
    double C_t, C_0, g, h_f, u, l, N;

    if (!PyArg_ParseTuple(args, "ddddddd", &C_t, &C_0, &g, &u, &h_f, &N, &l))
        return NULL;

    double eps = dielectric_constant_relative(C_t, C_0, g, u, h_f, N, l);
    return PyFloat_FromDouble(eps);
}


PyObject* Pydielectric_constant_relative_k(PyObject*, PyObject* args) {
    double C_t, C_0, k_f, l, N;

    if (!PyArg_ParseTuple(args, "ddddd", &C_t, &C_0, &k_f, &N, &l))
        return NULL;

    double eps = dielectric_constant_relative_k(C_t, C_0, k_f, N, l);
    return PyFloat_FromDouble(eps);
}


static PyMethodDef idcappy_methods[] = {
    { "k_thin", (PyCFunction)Pyk_thin, METH_VARARGS, R"pbdoc(
        k_thin(g: float, u: float, h: float) -> float
        
        Calculate the modulus for an elliptic integral for small to moderate material thickness.

        Can use any units as long as u, g, & h are the same units.

        :param g: gap spacing between fingers.

        :param u: unit cell of interdigital capacitor (finger width + gap spacing).

        :param h: material thickness.

        :return: elliptic integral modulus k.
    )pbdoc" },

    { "k_thick", (PyCFunction)Pyk_thick, METH_VARARGS, R"pbdoc(
        k_thick(g: float, u: float) -> float
        
        Calculate the modulus for an elliptic integral for large material thickness (like air).

        Can use any units as long as u & g are the same units.

        :param g: gap spacing between fingers.

        :param u: unit cell of interdigital capacitor (finger width + gap spacing).

        :return: elliptic integral modulus k.
    )pbdoc" },

    { "ellint_ratio", (PyCFunction)Pyellint_ratio, METH_O,  R"pbdoc(
        ellint_ratio(k: float) -> float

        Calculate the ratio K(k)/K'(k) (starts to fail for k < 1e-7).

        :param k: the modulus.

        :return: the ratio of a complete elliptic integral of the first kind to its complement.
    )pbdoc" },

    { "ellint_inv_ratio", (PyCFunction)Pyellint_inv_ratio, METH_O,  R"pbdoc(
        ellint_inv_ratio(k: float) -> float

        Calculate the ratio K'(k)/K(k).

        :param k: the modulus.

        :return: the ratio of the complement to the complete elliptic integral of the first kind to itself.
    )pbdoc" },

    { "ellint_ratio_approx", (PyCFunction)Pyellint_ratio_approx, METH_O,  R"pbdoc(
        ellint_ratio_approx(k: float) -> float

        Calculate the ratio K(k)/K'(k) for very small k (works well for k < 0.1; essential for k < 1e-7).

        :param k: the modulus.

        :return: the ratio of a complete elliptic integral of the first kind to its complement for very small k.
    )pbdoc" },

    { "capacitance_bare", (PyCFunction)Pycapacitance_bare, METH_VARARGS, R"pbdoc(
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

    { "capacitance_bare_k", (PyCFunction)Pycapacitance_bare_k, METH_VARARGS, R"pbdoc(
        capacitance_bare(k_a: float, k_s: float, N: float, l: float, eps_s: float) -> float

        Calculate the bare capacitance of an interdigital capacitor in pF based on moduli.

        :param k_a: modulus for air contribution.

        :param k_s: modulus for substrate contribution.

        :param N: total number of fingers.

        :param l: finger length in mm.

        :param eps_s: relative dielectric constant of the substrate.

        :return: bare capacitance in pF.
    )pbdoc" },

    { "capacitance_geometric", (PyCFunction)Pycapacitance_geometric, METH_VARARGS, R"pbdoc(
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

    { "capacitance_geometric_k", (PyCFunction)Pycapacitance_geometric_k, METH_VARARGS, R"pbdoc(
        capacitance_geometric(k_f: float, N: float, l: float) -> float

        Calculate the geometric capacitance for a film in pF based on the modulus.

        :param k_f: modulus for film contribution.

        :param N: total number of fingers.

        :param l: finger length in mm.

        :return: geometric capacitance in pF.
    )pbdoc" },

    { "capacitance_total", (PyCFunction)Pycapacitance_total, METH_VARARGS, R"pbdoc(
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

    { "capacitance_total_k", (PyCFunction)Pycapacitance_total_k, METH_VARARGS, R"pbdoc(
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

    { "dielectric_constant_relative", (PyCFunction)Pydielectric_constant_relative, METH_VARARGS, R"pbdoc(
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

    { "dielectric_constant_relative_k", (PyCFunction)Pydielectric_constant_relative_k, METH_VARARGS, R"pbdoc(
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