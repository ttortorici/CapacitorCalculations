#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <Windows.h>
#include <functions.h>


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


static PyMethodDef dielepy_methods[] = {
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


static PyModuleDef dielepy_module = {
    PyModuleDef_HEAD_INIT,
    "dielepy",                                  // Module name to use with Python import statements
    "Module Description",                       // Module description
    0,
    dielepy_methods                             // Structure that defines the methods of the module
};


PyMODINIT_FUNC PyInit_dielepy(void) {
    return PyModule_Create(&dielepy_module);
}