/**
python-console
Copyright (C) 2014  Alex Tsui

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "pyinterpreter.h"
#include <Python.h>
#include <iostream>
#include <map>
#include <memory>
#include "pyredirector.h"
#include "pybindings.h"

static PyThreadState *MainThreadState = NULL;

static PyThreadState *m_threadState = NULL;
static PyObject *glb = NULL;
static PyObject *loc = NULL;

static std::list<std::string> m_suggestions;

template <typename... Args> std::string string_format(const std::string &format, Args... args)
{
    size_t size = std::snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

void pyinterpreter_execute(const std::string &command, color_output_type color)
{
    PyEval_AcquireThread(m_threadState);
    color(0);

    PyObject *py_result;
    PyObject *dum;
    py_result = Py_CompileString(command.c_str(), "<stdin>", Py_single_input);
    if (py_result == 0) {
        if (PyErr_Occurred()) {
            color(1);
            PyErr_Print();
        }

        PyEval_ReleaseThread(m_threadState);
        return;
    }
    dum = PyEval_EvalCode(py_result, glb, loc);
    Py_XDECREF(dum);
    Py_XDECREF(py_result);
    if (PyErr_Occurred()) {
        color(1);
        PyErr_Print();
    }
    PyEval_ReleaseThread(m_threadState);
}

void pyinterpreter_execute_complete(const std::string &command, color_output_type color)
{
    PyEval_AcquireThread(m_threadState);
    color(0);
    int st = PyRun_SimpleString(command.c_str());
    if (st) {
        color(1);
        PyErr_Print();
    }
    PyEval_ReleaseThread(m_threadState);
}

void pyinterpreter_run_file(const std::string &python_file, color_output_type color)
{
    PyEval_AcquireThread(m_threadState);
    color(0);
    try {
        FILE *fp = fopen(python_file.c_str(), "r");
        if (fp == NULL) {
            fprintf(stderr, "Fatal error: file not found %s\n", python_file.c_str());
            exit(1);
        }
        PyRun_SimpleFile(fp, python_file.c_str());
        fclose(fp);
    } catch (boost::python::error_already_set const &) {
        // Parse and output the exception
        //std::string perror_str = parse_python_exception();
        //std::cout << "Error in Python: " << perror_str << std::endl;
        std::cout << "Error in Python: ";
    }

    PyEval_ReleaseThread(m_threadState);
}

const std::list<std::string> &pyinterpreter_suggest(const std::string &hint)
{
    PyEval_AcquireThread(m_threadState);
    m_suggestions.clear();
    int i = 0;
    redirector_take_output(m_threadState); // to clear any previous output
    std::string command = string_format("sys.completer.complete('%s', %d)\n", hint.c_str(), i);
    std::string res;
    do {
        PyObject *py_result;
        PyObject *dum;
        py_result = Py_CompileString(command.c_str(), "<stdin>", Py_single_input);
        dum = PyEval_EvalCode(py_result, glb, loc);
        Py_XDECREF(dum);
        Py_XDECREF(py_result);
        res = redirector_take_output(m_threadState);

        ++i;
        command = string_format("sys.completer.complete('%s', %d)\n", hint.c_str(), i);
        if (res.size()) {
            // throw away the newline
            res = res.substr(1, res.size() - 3);
            m_suggestions.push_back(res);
        }
    } while (res.size());

    PyEval_ReleaseThread(m_threadState);
    return m_suggestions;
}

void pyinterpreter_preinit()
{
    m_suggestions.clear();
    inittab_redirector();
}

void pyinterpreter_initialize()
{
    PyEval_InitThreads();
    MainThreadState = PyEval_SaveThread();
    PyEval_AcquireThread(MainThreadState);
    m_threadState = Py_NewInterpreter();

    PyObject *module = PyImport_ImportModule("__main__");
    loc = glb = PyModule_GetDict(module);

    PyRun_SimpleString("import sys\n"
                       "import redirector\n"
                       "sys.path.insert(0, \".\")\n" // add current path
                       "sys.stdout = redirector.redirector()\n"
                       "sys.stderr = sys.stdout\n"
                       "import rlcompleter\n"
                       "sys.completer = rlcompleter.Completer()\n");

    PyEval_ReleaseThread(m_threadState);
}

void pyinterpreter_finalize()
{
    m_suggestions.clear();

    PyEval_AcquireThread(m_threadState);
    Py_EndInterpreter(m_threadState);
    PyEval_ReleaseLock();

    PyEval_RestoreThread(MainThreadState);
}

void pyinterpreter_aquire()
{
    PyEval_AcquireThread(m_threadState);
}

void pyinterpreter_release()
{
    PyEval_ReleaseThread(m_threadState);
}
