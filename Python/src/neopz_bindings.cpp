#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/iostream.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// cpp STL
#include <map>                    // for map
#include <set>                    // for set

namespace py = pybind11;
using namespace py::literals;

// NeoPZ container classes
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pztrnsform.h"

// NeoPZ topology classes
#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzpyramid.h"
#include "tpzprism.h"
#include "tpzcube.h"

// Mesh
#include "pzgmesh.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeoelbc.h"
#include "pzgeoelside.h"
#include "tpzgeoelrefpattern.h"

#include "pzcmesh.h"

// Matrix
#include "pzmatrix.h"

//StrMatrix
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZStructMatrixBase.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"

// Material
#include "TPZMaterial.h"
#include "pzpoisson3d.h"
#include "TPZMatElasticity2D.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzpostprocanalysis.h"

#include "pzsolve.h"
#include "tpzautopointer.h"

// Elastoplasticity 
#include "TPZBndCondWithMem.h" 
#include "TPZBndCondWithMem_impl.h"
#include "TPZMatElastoPlastic2D.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"
#include "TPZElastoPlasticMem.h"
#include "TPZPlasticCriterion.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZBndCondWithMem.h"

// For SBFEM simulations
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "TPZBuildSBFem.h"
#include "tpzquadraticquad.h"


using namespace std;

PYBIND11_MAKE_OPAQUE(std::map<int,int>)

string to_string(const string &value) {return value;}

namespace {

    template<typename T>
    static void declareTPZVec(py::module &mod, std::string const &suffix) {
        using Class = TPZVec<T>;
        using PyClass = py::class_<Class, std::shared_ptr<Class>>;

        PyClass cls(mod, ("TPZVec" + suffix).c_str());

        cls.def(py::init());
        cls.def(py::init<int>());
        cls.def(py::init<int64_t, T>());
        cls.def("Resize", [](Class &vec, const int64_t &newsize) { return vec.Resize(newsize); });
        cls.def("Size", [](const Class &vec) { return vec.size(); });
        cls.def("__repr__",
                [](const Class &vec) {
                    std::string r("TPZVec [");
                    r += to_string(vec[0]);
                    for (int i = 1; i < vec.NElements(); i++) {
                        r += ", ";
                        r += to_string(vec[i]);
                    }
                    r += "]";
                    return r;
                }
        );

        cls.def("__getitem__",
                [](const Class &vec, int64_t position) {
                    if (position >= vec.size() || position < 0) throw py::index_error();
                    return vec[position];
                },
                py::is_operator()
        );
        cls.def("__setitem__",
                [](Class &vec, int64_t position, T value) {
                    if (position >= vec.size() || position < 0) throw py::index_error();
                    vec[position] = value;
                },
                py::is_operator()
        );

    } //declareTPZVec


    template <typename T>
    static void declareTPZManVector(py::module & mod, std::string const & suffix) {
        using Class = TPZManVector<T>;
        using PyClass = py::class_<Class, std::shared_ptr<Class>, TPZVec<T> >;

        PyClass cls(mod, ("TPZManVector" + suffix).c_str());

        cls.def(py::init());
        cls.def(py::init<int>());
        cls.def(py::init<int64_t, T>());
        cls.def("Resize", [](Class & vec, const int64_t& newsize) { return vec.Resize(newsize); });
        cls.def("Size", [](const Class & vec) { return vec.size(); });
        cls.def("__repr__",
                [](const Class & vec) {
                    std::string r("TPZManVector [");
                    r += to_string(vec[0]);
                    for (int i = 1; i < vec.NElements(); i++) {
                        r += ", ";
                        r += to_string(vec[i]);
                    }
                    r += "]";
                    return r;
                }
        );

        cls.def("__getitem__",
                [](const Class & vec, int64_t position) {
                    if (position >= vec.size() || position < 0) throw py::index_error();
                    return vec[position];
                },
                py::is_operator()
        );
        cls.def("__setitem__",
                [](Class & vec, int64_t position, T value) {
                    if (position >= vec.size() || position < 0) throw py::index_error();
                    vec[position] = value;
                },
                py::is_operator()
        );

    } //declareTPZManVector

    template <typename T>
    static void declareTPZStack(py::module & mod, std::string const & suffix) {
        using Class = TPZStack<T>;
        using PyClass = py::class_<Class, std::shared_ptr<Class>>;

        PyClass cls(mod, ("TPZStack" + suffix).c_str());

        cls.def(py::init());
        cls.def(py::init<int, T>());
        cls.def("Push", [](Class & stack, T object) {
            stack.Push(object);
        });
        cls.def("Pop", [](Class & stack) {
            return stack.Pop();
        });
        cls.def("Peek", [](Class & stack) {
            return stack.Peek();
        });
        cls.def("__repr__",
            [](const Class & stack) {
                std::string r("TPZStack [");
                for (int i = 0; i < stack.NElements(); i++) {
                    r += to_string(stack[i]);
                    if (i != stack.NElements() - 1) {
                        r += ", ";
                    }
                }
                r += "]";
                return r;
            }
        );
    } //declareTPZStack

}// dummyNamespace

PYBIND11_MODULE(NEOPZ, m) {
    m.doc() = R"pbdoc(
        -------------------------
        Python bindings for NeoPZ
        -------------------------
    )pbdoc";


//CONTAINERS
    declareTPZVec<int>(m, "_int");
    declareTPZVec<int64_t>(m, "_int64_t");
    declareTPZVec<double>(m, "_double");
    declareTPZVec<std::string>(m, "_string");

    declareTPZManVector<int>(m, "_int");
    declareTPZManVector<int64_t>(m, "_int64_t");
    declareTPZManVector<double>(m, "_double");
    declareTPZManVector<std::string>(m, "_string");

    declareTPZStack<int>(m, "_int");
    declareTPZStack<int64_t>(m, "_int64_t");
    declareTPZStack<double>(m, "_double");
    declareTPZStack<std::string>(m, "_string");

    // TPZFMatrix<double> bindings
    py::class_<TPZFMatrix<double>>(m, "TPZFMatrix")
        .def(py::init())
        .def(py::init<int64_t, int64_t>())
        .def(py::init<int64_t, int64_t, double>())
        .def("GetVal", [](TPZFMatrix<double>& matrix, const int64_t& row, const int64_t& col) {
            if (row >= matrix.Rows() || row < 0) throw py::index_error();
            if (col >= matrix.Cols() || col < 0) throw py::index_error();
            return matrix.GetVal(row, col);
        })
        .def("SetItem",
             [](TPZFMatrix<double>& mat, int64_t rows, int64_t cols, double value) {
                 if (rows >= mat.Rows() || rows < 0) throw py::index_error();
                 if (cols >= mat.Cols() || cols < 0) throw py::index_error();
                 mat(rows,cols) = value;
             }//,             py::is_operator()
        )
        .def("__repr__",
             [](TPZFMatrix<double>& matrix) {
                 std::string r("TPZFMatrix ");
                 r += "'(";
                 r += std::to_string(matrix.Rows());
                 r += " x ";
                 r += std::to_string(matrix.Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < matrix.Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < matrix.Cols(); col++) {
                         r += std::to_string(matrix.GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]";
                 return r;
             }
        )   

        .def("Zero", & TPZFMatrix<double>::Zero)
        .def("Cols", & TPZFMatrix<double>::Cols)
        .def("Rows", & TPZFMatrix<double>::Rows)
        .def("Resize", & TPZFMatrix<double>::Resize)
        .def("SetSize", & TPZFMatrix<double>::SetSize)

        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self *= float())

            //Global Functions
        .def("Norm", [](TPZFMatrix<double> &matrix){ return Norm(matrix);})
        // .def("Norm", &TPZFMatrix<double>::Norm)

         ;


    py::class_<TPZMatrix<double> >(m, "TPZMatrix")
        // .def(py::init<>())
        .def("SetIsDecomposed", & TPZMatrix<double>::SetIsDecomposed)
        
        .def("__repr__", [](TPZMatrix<double> &matrix){
                std::ofstream printstream;
                matrix.Print("Matrix = ", printstream);
                return printstream;
            }
        )
    ;

    // TPZAdmChunkVectorNodes bindings
    py::class_<TPZAdmChunkVector<TPZGeoNode>>(m, "TPZAdmChunkVectorNodes")
        .def(py::init<int>())
        .def("NElements", &TPZAdmChunkVector<TPZGeoNode>::NElements)
        .def("Resize", &TPZAdmChunkVector<TPZGeoNode>::Resize)
        .def("__getitem__",
             [](const TPZAdmChunkVector<TPZGeoNode>& vec, int64_t position) {
                 if (position >= vec.NElements() || position < 0) throw py::index_error();
                 return vec[position];
             },
             py::is_operator()
        )
    ;

    // TPZTransform bindings
    py::class_<TPZTransform<double>>(m, "TPZTransform")
        .def(py::init())
        .def(py::init<int>())
        .def(py::init<int, int>())
        .def("Mult", py::overload_cast<>(&TPZTransform<double>::Mult))
        .def("Sum", py::overload_cast<>(&TPZTransform<double>::Sum))
        .def("Mult", py::overload_cast<>(&TPZTransform<double>::Mult, py::const_))
        .def("Sum", py::overload_cast<>(&TPZTransform<double>::Sum, py::const_))
        .def("SetMatrix", &TPZTransform<double>::SetMatrix)
        .def("Multiply", &TPZTransform<double>::Multiply)
        .def("Apply", &TPZTransform<double>::Apply)
        .def("__repr__",
             [](const TPZTransform<double>& trans) {
                //  std::stringstream repr;
                //  trans.Print(repr);
                //  return repr.str();
                 std::string r("TPZTransform\n");
                 r += "  Mult ";
                 r += "'(";
                 r += std::to_string(trans.Mult().Rows());
                 r += " x ";
                 r += std::to_string(trans.Mult().Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < trans.Mult().Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < trans.Mult().Cols(); col++) {
                         r += std::to_string(trans.Mult().GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]\n";
                 r += "  Sum";
                 r += "'(";
                 r += std::to_string(trans.Sum().Rows());
                 r += " x ";
                 r += std::to_string(trans.Sum().Cols());
                 r += ")' = [\n";
                 for (int64_t row = 0; row < trans.Sum().Rows(); row++) {
                     r += "\t";
                     for (int64_t col = 0; col < trans.Sum().Cols(); col++) {
                         r += std::to_string(trans.Sum().GetVal(row, col));
                         r += "  ";
                     }
                     r += "\n";
                 }
                 r += "]";
                 return r;
             }
        )
    ;

    // TPZPoint bindings
    py::class_<pztopology::TPZPoint>(m, "TPZPoint")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPoint::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPoint::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPoint::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPoint::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPoint::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPoint::NumSides))
        .def_static("CenterPoint", &pztopology::TPZPoint::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPoint& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPoint::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPoint::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPoint::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZPoint::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPoint::CreateSideIntegrationRule)
    ;

    // TPZLine bindings
    py::class_<pztopology::TPZLine>(m, "TPZLine")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZLine::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZLine::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZLine::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZLine::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZLine::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZLine::NumSides))
        .def_static("CenterPoint", &pztopology::TPZLine::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZLine& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZLine::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZLine::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZLine::TransformElementToSide)
        .def_static("IsInParametricDomain", &pztopology::TPZLine::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZLine::CreateSideIntegrationRule)
    ;

    // TPZTriangle bindings
    py::class_<pztopology::TPZTriangle>(m, "TPZTriangle")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZTriangle::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTriangle::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTriangle::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTriangle::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTriangle::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZTriangle::NumSides))
        .def_static("CenterPoint", &pztopology::TPZTriangle::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTriangle& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTriangle::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTriangle::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTriangle::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZTriangle::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTriangle::CreateSideIntegrationRule)
    ;


    // TPZQuadrilateral bindings
    py::class_<pztopology::TPZQuadrilateral>(m, "TPZQuadrilateral")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZQuadrilateral::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZQuadrilateral::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZQuadrilateral::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZQuadrilateral::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZQuadrilateral::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZQuadrilateral::NumSides))
        .def_static("CenterPoint", &pztopology::TPZQuadrilateral::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZQuadrilateral& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZQuadrilateral::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZQuadrilateral::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZQuadrilateral::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZQuadrilateral::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZQuadrilateral::CreateSideIntegrationRule)
    ;

    // TPZTetrahedron bindings
    py::class_<pztopology::TPZTetrahedron>(m, "TPZTetrahedron")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZTetrahedron::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZTetrahedron::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZTetrahedron::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZTetrahedron::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZTetrahedron::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZTetrahedron::NumSides))
        .def_static("CenterPoint", &pztopology::TPZTetrahedron::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZTetrahedron& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZTetrahedron::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZTetrahedron::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZTetrahedron::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZTetrahedron::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZTetrahedron::CreateSideIntegrationRule)
    ;

    // TPZPyramid bindings
    py::class_<pztopology::TPZPyramid>(m, "TPZPyramid")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPyramid::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPyramid::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPyramid::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPyramid::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPyramid::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPyramid::NumSides))
        .def_static("CenterPoint", &pztopology::TPZPyramid::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPyramid& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPyramid::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPyramid::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPyramid::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZPyramid::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPyramid::CreateSideIntegrationRule)
    ;

    // TPZPrism bindings
    py::class_<pztopology::TPZPrism>(m, "TPZPrism")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZPrism::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZPrism::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZPrism::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZPrism::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZPrism::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZPrism::NumSides))
        .def_static("CenterPoint", &pztopology::TPZPrism::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZPrism& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZPrism::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZPrism::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZPrism::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZPrism::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZPrism::CreateSideIntegrationRule)
    ;

    // TPZCube bindings
    py::class_<pztopology::TPZCube>(m, "TPZCube")
        .def(py::init())
        .def_static("SideDimension", &pztopology::TPZCube::SideDimension)
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("LowerDimensionSides", py::overload_cast<int, TPZStack<int> &, int>(&pztopology::TPZCube::LowerDimensionSides))
        .def_static("HigherDimensionSides", &pztopology::TPZCube::HigherDimensionSides)
        .def_static("NSideNodes", &pztopology::TPZCube::NSideNodes)
        .def_static("SideNodeLocId", &pztopology::TPZCube::SideNodeLocId)
        .def_static("NumSides", py::overload_cast<>(&pztopology::TPZCube::NumSides))
        .def_static("CenterPoint", &pztopology::TPZCube::CenterPoint)
        .def_static("RefElVolume", [](pztopology::TPZCube& topology) { return topology.RefElVolume(); })
        .def_static("SideToSideTransform", &pztopology::TPZCube::SideToSideTransform)
        .def_static("TransformSideToElement", &pztopology::TPZCube::TransformSideToElement)
        .def_static("TransformElementToSide", &pztopology::TPZCube::TransformElementToSide)
        .def_static("IsInParametricDomain", (bool (*) (const TPZVec<REAL>&, REAL)) &pztopology::TPZCube::IsInParametricDomain)
        .def_static("CreateSideIntegrationRule", &pztopology::TPZCube::CreateSideIntegrationRule)
    ;

    // TPZGeoMesh bindings
    py::class_<TPZGeoMesh>(m, "TPZGeoMesh")
        .def(py::init())
        .def("Print", [](TPZGeoMesh &GeoMesh){ return GeoMesh.Print();})
        .def("BuildConnectivity", &TPZGeoMesh::BuildConnectivity)
        .def("NElements", &TPZGeoMesh::NElements)
        .def("Dimension", &TPZGeoMesh::Dimension)
        .def("NNodes", &TPZGeoMesh::NNodes)
        .def("NodeVec", py::overload_cast<>(&TPZGeoMesh::NodeVec))
        .def("NElements", &TPZGeoMesh::NElements)
        .def("Element", &TPZGeoMesh::Element)
        .def("ReadUNSWSBGeoFile",[](TPZGeoMesh & self, const std::string &filename, int ESkeleton, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices)
    {

        int maxvol = -1;

        std::ifstream file(filename);

        map<set<int64_t> , int64_t> midnode;
        string buf;
        getline(file,buf);
        if(!file) DebugStop();
        int64_t nnodes, nvolumes;
        file >> nnodes >> nvolumes;
        elpartition.Resize(nvolumes*6, -1);
        TPZGeoMesh *gmesh = new TPZGeoMesh;
        gmesh->SetDimension(3);
        gmesh->NodeVec().Resize(nnodes);
        for (int64_t in=0; in<nnodes; in++) {
            TPZManVector<REAL,3> xco(3);
            for (int i=0; i<3; i++) {
                file >> xco[i];
            }
            gmesh->NodeVec()[in].Initialize(xco, *gmesh);
        }
#ifdef PZDEBUG
        std::set<int64_t> badvolumes;
#endif
        int64_t nothing;
        file >> nothing;
        for (int64_t iv=0; iv<nvolumes; iv++) {
#ifdef PZDEBUG
            map<set<int64_t>,int64_t> nodepairs;
#endif
            int nfaces;
            file >> nfaces;
            for (int face = 0; face < nfaces; face++) {
                int elnnodes;
                file >> elnnodes;

                TPZManVector<int64_t,10> nodes(elnnodes);
                for (int i=0; i<elnnodes; i++) {
                    file >> nodes[i];
                    nodes[i]--;
#ifdef PZDEBUG
                    if (i>0) {
                        set<int64_t> edge;
                        edge.insert(nodes[i-1]);
                        edge.insert(nodes[i]);
                        nodepairs[edge]++;
                    }
                    if (i==elnnodes-1) {
                        set<int64_t> edge;
                        edge.insert(nodes[0]);
                        edge.insert(nodes[i]);
                        nodepairs[edge]++;
                    }
#endif
                }

                // tototototo
                if (maxvol != -1 && iv >= maxvol) {
                    continue;
                }
                if (elnnodes == 1)
                {
                    int64_t index;
                    MElementType eltype = EPoint;
                    gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                    elpartition[index] = iv;

                }
                else if (elnnodes == 2)
                {
                    int64_t index;
                    MElementType eltype = EOned;
                    gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                    elpartition[index] = iv;

                }
                else if (elnnodes == 3 || elnnodes == 4)
                {
                    int64_t index;
                    MElementType eltype = ETriangle;
                    if (elnnodes == 4) {
                        eltype = EQuadrilateral;
                    }
                    gmesh->CreateGeoElement(eltype, nodes, ESkeleton, index);
                    elpartition[index] = iv;
                }
                else if(elnnodes == 8)
                {
                    int64_t index;
                    new TPZGeoElRefPattern<pzgeom::TPZQuadraticQuad> (nodes, ESkeleton, *gmesh,  index);
                    elpartition[index] = iv;

                }
                else if(elnnodes > 4)
                {
                    set<int64_t>  elnodes;
                    TPZManVector<REAL,3> midxco(3,0.);
                    for (int i=0; i<elnnodes; i++) {
                        elnodes.insert(nodes[i]);
                        TPZManVector<REAL,3> x(3);
                        gmesh->NodeVec()[nodes[i]].GetCoordinates(x);
                        //                    std::cout << "x " << x << endl;
                        for(int j=0; j<3; j++) midxco[j] += x[j]/elnnodes;
                    }
                    int64_t midindex = -1;
                    if (midnode.find(elnodes) == midnode.end()) {
                        midindex = gmesh->NodeVec().AllocateNewElement();
                        gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);
                        midnode[elnodes] = midindex;
                    }
                    else
                    {
                        midindex = midnode[elnodes];
                    }
                    for (int triangle = 0; triangle <elnnodes; triangle++) {
                        TPZManVector<int64_t,3> nodeindices(3);
                        for (int in=0; in<2; in++) {
                            nodeindices[in] = nodes[(triangle+in)%elnnodes];
                        }
                        nodeindices[2] = midindex;
                        int64_t index;
                        gmesh->CreateGeoElement(ETriangle, nodeindices, ESkeleton, index);
                        elpartition[index] = iv;
                    }
                }
                else
                {
                    DebugStop();
                }
            }
#ifdef PZDEBUG
            bool suspicious = false;
            for (auto it = nodepairs.begin(); it != nodepairs.end(); it++) {
                if(it->second != 2) suspicious = true;
            }
            if (suspicious == true) {
                std::cout << "volume " << iv << " has no closure\n";
                badvolumes.insert(iv);
            }
#endif
            if (elpartition.size() < gmesh->NElements()+100) {
                elpartition.Resize(elpartition.size()*2, -1);
            }
        }
        // totototototo
        if (maxvol != -1) {
            nvolumes = maxvol;
        }
        int64_t nmidnodes = midnode.size();
        gmesh->NodeVec().Resize(nvolumes+nmidnodes+nnodes);
        scalingcenterindices.Resize(nvolumes, -1);
        for (int64_t in=0; in<nvolumes; in++) {
            TPZManVector<REAL,3> xco(3);
            for (int i=0; i<3; i++) {
                file >> xco[i];
            }
            gmesh->NodeVec()[nnodes+nmidnodes+in].Initialize(xco, *gmesh);
            scalingcenterindices[in] = nnodes+nmidnodes+in;
        }
        {
            ofstream mirror("gmesh.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, mirror);
        }
#ifdef PZDEBUG
        if (badvolumes.size()) {
            int64_t nel = gmesh->NElements();
            TPZManVector<REAL> elval(nel,0);
            for (int64_t el=0; el<nel; el++) {
                if (badvolumes.find(elpartition[el]) != badvolumes.end()) {
                    elval[el] = 10.;
                }
            }
            {
                ofstream badel("gmesh_bad.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, badel, elval);
            }
        }
#endif
        elpartition.Resize(gmesh->NElements(), -1);
        std::cout << "Building element connectivity\n";
        gmesh->BuildConnectivity();
        return gmesh;
    })
    ;

    // TPZGeoNode bindings
    py::class_<TPZGeoNode>(m, "TPZGeoNode")
        .def(py::init())
        .def("GetCoordinates", &TPZGeoNode::GetCoordinates)
    ;

    // TPZGeoEl
    py::class_<TPZGeoEl, std::unique_ptr<TPZGeoEl, py::nodelete>  >(m, "TPZGeoEl")
        .def("NSides", &TPZGeoEl::NSides)
        .def("NSideNodes", &TPZGeoEl::NSideNodes)
        .def("SideNodeIndex", &TPZGeoEl::SideNodeIndex)
        .def("SideDimension", &TPZGeoEl::SideDimension)
        .def("Dimension", &TPZGeoEl::Dimension)
        .def("Neighbour", &TPZGeoEl::Neighbour)
    ;

    // TPZGeoElBC
    py::class_<TPZGeoElBC >(m, "TPZGeoElBC")
        .def(py::init<TPZGeoEl*, int, int>())
    ;

    // TPZGeoElSide
    py::class_<TPZGeoElSide, std::unique_ptr<TPZGeoElSide, py::nodelete> >(m, "TPZGeoElSide")
        .def(py::init())
        .def(py::init<TPZGeoEl*, int>())
        .def("Neighbour", &TPZGeoElSide::Neighbour)
    ;

    // TPZGMshReader
    py::class_<TPZGmshReader>(m, "TPZGmshReader")
        .def(py::init())
        .def("GeometricGmshMesh3", &TPZGmshReader::GeometricGmshMesh3, "Reads geometric mesh from GMsh (3.x) .msh file.") 
        .def("GeometricGmshMesh4", &TPZGmshReader::GeometricGmshMesh4, "Reads geometric mesh from GMsh (4.x) .msh file.") 
    ;

    //TPZMaterial
    py::class_<TPZMaterial, std::unique_ptr<TPZMaterial, py::nodelete>>(m, "TPZMaterial")
        .def("CreateBC", &TPZMaterial::CreateBC)
        .def("SetId", &TPZMaterial::SetId)
        .def("Id", &TPZMaterial::Id)
    ;

    py::class_<TPZBndCond, TPZMaterial , std::unique_ptr<TPZBndCond, py::nodelete>>(m, "TPZBndCond")
        .def(py::init())
        .def(py::init<int>())
        .def(py::init< TPZMaterial * ,int ,int  , TPZFMatrix<STATE> & ,TPZFMatrix<STATE> & >())
    ;


    py::class_<TPZBndCondWithMem<TPZElastoPlasticMem>, TPZBndCond, std::unique_ptr<TPZBndCondWithMem<TPZElastoPlasticMem>, py::nodelete>>(m, "TPZBndCondWithMem")
        .def(py::init())
        .def(py::init< TPZMaterial * ,int ,int  , TPZFMatrix<STATE> & ,TPZFMatrix<STATE> & >())
    ;


    py::class_<TPZMatPoisson3d, TPZMaterial, std::unique_ptr<TPZMatPoisson3d, py::nodelete>>(m, "TPZMatPoisson3d")
        .def(py::init<int, int>())
    ;

    py::class_<TPZMatElastoPlastic< TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >, TPZMaterial, std::unique_ptr<TPZMatElastoPlastic< TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >, py::nodelete>>(m, "TPZMatElastoPlasticMC")
        .def(py::init())
    ;

    py::class_<TPZMatWithMem<TPZElastoPlasticMem>, TPZMaterial, std::unique_ptr<TPZMatWithMem<TPZElastoPlasticMem>, py::nodelete>>(m, "TPZMatWithMem")
    // .def("SetUpdateMem", & TPZMatWithMem::SetUpdateMem)
    // .def("SetUpdateMem", [](TPZMatWithMem<TPZElastoPlasticMem> & plasticupdate){ return plasticupdate.SetUpdateMem;})
    ;

    py::class_<TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >,TPZMatElastoPlastic< TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >, TPZMatWithMem<TPZElastoPlasticMem>, TPZMaterial, std::unique_ptr<TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >, py::nodelete>>(m, "TPZMatElastoPlastic2DMC")
        .def(py::init<int, int>())
        .def("SetPlasticityModel", &TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >::SetPlasticityModel)
        .def("SetDefaultMem", &TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >::SetDefaultMem)
    ;

    py::class_<TPZElasticResponse>(m, "TPZElasticResponse")
        .def(py::init())
        .def("SetEngineeringData", & TPZElasticResponse::SetEngineeringData)
    ;

    py::class_<TPZTensor<REAL>>(m, "TPZTensor")
        .def(py::init())
        .def("Zero", & TPZTensor<REAL>::Zero)
    ;

    py::class_<TPZPlasticState<STATE>>(m, "TPZPlasticState")
        .def(py::init())
    ;

    py::class_<TPZElastoPlasticMem>(m, "TPZElastoPlasticMem")
        .def(py::init())
        .def("SetElasticResponse", [](TPZElastoPlasticMem &self, TPZElasticResponse &ER){
            self.m_ER = ER;
            return;
        })
        .def("SetStress", [](TPZElastoPlasticMem &self, TPZTensor<REAL> &stress){
            self.m_sigma = stress;
            return;
        })
        .def("SetPlasticState", [](TPZElastoPlasticMem &self, TPZPlasticState<REAL> &plastic_state){
            self.m_elastoplastic_state = plastic_state;
            return;
        })
    ;

    
    py::class_<TPZPlasticCriterion, std::unique_ptr<TPZPlasticCriterion, py::nodelete>>(m, "TPZPlasticCriterion")
    ;

    py::class_<TPZYCMohrCoulombPV, TPZPlasticCriterion>(m, "TPZYCMohrCoulombPV")
        .def(py::init())
        .def("SetUp", & TPZYCMohrCoulombPV::SetUp)
    ;

    py::class_<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, std::unique_ptr<TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, py::nodelete>>(m, "TPZPlasticStepPVMC")
        .def(py::init())
        .def("SetElasticResponse", & TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>::SetElasticResponse)
        // .def("GetYC", & TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>::GetYC)
        .def("YC", [](TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> & plasticstep){ return plasticstep.fYC;})
        .def("fN", [](TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> & plasticstep){ return plasticstep.fN;})
    ;
 
    py::class_<TPZMatElasticity2D, TPZMaterial >(m, "TPZMatElasticity2D")
        .def(py::init<int>())
    ;
    
    py::class_<TPZElasticity3D, TPZMaterial, std::unique_ptr<TPZElasticity3D, py::nodelete>>(m, "TPZElasticity3D")
        .def(py::init<int>())
        .def("SetMaterialDataHook", & TPZElasticity3D::SetMaterialDataHook)
    ;
    

    py::class_<TPZCreateApproximationSpace, std::unique_ptr<TPZCreateApproximationSpace, py::nodelete>>(m, "TPZCreateApproximationSpace")
        .def(py::init())
        .def("NeedsMemory", [](TPZCreateApproximationSpace &approximateSpace){ return approximateSpace.NeedsMemory();})
        .def("CreateWithMemory", & TPZCreateApproximationSpace::CreateWithMemory) 
    ;


    py::class_<TPZCompMesh , std::unique_ptr<TPZCompMesh, py::nodelete>>(m, "TPZCompMesh")
        .def(py::init())
        .def(py::init<TPZGeoMesh *>())
        .def("AutoBuild", [](TPZCompMesh &compmesh){ return compmesh.AutoBuild();})
        .def("SetDimModel", &TPZCompMesh::SetDimModel )
        .def("InsertMaterialObject", [](TPZCompMesh &compmesh, TPZMaterial *mat){ return compmesh.InsertMaterialObject(mat);} )
        .def("SetAllCreateFunctionsContinuous", &TPZCompMesh::SetAllCreateFunctionsContinuous)
        .def("SetAllCreateFunctionsContinuousWithMem", &TPZCompMesh::SetAllCreateFunctionsContinuousWithMem)  
        .def("NMaterials", &TPZCompMesh::NMaterials )
        .def("NElements", &TPZCompMesh::NElements)
        .def("Print", [](TPZCompMesh &compmesh){ return compmesh.Print();})
        .def("SetDefaultOrder",&TPZCompMesh::SetDefaultOrder)
        .def("FindMaterial", &TPZCompMesh::FindMaterial)
        .def("NEquations", &TPZCompMesh::NEquations)                
        .def("ApproxSpace", [](TPZCompMesh &compmesh){ return compmesh.ApproxSpace();})
        .def("Solution", &TPZCompMesh::Solution) 
        .def("__repr__",
             [](TPZCompMesh & comp) {
                 std::ofstream printstream;
                 comp.Print(printstream);
                 return printstream;
             }
        )
    ;
    

    
    

    py::enum_<DecomposeType>(m, "DecomposeType")
        .value("ECholesky", DecomposeType::ECholesky)
        .value("ELDLt", DecomposeType::ELDLt)
        .value("ELU", DecomposeType::ELU)
        .export_values()
    ;

    py::class_<TPZStructMatrix >(m, "TPZStructMatrix")
        // .def("__repr__",
        //      [](TPZStructMatrix & matrix) {
        //          std::ofstream printstream;
        //          matrix.Print(printstream);
        //          return printstream;
        //      }
        // )
    ;

    py::class_<TPZFStructMatrix,TPZStructMatrix >(m, "TPZFStructMatrix")

        .def(py::init<TPZCompMesh *>())
    ;


    py::class_<TPZSkylineStructMatrix, TPZStructMatrix >(m, "TPZSkylineStructMatrix")
        
        .def(py::init<TPZCompMesh *>())
    ;


    py::class_<TPZSkylineNSymStructMatrix, TPZSkylineStructMatrix >(m, "TPZSkylineNSymStructMatrix")
        .def(py::init<TPZCompMesh *>())
    ;


    py::class_<TPZFrontStructMatrix<TPZFrontSym<STATE> >, TPZStructMatrix>(m, "TPZFrontStructMatrix")
    .def(py::init<TPZCompMesh *>())
    ;

    py::class_<TPZParFrontStructMatrix<TPZFrontSym<STATE> >, TPZFrontStructMatrix<TPZFrontSym<STATE> >>(m, "TPZParFrontStructMatrix")
    .def(py::init<TPZCompMesh *>())
    ;


    py::class_<TPZStructMatrixBase >(m, "TPZStructMatrixBase")
    .def("SetNumThreads", &TPZStructMatrixBase::SetNumThreads)
    ;

    py::class_<TPZSymetricSpStructMatrix, TPZStructMatrix >(m, "TPZSymetricSpStructMatrix")
        .def(py::init<TPZCompMesh *>())
        .def("Create", [](TPZSymetricSpStructMatrix & spmatrix) {
        })
        .def("SetupMatrixData", [](TPZSymetricSpStructMatrix & spmatrix,TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex) {
            return spmatrix.SetupMatrixData(elgraph, elgraphindex);
    })   
    ;

    py::class_<TPZAnalysis >(m, "TPZAnalysis")
        .def(py::init())
        .def(py::init<TPZCompMesh *, bool>())
        .def("SetStructuralMatrix",py::overload_cast<TPZStructMatrix &>(&TPZAnalysis::SetStructuralMatrix))
        .def("SetSolver", &TPZAnalysis::SetSolver)
        // .def("Solver", &TPZAnalysis::Solver)
        // .def("Solver", [](TPZAnalysis &self) ->TPZMatrixSolver<STATE>& {
        //     return self.Solver();
        // })
        .def("PrintMatrix", [](TPZAnalysis &self){
            self.Solver().Matrix()->Print("K = ",std::cout,EMathematicaInput);
            return 0;
        })

          .def("PrintRhs", [](TPZAnalysis &self){
            self.Rhs().Print("R = ", std::cout,EMathematicaInput);
            return 0;
        })
        .def("SetStructMatrixDecomposed", [](TPZAnalysis &self, bool key = true){
            self.Solver().Matrix()->SetIsDecomposed(key);
            return key;
        })
        .def("LoadSolution", py::overload_cast<const TPZFMatrix<STATE> & >(&TPZAnalysis::LoadSolution))
        .def("Solution", [](TPZAnalysis &sol){ return sol.Solution();})
        .def("Assemble", &TPZAnalysis::Assemble)
        .def("AssembleResidual", &TPZAnalysis::AssembleResidual)
        .def("Solve", &TPZAnalysis::Solve)
        .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const std::string &  >(&TPZAnalysis::DefineGraphMesh))
        .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const TPZVec<std::string>&, const std::string &  >(&TPZAnalysis::DefineGraphMesh))
        .def("PostProcess", py::overload_cast<int, int>(&TPZAnalysis::PostProcess))
        .def("Mesh", &TPZAnalysis::Mesh)
        .def("Rhs",&TPZAnalysis::Rhs)
//        .def("SetRhs", &TPZAnalysis::SetRhs)
        .def("AcceptPseudoTimeStepSolution",[](TPZAnalysis & self){
            TPZCompMesh *cmesh = self.Mesh();
            bool update = true;
            {
                std::map<int, TPZMaterial *> &refMatVec = cmesh->MaterialVec();
                TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
                for(auto mit = refMatVec.begin(); mit != refMatVec.end(); mit++){
                    pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(mit->second);
                    if(pMatWithMem != NULL){
                        pMatWithMem->SetUpdateMem(update);
                    }
                }
            }
            self.AssembleResidual();
            update = false;
            {
                std::map<int, TPZMaterial *> &refMatVec = cmesh->MaterialVec();
                TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
                for(auto mit = refMatVec.begin(); mit != refMatVec.end(); mit++){
                    pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>(mit->second);
                    if(pMatWithMem != NULL){
                        pMatWithMem->SetUpdateMem(update);
                    }
                }
            }
        })
    ;


    py::class_<TPZAutoPointer<STATE> >(m, "TPZAutoPointer")

    ;


    py::class_<TPZSolver<STATE> >(m, "TPZSolver")
    ;

    py::class_<TPZMatrixSolver<STATE>, TPZSolver<STATE> >(m, "TPZMatrixSolver")
 
        .def("Matrix", [](TPZMatrixSolver<STATE> & self) ->TPZMatrix<STATE>&  {
            // TPZMatrix<STATE> matrix 
            // matrix = *(self.Matrix().operator->());
            return *(self.Matrix().operator->());
        })

    ;

    
    py::class_<TPZStepSolver<STATE>, TPZMatrixSolver<STATE> >(m, "TPZStepSolver")
        .def(py::init())
        .def("SetDirect", &TPZStepSolver<STATE>::SetDirect)
    ;

    py::class_<TPZPostProcAnalysis, TPZAnalysis, std::unique_ptr<TPZPostProcAnalysis, py::nodelete> >(m, "TPZPostProcAnalysis")
        .def(py::init())
        .def("SetCompMesh", &TPZPostProcAnalysis::SetCompMesh)
        .def("SetPostProcessVariables", &TPZPostProcAnalysis::SetPostProcessVariables)
        .def("SetStructuralMatrix",py::overload_cast<TPZStructMatrix &>(&TPZAnalysis::SetStructuralMatrix))
        .def("TransferSolution", &TPZPostProcAnalysis::TransferSolution)
        // .def("DefineGraphMesh", py::overload_cast<int,const TPZVec<std::string> &,const TPZVec<std::string>&, const TPZVec<std::string>&, const std::string &  >(&TPZPostProcAnalysis::DefineGraphMesh))
        // .def("PostProcess", py::overload_cast<int, int>(&TPZPostProcAnalysis::PostProcess))
    ;

//    py::class_<TPZSBFemVolume, std::unique_ptr<TPZSBFemVolume, py::nodelete> >(m, "TPZSBFemVolume")
//        .def(py::init())
//    ;

    // TPZGeoMesh bindings
    py::class_<TPZVTKGeoMesh, std::unique_ptr<TPZVTKGeoMesh, py::nodelete>>(m, "TPZVTKGeoMesh")
        .def(py::init())
        .def_static("PrintGMeshVTK",  py::overload_cast<TPZGeoMesh*, const char *, int>(&TPZVTKGeoMesh::PrintGMeshVTK))
    ;

    py::bind_map<std::map<int, int>>(m, "MapIntInt");

    // TPZBuildSBFem bindings
    py::class_<TPZBuildSBFem, std::unique_ptr<TPZBuildSBFem, py::nodelete>>(m, "TPZBuildSBFem")
        .def(py::init<TPZGeoMesh*, int, std::map<int,int> &>())
        .def("SetPartitions", &TPZBuildSBFem::SetPartitions)
        .def("BuildComputationalMeshFromSkeleton", &TPZBuildSBFem::BuildComputationalMeshFromSkeleton)
    ;


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
