//
// Created by jiaruiyan on 1/19/21.
//

#include "AppReader.h"
#include "Triplet.h"
#include <unordered_set>
#include <oneapi/tbb.h>


class TriSet{
    const Eigen::MatrixXi* m_Tet;
public:
    std::unordered_set<Triplet, TripletHashFunc> my_tri_set;
    void operator()(const tbb::blocked_range<size_t>& r){
        const Eigen::MatrixXi& Tet_ref = *m_Tet;
        size_t end = r.end();
        for(size_t i = r.begin(); i != end; ++i){
            // i is the element idx.
            const Eigen::RowVector4i elemVInd = Tet_ref.row(i);
            my_tri_set.insert(Triplet(elemVInd[0], elemVInd[2], elemVInd[1]));
            my_tri_set.insert(Triplet(elemVInd[0], elemVInd[3], elemVInd[2]));
            my_tri_set.insert(Triplet(elemVInd[0], elemVInd[1], elemVInd[3]));
            my_tri_set.insert(Triplet(elemVInd[1], elemVInd[2], elemVInd[3]));
        }
    }

    TriSet(TriSet& x, tbb::split) : m_Tet(x.m_Tet), my_tri_set(std::unordered_set<Triplet, TripletHashFunc>()){}

    TriSet(const Eigen::MatrixXi* Tet)
        : m_Tet(Tet), my_tri_set(std::unordered_set<Triplet, TripletHashFunc>())
    {}

    void join(const TriSet& ts){
        std::unordered_set<Triplet, TripletHashFunc> tmp(my_tri_set);
        std::set_union(tmp.begin(), tmp.end(),
                       ts.my_tri_set.begin(), ts.my_tri_set.end(),
                       std::inserter(my_tri_set, my_tri_set.begin()));
    }

};

void AppReader::read_test_case(int id,
                    Eigen::MatrixXd &X,
                    Eigen::MatrixXi &Tet,
                    Eigen::MatrixXi &BTri){
    // Unused API variables
    Eigen::MatrixXi Tri;
    Eigen::VectorXi TriTag;
    Eigen::VectorXi TetTag;

    // Call API
    // std::string meshStr = "box3D_v518_t2112.msh";
    std::string meshStr = "armadillolow_17698v_72161t.msh";
    const std::string readfile = this->m_default_mesh_path + meshStr;
    // X and Tet are only two useful parameters here.
    if (!igl::readMSH(readfile, X, Tri, Tet, TriTag, TetTag)){
        throw std::invalid_argument("Cannot read .msh file!");
    }

    // Check read mesh
    assert(X.rows() > 0 && X.cols() == 3);
    assert(Tet.rows() > 0 && Tet.cols() == 4);

    // Find boundary triangles -- It should have a wiser way to do this.
    // Can blocked_range<iteraotr> work ?
    // Put all triangles used by tet into a set.
    TriSet ts(&Tet);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, TetTag.rows()), ts);
    std::unordered_set<Triplet, TripletHashFunc>& tri_set_ref = ts.my_tri_set;
    // Fill the BTri.
    std::map<int, int> BTriMap;
    int mapKey = 0;
    for(auto itr = tri_set_ref.begin(); itr != tri_set_ref.end(); itr++){
        const int* triVInd = itr->key;
        // find dual triangle with reversed indices:
        auto finder = tri_set_ref.find(Triplet(triVInd[2], triVInd[1], triVInd[0]));
        if(finder == tri_set_ref.end()){
            finder = tri_set_ref.find(Triplet(triVInd[1], triVInd[0], triVInd[2]));
            if (finder == tri_set_ref.end()){
                finder = tri_set_ref.find(Triplet(triVInd[0], triVInd[2], triVInd[1]));
                if(finder == tri_set_ref.end()){
                    int oldSize = BTri.rows();
                    BTri.conservativeResize(oldSize + 1, 3);
                    BTri(oldSize, 0) = triVInd[0];
                    BTri(oldSize, 1) = triVInd[1];
                    BTri(oldSize, 2) = triVInd[2];
                }
            }
        }
    }

}

