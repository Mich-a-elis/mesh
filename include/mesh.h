#pragma once

#define NC "\e[0m"
#define RED "\e[0;31m"
#define BLU "\e[0;34m"
#define H_BLU "\e[0;94m"
#define GRN "\e[0;32m"
#define PURP "\e[0;35m"
#define CYN "\e[0;36m"
#define REDB "\e[041m"
#define YELLOW "\e[0;33m"

#include <iostream>
#include <vector>
#include <filesystem>
#include <memory>
#include <string>
#include <cassert>
#include <map>
#include <fstream>
#include <algorithm>
#include <complex>
#include <climits>
#include <cmath>

#include "vect.h"

namespace mesh
{
    class Vertex;
    class Edge;
    class Face;
    class HalfEdge;

    typedef std::vector<HalfEdge>::iterator HalfEdgeIter;
    typedef std::vector<HalfEdge>::const_iterator HalfEdgeCIter;

    class Vertex
    {
    public:
        Vertex(){};
        Vertex(const Vect &coords) { this->coords = coords; }
        ~Vertex(){};

        std::vector<std::shared_ptr<Face>> get_faces();

        bool on_bound = false;
        std::map<int, HalfEdgeIter> tag2he;
        std::map<int, int> tag2num;
        std::map<int, Vect> tag2oZ;
        std::map<int, Vect> tag2oX;
        Vect coords;
        std::vector<std::weak_ptr<Face>> faces;
        std::string nature = "None"; // None, sing, intersection
    };
    inline std::ostream &operator<<(std::ostream &os, const Vertex &vert);
    std::vector<std::shared_ptr<mesh::Face>> verts_faces_intersection(const std::shared_ptr<mesh::Vertex> &vert_a, const std::shared_ptr<mesh::Vertex> &vert_b);

    class Edge
    {
    public:
        Edge(const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b);
        ~Edge(){};

        std::vector<std::shared_ptr<Face>> get_faces();

        double lenght;
        bool is_wall = false;
        std::map<int, HalfEdgeIter> tag2he;
        std::map<int, int> tag2num;
        std::vector<std::shared_ptr<Vertex>> verts;
    };
    std::ostream &operator<<(std::ostream &os, const Edge &edge);

    class Face : public std::enable_shared_from_this<Face>
    {
    public:
        Face(const std::vector<std::shared_ptr<Vertex>> &verts);
        ~Face(){};

        // HalfEdgeIter halfEdge() { return face_halfEdge; }
        std::pair<bool, Vect> in_face(const Vect &pt);

        HalfEdgeIter he;
        std::map<int, int> tag2num;
        int flag = -1;
        int tag = -1;
        double area;
        Vect bary;
        std::vector<std::shared_ptr<Vertex>> verts;
        std::vector<std::shared_ptr<Edge>> edges;
        Vect normal;
        Vect oX;
    };
    std::ostream &operator<<(std::ostream &os, const Face &face);

    inline std::vector<std::shared_ptr<Vertex>> get_verts_from_faces_intersect(const std::shared_ptr<Face> &face_1, const std::shared_ptr<Face> &face_2)
    {
        std::vector<std::shared_ptr<Vertex>> res;
        for (const auto &v : face_1->verts)
        {
            for (const auto &vv : face_2->verts)
            {
                if (v == vv)
                    res.push_back(v);
            }
        }
        return res;
    }

    class HalfEdge
    {
    public:
        HalfEdge(){};
        ~HalfEdge(){};

        int num = -1;
        double ang_corner = std::numeric_limits<double>::quiet_NaN(); // Corner angle at the origin of half-edge in the face associated with the half-edge."
        std::complex<double> localVec;
        HalfEdgeIter next;                        // next halfedge in face
        HalfEdgeIter prev;                        // prev halfedge in face
        HalfEdgeIter twin;                        // twin halfedge across edge
        std::shared_ptr<Vertex> origin = nullptr; // originating vertex
        std::shared_ptr<Vertex> target = nullptr; // target vertex
        std::shared_ptr<Edge> egde = nullptr;     // incident edge
        std::shared_ptr<Face> face = nullptr;     // incident face
    };

    struct SubMesh
    {
        int tag = -1;
        std::vector<std::shared_ptr<Vertex>> verts;
        std::vector<std::shared_ptr<Edge>> edges;
        std::vector<std::shared_ptr<Face>> faces;
    };

    struct Data_Visu
    {
        std::vector<double> face_scalar = std::vector<double>();
        std::vector<double> vert_scalar = std::vector<double>();
        std::vector<Vect> vert_vector = std::vector<Vect>();
        std::vector<Vect> vert_cross = std::vector<Vect>();
        std::string output_directory = "";
    };

    class Mesh : public std::enable_shared_from_this<Mesh>
    {
    public:
        Mesh(){};
        Mesh(const std::filesystem::path &mesh_path);
        // Mesh(const Mesh &mesh);
        // const Mesh &operator=(const Mesh &mesh);
        ~Mesh(){};

        void read_mesh(const std::filesystem::path &mesh_path);

        void make_mesh();
        void check_consistency();

        void write_mesh(const std::string &title, mesh::Data_Visu data_visu = Data_Visu());

        int tag = -1;
        std::filesystem::path mesh_path;
        std::vector<int> verts_wall_num;
        std::vector<std::shared_ptr<Vertex>> verts;
        std::vector<std::shared_ptr<Edge>> edges;
        std::vector<std::shared_ptr<Face>> faces;
        std::vector<HalfEdge> halfEdges;
        std::vector<std::shared_ptr<Mesh>> subMesh;
    };

    void tokenize(std::string const &str, const char delim, std::vector<std::string> &out);
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> read_mesh(const std::filesystem::path &mesh_path);
} // namespace mesh
