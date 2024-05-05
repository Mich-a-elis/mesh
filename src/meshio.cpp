#include <iostream>
#include <vector>
#include <filesystem>
#include <string>
#include <cassert>
#include <fstream>
#include <tuple>

#include <omp.h>
#include <set>

#include "mesh.h"

using namespace std;

namespace mesh
{
    void tokenize(std::string const &str, const char delim, std::vector<std::string> &out)
    {
        size_t start;
        size_t end = 0;

        while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
        {
            end = str.find(delim, start);
            out.push_back(str.substr(start, end - start));
        }
    }

    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> read_mesh(const std::filesystem::path &mesh_path)
    {
        std::cout << "Reading mesh: " << mesh_path << std::endl;

        assert(std::filesystem::exists(mesh_path) ? true : (std::cerr << "Error: couldn't open file " << mesh_path << " for input.\n", false));

        std::vector<double> verts_coords;
        std::vector<int> wall;
        std::vector<int> face_indVerts;

        std::vector<std::string> list_lines;
        // lecture complete du fichier
        {
            std::ifstream file(mesh_path);
            assert(file.is_open() ? true : (std::cerr << "Error: couldn't open file " << mesh_path << " for input.\n", false));
            std::string contents;
            file.seekg(0, std::ios::end);
            contents.resize(file.tellg());
            file.seekg(0, std::ios::beg);
            file.read(&contents[0], contents.size());
            const char delim = '\n';
            tokenize(contents, delim, list_lines);
            file.close();
        }

        if (mesh_path.extension().string() == ".msh")
        {
            int ind_line = 0;
            std::string temp_string;

            // lecture des points
            ind_line++;
            double version;
            int file_type, data_size;
            std::istringstream(list_lines[ind_line]) >> version >> file_type >> data_size;
            assert(version == 2.2);
            assert(ind_line < list_lines.size());

            do
            {
                ind_line++;
                std::istringstream(list_lines[ind_line]) >> temp_string;
            } while (ind_line < list_lines.size() && temp_string != "$Nodes");
            assert(ind_line < list_lines.size() && temp_string == "$Nodes");

            ind_line++;
            int nv = 0;
            std::istringstream(list_lines[ind_line]) >> nv;
            verts_coords.resize(3 * nv);
            ind_line++;
#pragma omp parallel // num_threads(4)
            {
                int iNode;
                double x, y, z;
                std::istringstream line;
#pragma omp for schedule(auto)
                for (int v = 0; v < nv; ++v)
                {
                    assert(list_lines[v + ind_line].size() != 0);
                    line.clear();
                    line.str(list_lines[v + ind_line]);
                    line >> iNode >> x >> y >> z;
                    // this->verts[v] = make_shared<Vertex>();
                    // this->verts[v]->coords = Vect(x, y, z);
                    // this->verts[v]->tag2num[this->tag] = v;
                    verts_coords[3 * v] = x;
                    verts_coords[3 * v + 1] = y;
                    verts_coords[3 * v + 2] = z;
                }
            }
            ind_line += verts_coords.size() / 3 - 1;
            ind_line += 2;
            if (ind_line >= list_lines.size())
                return std::tuple(verts_coords, face_indVerts, wall);
            std::istringstream(list_lines[ind_line]) >> temp_string;
            assert(temp_string == "$Elements");

            int nbElements = 0;
            ind_line++;
            std::istringstream(list_lines[ind_line]) >> nbElements;

            if (nbElements <= 0)
                return std::tuple(verts_coords, face_indVerts, wall);

            {
                std::istringstream line;
                int iElt, elt_type, nbFlags, flag, un_used, a, b, c;

                ind_line++;
                assert(list_lines[ind_line].size() != 0);
                line.clear();
                line.str(list_lines[ind_line]);
                line >> iElt >> elt_type >> nbFlags;
                while (elt_type == 1)
                {
                    assert(nbFlags == 2);
                    line >> flag >> un_used >> a >> b;
                    // this->edges.push_back(make_shared<Edge>(this->verts[a - 1], this->verts[b - 1]));
                    // this->edges.back()->tag2num[this->tag] = this->edges.size() - 1;
                    // this->edges.back()->is_wall = true;
                    wall.push_back(a - 1);
                    wall.push_back(b - 1);
                    ind_line++;
                    if (ind_line >= list_lines.size())
                        break;
                    assert(list_lines[ind_line].size() != 0);
                    line.clear();
                    line.str(list_lines[ind_line]);
                    line >> iElt >> elt_type >> nbFlags;
                }
            }

            int nbFace = nbElements - wall.size() / 2;
            face_indVerts.resize(3 * nbFace);
            // #pragma omp parallel // num_threads(4)
            {
                std::istringstream line;
                int iElt, elt_type, nbFlags, flag, un_used, a, b, c;
                // #pragma omp for schedule(auto)
                for (int f = 0; f < nbFace; ++f)
                {
                    assert(list_lines[ind_line + f].size() != 0);
                    line.clear();
                    line.str(list_lines[ind_line + f]);
                    line >> iElt >> elt_type >> nbFlags;
                    assert(elt_type == 2);
                    assert(nbFlags == 2);
                    line >> flag >> un_used >> a >> b >> c;
                    // this->faces[f] = std::make_shared<Face>(std::vector<shared_ptr<Vertex>>({this->verts[a - 1],
                    //                                                                          this->verts[b - 1],
                    //                                                                          this->verts[c - 1]}));
                    // this->faces[f]->tag2num[this->tag] = this->faces.size() - 1;
                    // this->faces[f]->flag = flag;
                    face_indVerts[3 * f] = a - 1;
                    face_indVerts[3 * f + 1] = b - 1;
                    face_indVerts[3 * f + 2] = c - 1;
                }
            }
        }
        else if (mesh_path.extension().string() == ".obj")
        {
            std::string unused;
            std::istringstream line;
            double x, y, z;
            int a, b, c;
            for (int k = 0; k < list_lines.size(); ++k)
            {
                line.clear();
                line.str(list_lines[k]);
                line >> unused;
                if (unused == "v")
                {
                    line >> x >> y >> z;
                    verts_coords.push_back(x);
                    verts_coords.push_back(y);
                    verts_coords.push_back(z);
                }
                else if (unused == "f")
                {
                    line >> a >> b >> c;
                    face_indVerts.push_back(a - 1);
                    face_indVerts.push_back(b - 1);
                    face_indVerts.push_back(c - 1);
                }
            }
        }
        else
            assert(false && "Error: input file format\n");

        return std::tuple(verts_coords, face_indVerts, wall);
    }
}