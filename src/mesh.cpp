#define _USE_MATH_DEFINES

#include <omp.h>
#include <set>

#include "mesh.h"

using namespace std;

namespace mesh
{
    std::vector<std::shared_ptr<Face>> Vertex::get_faces()
    {
        vector<shared_ptr<Face>> res_faces(this->faces.size());
        for (int f = 0; f < this->faces.size(); ++f)
            res_faces[f] = this->faces[f].lock();
        return res_faces;
    }

    inline std::ostream &operator<<(std::ostream &os, const Vertex &vert)
    {
        os << vert.coords;
        return os;
    }

    std::vector<std::shared_ptr<mesh::Face>> verts_faces_intersection(const std::shared_ptr<mesh::Vertex> &vert_a, const std::shared_ptr<mesh::Vertex> &vert_b)
    {
        vector<shared_ptr<Face>> res_faces;
        vector<shared_ptr<Face>> vert_a_faces = vert_a->get_faces();
        vector<shared_ptr<Face>> vert_b_faces = vert_b->get_faces();
        std::sort(vert_a_faces.begin(), vert_a_faces.end());
        std::sort(vert_b_faces.begin(), vert_b_faces.end());
        std::set_intersection(vert_a_faces.begin(), vert_a_faces.end(),
                              vert_b_faces.begin(), vert_b_faces.end(),
                              back_inserter(res_faces));
        return res_faces;
    }

    Edge::Edge(const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b)
    {
        this->verts = vector<shared_ptr<Vertex>>({a, b});
        this->lenght = (this->verts.front()->coords - this->verts.back()->coords).norm();
        assert(this->lenght > 1e-12);
    }

    std::vector<std::shared_ptr<Face>> Edge::get_faces()
    {
        vector<shared_ptr<Face>> faces;
        for (const auto &[tag, he] : this->tag2he)
        {
            if (he->face != nullptr)
                faces.push_back(he->face);
            if (he->twin->face != nullptr)
                faces.push_back(he->twin->face);
        }
        return faces;
    }

    std::ostream &operator<<(std::ostream &os, const Edge &edge)
    {
        os << (edge.verts.front()->coords) << " " << (edge.verts.back()->coords);
        return os;
    }

    Face::Face(const std::vector<std::shared_ptr<Vertex>> &verts)
    {
        this->verts = verts;
        //  double l_0 = (*this->verts[2]->coords - *this->verts[1]->coords).norm();
        //  double l_1 = (*this->verts[0]->coords - *this->verts[2]->coords).norm();
        //  double l_2 = (*this->verts[1]->coords - *this->verts[0]->coords).norm();
        //  double p_2 = (l_0 + l_1 + l_2) / 2;
        //   Formule de Heron... bizarerie cosmique lorsqu'un côté est très grand par rapport aux aux autres
        //   ou par exemple lorsque l_0=l_1+l_2 on trouve area=0
        //  this->area = sqrt(p_2 * (p_2 - l_0) * (p_2 - l_1) * (p_2 - l_2));
        //  cout << "verif " << l_0 << " " << l_1 << " " << l_2 << " " << p_2 << endl;
        //  cout << "verts " << *verts[0]->coords << " " << *verts[1]->coords << " " << *verts[2]->coords << endl;
        this->area = 0.5 * ((this->verts[1]->coords - this->verts[0]->coords) ^ (this->verts[2]->coords - this->verts[0]->coords)).norm();
        // cout << "area " << this->area << endl;
        assert(this->area > 1e-15);
        //  this->radius_circum = (l_0 * l_1 * l_2) / (4 * this->aire);
        //  assert(this->radius_circum > 1e-10);
        this->bary = ((this->verts[0]->coords + this->verts[1]->coords + this->verts[2]->coords) / 3);
    }

    pair<bool, Vect> Face::in_face(const Vect &pt)
    {
        return in_triangle(vector<Vect>({this->verts[0]->coords, this->verts[1]->coords, this->verts[2]->coords}), pt, this->normal);
    }

    Mesh::Mesh(const std::filesystem::path &mesh_path)
    {
        this->read_mesh(mesh_path);
        this->make_mesh();
    }

    //     const Mesh &Mesh ::operator=(const Mesh &mesh)
    //     {
    //         // map<HalfEdgeCIter, HalfEdgeIter, HalfEdgeCIterCompare> halfedgeOldToNew;
    //         // map<VertexCIter, VertexIter, VertexCIterCompare> vertexOldToNew;
    //         // map<EdgeCIter, EdgeIter, EdgeCIterCompare> edgeOldToNew;
    //         // map<FaceCIter, FaceIter, FaceCIterCompare> faceOldToNew;

    //         // // copy geometry from the original mesh and create a
    //         // // map from pointers in the original mesh to
    //         // // those in the new mesh
    //         // halfedges.clear();
    //         // for (HalfEdgeCIter i = mesh.halfedges.begin(); i != mesh.halfedges.end(); i++)
    //         //     halfedgeOldToNew[i] = halfedges.insert(halfedges.end(), *i);
    //         // vertices.clear();
    //         // for (VertexCIter i = mesh.vertices.begin(); i != mesh.vertices.end(); i++)
    //         //     vertexOldToNew[i] = vertices.insert(vertices.end(), *i);
    //         // edges.clear();
    //         // for (EdgeCIter i = mesh.edges.begin(); i != mesh.edges.end(); i++)
    //         //     edgeOldToNew[i] = edges.insert(edges.end(), *i);
    //         // faces.clear();
    //         // for (FaceCIter i = mesh.faces.begin(); i != mesh.faces.end(); i++)
    //         //     faceOldToNew[i] = faces.insert(faces.end(), *i);

    //         // // ``search and replace'' old pointers with new ones
    //         // for (HalfEdgeIter i = halfedges.begin(); i != halfedges.end(); i++)
    //         // {
    //         //     i->next = halfedgeOldToNew[i->next];
    //         //     i->flip = halfedgeOldToNew[i->flip];
    //         //     i->from = vertexOldToNew[i->from];
    //         //     i->edge = edgeOldToNew[i->edge];
    //         //     i->face = faceOldToNew[i->face];
    //         // }

    //         // for (VertexIter i = vertices.begin(); i != vertices.end(); i++)
    //         // {
    //         //     i->out = halfedgeOldToNew[i->out];
    //         // }

    //         // for (EdgeIter i = edges.begin(); i != edges.end(); i++)
    //         //     i->he = halfedgeOldToNew[i->he];
    //         // for (FaceIter i = faces.begin(); i != faces.end(); i++)
    //         //     i->he = halfedgeOldToNew[i->he];

    //         // return *this;
    //     }

    void Mesh::make_mesh()
    {
        cout << H_BLU << "Msg... " << NC << "Make mesh  " << endl;
        assert(this->verts.size() >= 3);
        assert(this->faces.size() >= 1);

        vector<vector<shared_ptr<Edge>>> vert_2_edges(this->verts.size());
        for (int e = 0; e < this->edges.size(); ++e)
        {
            // assert(this->edges[e]->is_wall == true);
            vert_2_edges[this->edges[e]->verts.front()->tag2num[this->tag]].push_back(this->edges[e]);
            vert_2_edges[this->edges[e]->verts.back()->tag2num[this->tag]].push_back(this->edges[e]);
        }

        vector<vector<shared_ptr<Face>>> edge_2_faces(3 * this->faces.size());

        for (int f = 0; f < this->faces.size(); ++f)
        {
            this->faces[f]->normal = (this->faces[f]->verts[1]->coords - this->faces[f]->verts[0]->coords) ^
                                     (this->faces[f]->verts[2]->coords - this->faces[f]->verts[0]->coords);
            this->faces[f]->normal.normalize();
            this->faces[f]->oX = this->faces[f]->verts[1]->coords - this->faces[f]->verts[0]->coords;
            this->faces[f]->oX.normalize();
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            this->faces[f]->edges.resize(3, nullptr);
            for (int k = 0; k < 3; ++k)
            {
                for (const auto &edge : vert_2_edges[this->faces[f]->verts[(k + 1) % 3]->tag2num[this->tag]])
                    if (edge->verts.front() == this->faces[f]->verts[(k + 2) % 3] ||
                        edge->verts.back() == this->faces[f]->verts[(k + 2) % 3])
                    {
                        this->faces[f]->edges[k] = edge;
                        break;
                    }
                if (this->faces[f]->edges[k] == nullptr)
                {
                    this->edges.push_back(make_shared<Edge>(this->faces[f]->verts[(k + 1) % 3], this->faces[f]->verts[(k + 2) % 3]));
                    this->edges.back()->tag2num[this->tag] = this->edges.size() - 1;
                    this->faces[f]->edges[k] = this->edges.back();
                    /*this->edges.back()->verts.front()->edges.push_back(this->edges.back());
                    this->edges.back()->verts.back()->edges.push_back(this->edges.back());
                    this->edges.back()->verts.front()->verts.push_back(this->edges.back()->verts.back());
                    this->edges.back()->verts.back()->verts.push_back(this->edges.back()->verts.front());*/
                    vert_2_edges[this->faces[f]->verts[(k + 1) % 3]->tag2num[this->tag]].push_back(this->edges.back());
                    vert_2_edges[this->faces[f]->verts[(k + 2) % 3]->tag2num[this->tag]].push_back(this->edges.back());
                }
                edge_2_faces[this->faces[f]->edges[k]->tag2num.at(this->tag)].push_back(this->faces[f]);
                // this->faces[f]->edges[k]->faces.push_back(this->faces[f]);
                // this->faces[f]->verts[k]->faces.push_back(this->faces[f]);
            }
        }

        set<shared_ptr<Edge>> edges_wall_set;
        for (int v = 0; v < this->verts_wall_num.size(); ++v)
        {
            assert(verts_wall_num[v] >= 0 && verts_wall_num[v] < this->verts.size());
            edges_wall_set.insert(vert_2_edges[verts_wall_num[v]].begin(), vert_2_edges[verts_wall_num[v]].end());
        }
        for (const auto &edge : edges_wall_set)
            edge->is_wall = true;

        // for (int e = 0; e < edges.size(); ++e)
        // {
        //     // cout << "e " << e << endl;
        //     if (edge_2_faces[e].size() != 2)
        //     {
        //         this->edges[e]->is_wall = true;
        //     }
        // }

        int max_face_tag = -2;
        // cout << "Tag---";
        {
            // cout << H_BLU << "Msg... " << NC << " tag_face  " << endl;
            vector<shared_ptr<Face>> container(this->faces.size(), nullptr);
            int deb = 0;
            int fin = deb;
            int mark = 0;
            int nb = 1;
            this->faces[0]->tag = mark;
            container[deb] = this->faces[0];

            while (nb > 0)
            {
                while (nb != 0)
                {
                    nb = 0;
                    for (int t = deb; t <= fin; ++t)
                    {
                        for (int k = 0; k < 3; ++k)
                        {
                            if (container[t]->edges[k]->is_wall == false)
                            {
                                shared_ptr<Face> adj_face = nullptr;
                                if (edge_2_faces[container[t]->edges[k]->tag2num[this->tag]].front() == container[t])
                                    adj_face = edge_2_faces[container[t]->edges[k]->tag2num[this->tag]].back();
                                else if (edge_2_faces[container[t]->edges[k]->tag2num[this->tag]].back() == container[t])
                                    adj_face = edge_2_faces[container[t]->edges[k]->tag2num[this->tag]].front();
                                else
                                    assert(false);
                                if (adj_face->tag == -1)
                                {
                                    nb++;
                                    adj_face->tag = mark;
                                    container[fin + nb] = adj_face;
                                }
                            }
                        }
                    }
                    deb = fin;
                    fin = fin + nb;
                }

                assert(nb == 0);
                for (int f = 0; f < this->faces.size(); ++f)
                    if (this->faces[f]->tag == -1)
                    {
                        assert(deb == fin);
                        deb++;
                        fin = deb;
                        mark = mark + 1;
                        nb = 1;
                        this->faces[f]->tag = mark;
                        container[deb] = this->faces[f];
                        break;
                    }
            }
            max_face_tag = mark;
        }
        assert(max_face_tag >= 0);

        // cout << "Manifold mesh---";
        // cout << H_BLU << "Msg... " << NC << " build_manifold_mesh  " << endl;
        this->halfEdges.reserve(2 * this->edges.size() + this->faces.size());
        this->halfEdges.resize(3 * this->faces.size());
        // assert(max_face_tag != 0);
        this->subMesh.resize(max_face_tag + 1);
        for (int s = 0; s < this->subMesh.size(); ++s)
        {
            this->subMesh[s] = make_shared<Mesh>();
            this->subMesh[s]->tag = s;
        }

        vector<int> incr_per_tag(this->subMesh.size(), -1);
        for (int f = 0; f < this->faces.size(); ++f)
        {
            this->subMesh[this->faces[f]->tag]->faces.push_back(this->faces[f]);
            this->faces[f]->tag2num[this->faces[f]->tag] = this->subMesh[this->faces[f]->tag]->faces.size() - 1;
            for (int k = 0; k < 3; ++k)
            {
                this->halfEdges[3 * f + k].prev = this->halfEdges.begin() + 3 * f + (k + 2) % 3; // ((((3 * f + k - 1) - 3 * f) % 3) + 3) % 3 + 3 * f;
                this->halfEdges[3 * f + k].next = this->halfEdges.begin() + 3 * f + (k + 1) % 3; // ((((3 * f + k + 1) - 3 * f) % 3) + 3) % 3 + 3 * f;
                this->halfEdges[3 * f + k].origin = this->faces[f]->verts[(k + 1) % 3];
                this->halfEdges[3 * f + k].target = this->faces[f]->verts[(k + 2) % 3];
                this->halfEdges[3 * f + k].egde = this->faces[f]->edges[k];
                this->halfEdges[3 * f + k].face = this->faces[f];
                incr_per_tag[this->faces[f]->tag]++;
                this->halfEdges[3 * f + k].num = incr_per_tag[this->faces[f]->tag];

                if (this->faces[f]->verts[k]->tag2num.find(this->faces[f]->tag) == this->faces[f]->verts[k]->tag2num.end())
                {
                    this->subMesh[this->faces[f]->tag]->verts.push_back(this->faces[f]->verts[k]);
                    this->faces[f]->verts[k]->tag2num[this->faces[f]->tag] = this->subMesh[this->faces[f]->tag]->verts.size() - 1;
                    this->faces[f]->verts[k]->tag2he[this->faces[f]->tag] = this->halfEdges[3 * f + k].prev;
                }

                if (this->faces[f]->edges[k]->tag2num.find(this->faces[f]->tag) == this->faces[f]->edges[k]->tag2num.end())
                {
                    this->subMesh[this->faces[f]->tag]->edges.push_back(this->faces[f]->edges[k]);
                    this->faces[f]->edges[k]->tag2num[this->faces[f]->tag] = this->subMesh[this->faces[f]->tag]->edges.size() - 1;
                    this->faces[f]->edges[k]->tag2he[this->faces[f]->tag] = this->halfEdges.begin() + 3 * f + k;

                    if (edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].size() != 2 ||
                        (edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].size() == 2 &&
                         edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].front()->tag != edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].back()->tag))
                    {
                        this->halfEdges.push_back(HalfEdge());
                        incr_per_tag[this->faces[f]->tag]++;
                        this->halfEdges.back().num = incr_per_tag[this->faces[f]->tag];
                        this->halfEdges.back().twin = this->halfEdges.begin() + 3 * f + k;
                        this->halfEdges.back().twin->twin = this->halfEdges.end() - 1;
                        assert(this->halfEdges.end() - 1 == this->halfEdges.back().twin->twin);
                        this->halfEdges.back().origin = this->halfEdges.back().twin->target;
                        this->halfEdges.back().target = this->halfEdges.back().twin->origin;
                        this->halfEdges.back().egde = this->faces[f]->edges[k];
                    }
                }
                else if (edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].size() == 2 &&
                         edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].front()->tag == edge_2_faces[this->faces[f]->edges[k]->tag2num[this->tag]].back()->tag)
                {
                    this->halfEdges[3 * f + k].twin = this->faces[f]->edges[k]->tag2he.at(this->faces[f]->tag);
                    this->halfEdges[3 * f + k].twin->twin = this->halfEdges.begin() + 3 * f + k;
                    if (!(this->halfEdges[3 * f + k].origin == this->halfEdges[3 * f + k].twin->target &&
                          this->halfEdges[3 * f + k].target == this->halfEdges[3 * f + k].twin->origin))
                    {
                        cerr << RED << "Problème d'orientation dans le maillage input" << endl;
                        cerr << "Certains triangles du maillage ne sont pas orienté dans le même sens" << NC << endl;
                        throw std::logic_error("Orientation problem.");
                    }
                }
                else
                    assert(false);
            }
            this->faces[f]->he = this->halfEdges.begin() + 3 * f;
        }
        assert(this->halfEdges.size() <= 2 * this->edges.size() + this->faces.size());

// this->vertOnBound.resize(this->verts.size(), 0);
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            for (int h = 0; h < this->halfEdges.size(); ++h)
            {
                if (this->halfEdges[h].face == nullptr)
                {
                    this->halfEdges[h].egde->is_wall = true;
                    this->halfEdges[h].egde->verts.front()->on_bound = true;
                    this->halfEdges[h].egde->verts.back()->on_bound = true;
                    // this->vertOnBound[this->halfEdges[h].origin->num(this->tag)] = 1;
                    this->halfEdges[h].twin->origin->tag2he[this->halfEdges[h].twin->face->tag] = (this->halfEdges.begin() + h)->twin;
                    auto next = this->halfEdges[h].twin;
                    do
                    {
                        next = next->prev->twin;
                    } while (next->face != nullptr);
                    this->halfEdges[h].next = next;
                    next->prev = this->halfEdges.begin() + h;
                }
            }
        }

        // cout << "nb subTh " << this->subMesh.size() << endl;

#pragma omp parallel
        {
#pragma omp for schedule(auto)
            for (int v = 0; v < this->verts.size(); ++v)
            {
                this->verts[v]->tag2oZ[this->tag] = Vect(0, 0, 0);
                double big_sum_area = 0;
                for (const auto &[tag, he] : this->verts[v]->tag2he)
                {
                    double sum_area = 0;
                    double sum_ang = 0;
                    this->verts[v]->tag2oZ[tag] = Vect(0, 0, 0);
                    auto hee = he;
                    do
                    {
                        if (hee->face != nullptr)
                            this->verts[v]->faces.push_back(hee->face);
                        //////////////////////////////////////
                        this->verts[v]->tag2oZ[tag] += hee->face->normal * hee->face->area;
                        this->verts[v]->tag2oZ[this->tag] += hee->face->normal * hee->face->area;
                        sum_area += hee->face->area;
                        big_sum_area += hee->face->area;
                        double lA = hee->egde->lenght;
                        double lOpp = hee->next->egde->lenght;
                        double lB = hee->next->next->egde->lenght;
                        double temp = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
                        temp = temp > 1.0 ? 1.0 : temp < -1.0 ? -1.0
                                                              : temp; // Clamp between [-1,1]
                        hee->ang_corner = std::acos(temp);
                        sum_ang += hee->ang_corner;
                        hee = hee->prev->twin;
                    } while (hee != he && hee->face != nullptr);
                    this->verts[v]->tag2oZ[tag] = this->verts[v]->tag2oZ[tag] * (1 / sum_area);
                    this->verts[v]->tag2oZ[tag].normalize();

                    hee = he;
                    double temp_ang = 0;
                    this->verts[v]->tag2oX[tag] = Vect(0, 0, 0);
                    do
                    {
                        hee->localVec = hee->egde->lenght * std::complex<double>(cos(temp_ang), sin(temp_ang));
                        if (this->verts[v]->on_bound == true)
                            temp_ang += (M_PI / sum_ang) * hee->ang_corner;
                        else
                            temp_ang += (2 * M_PI / sum_ang) * hee->ang_corner;

                        this->verts[v]->tag2oX[tag] += (hee->target->coords - hee->origin->coords)
                                                           .proj(this->verts[v]->tag2oZ.at(tag))
                                                           .rotate(this->verts[v]->tag2oZ.at(tag), -atan2(hee->localVec.imag(), hee->localVec.real()));

                        hee = hee->prev->twin;
                    } while (hee != he);
                    this->verts[v]->tag2oX[tag].normalize();

                    assert(this->verts[v] == he->origin);
                }
                this->verts[v]->tag2oZ[this->tag] = this->verts[v]->tag2oZ[this->tag] * (1 / big_sum_area);
                this->verts[v]->tag2oZ[this->tag].normalize();
            }

#pragma omp for schedule(auto)
            for (int h = 0; h < this->halfEdges.size(); ++h)
            {
                auto he = this->halfEdges.begin() + h;
                assert(he->num != -1);
                assert(he == he->twin->twin);
                assert(he->face == he->next->face);
                assert(he->face == he->prev->face);
                assert(he == he->prev->next);
                assert(he == he->next->prev);
                assert(he->origin != nullptr);
                assert(he->egde != nullptr);
            }

#pragma omp for schedule(auto)
            for (int e = 0; e < this->edges.size(); ++e)
            {
                for (const auto &[tag, he] : this->edges[e]->tag2he)
                {
                    assert(he == he->twin->twin);
                    if (he->origin == this->edges[e]->verts.front())
                        assert(he->twin->origin == this->edges[e]->verts.back());
                    else if (he->origin == this->edges[e]->verts.back())
                        assert(he->twin->origin == this->edges[e]->verts.front());
                    else
                        assert(false);
                }
            }

#pragma omp for schedule(auto)
            for (int f = 0; f < this->faces.size(); ++f)
                assert(this->faces[f] == this->faces[f]->he->face);
        }
    }
} // namespace mesh
