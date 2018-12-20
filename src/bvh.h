
#ifndef _BVH_H
#define _BVH_H

#include <cassert>

#include <vector>
#include <algorithm>
#include <chrono>

#include "vec.h"
#include "ray.h"
#include "mesh.h"


// utilitaire : renvoie les coordonnees min / max de 2 points / vecteurs
static inline 
Point min( const Point& a, const Point& b ) { return Point( std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) ); }
static inline 
Point max( const Point& a, const Point& b ) { return Point( std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) ); }


struct BBox
{
    Point pmin;
    Point pmax;
    
    BBox( ) = default;
    
    BBox( const Point& p ) : pmin(p), pmax(p) {}
    BBox( const Point& _pmin, const Point& _pmax ) : pmin(_pmin), pmax(_pmax) {}
    BBox( const BBox& b ) : pmin(b.pmin), pmax(b.pmax) {}
    BBox& operator= ( const BBox& b ) { pmin= b.pmin; pmax= b.pmax; return *this; }
    
    BBox& insert( const Point& p ) { pmin= min(pmin, p); pmax= max(pmax, p); return *this; }
    BBox& insert( const BBox& b ) { pmin= min(pmin, b.pmin); pmax= max(pmax, b.pmax); return *this; }
    
    Point centroid( ) const { return (pmin+pmax)/2; }
    
    bool empty( ) const { Vector d(pmin, pmax); return d.x >= 0 && d.y >= 0 && d.z >= 0; }
    
    float area( ) const { Vector d(pmin, pmax); return 2 * d.x*d.y + 2 * d.x*d.z + 2 * d.y*d.z; }
    
    int axis( ) const 
    { 
        Vector d(pmin, pmax);
        if(d.x > d.y && d.x > d.z) return 0;
        else if(d.y > d.z) return 1;
        else  return 2;
    }
};

static inline
BBox merge( const BBox& a, const BBox& b ) { return BBox(min(a.pmin, b.pmin), max(a.pmax, b.pmax)); }
static inline
BBox overlap( const BBox& a, const BBox& b ) { return BBox(max(a.pmin, b.pmin), min(a.pmax, b.pmax)); }


struct NodeHit
{
    float tmin, tmax;
    
    NodeHit( ) : tmin(0), tmax(-1) {}   // pas d'intersection
    NodeHit( const float _tmin, const float _tmax ) : tmin(_tmin), tmax(_tmax) {}
    
    operator bool( ) const { return (tmin <= tmax); }      // renvoie vrai si l'intersection est initialisee...
};

struct Node
{
    BBox bounds;
    int left, right;
    
    Node( const BBox& _bounds, const int _v0, const int _v1 ) : bounds(_bounds), left(_v0), right(_v1) {}
    
    //! verifie que le noeud est une feuille.
    bool leaf( ) const { return right < 0; }
    //! renvoie le premier triangle reference par la feuille.
    int begin( ) const { assert(leaf()); return left; }
    //! renvoie le dernier triangle reference par la feuille.
    int end( ) const { assert(leaf()); return left - right; }
    
    //! verifie que le noeud est un noeud interne.
    bool node( ) const { return !leaf(); }
    
    NodeHit intersect( const Ray& ray, const Vector& invd, const float htmax ) const
    {
        Point rmin= bounds.pmin;
        Point rmax= bounds.pmax;
        if(ray.d.x < 0) std::swap(rmin.x, rmax.x);
        if(ray.d.y < 0) std::swap(rmin.y, rmax.y);
        if(ray.d.z < 0) std::swap(rmin.z, rmax.z);
        Vector dmin= (rmin - ray.o) * invd;
        Vector dmax= (rmax - ray.o) * invd;
        
        float tmin= std::max(dmin.z, std::max(dmin.y, std::max(dmin.x, 0.f)));
        float tmax= std::min(dmax.z, std::min(dmax.y, std::min(dmax.x, htmax)));
        return NodeHit(tmin, tmax);
    }
};

//! constructeur nomme. initialise un noeud.
static inline
Node make_node( const BBox& bounds, const int left, const int right )
{
    Node tmp= Node(bounds, left, right);
    assert(tmp.node()); 
    return tmp;
}

//! constructeur nomme. initialise une feuille.
static inline
Node make_leaf( const BBox& bounds, const int begin, const int end )
{
    Node tmp= Node(bounds, begin, - (end - begin));
    assert(tmp.leaf()); 
    return tmp;
}


struct Triangle
{
    Point p;
    Vector e1, e2;
    int id;
    
    Triangle( const Point& _a, const Point& _b, const Point& _c, const int _id ) : p(_a), e1(Vector(_a, _b)), e2(Vector(_a, _c)), id(_id) {}
    
    //! renvoie l'englobant du triangle.
    BBox bounds( ) const { return BBox(p).insert(p+e1).insert(p+e2); }
    
    /* calcule l'intersection ray/triangle
        cf "fast, minimum storage ray-triangle intersection" 
        http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
        
        renvoie faux s'il n'y a pas d'intersection valide (une intersection peut exister mais peut ne pas se trouver dans l'intervalle [0 htmax] du rayon.)
        renvoie vrai + les coordonnees barycentriques (u, v) du point d'intersection + sa position le long du rayon (t).
        convention barycentrique : p(u, v)= (1 - u - v) * a + u * b + v * c
    */
    Hit intersect( const Ray &ray, const float htmax ) const
    {
        Vector pvec= cross(ray.d, e2);
        float det= dot(e1, pvec);
        
        float inv_det= 1 / det;
        Vector tvec(p, ray.o);

        float u= dot(tvec, pvec) * inv_det;
        if(u < 0 || u > 1) return Hit();

        Vector qvec= cross(tvec, e1);
        float v= dot(ray.d, qvec) * inv_det;
        if(v < 0 || u + v > 1) return Hit();

        float t= dot(e2, qvec) * inv_det;
        if(t > htmax || t < 0) return Hit();
        
        return Hit(id, t, u, v);           // p(u, v)= (1 - u - v) * a + u * b + v * c
    }
};



struct Builder
{
    Builder( ) = default;
    virtual ~Builder( ) = default;
    
    virtual int operator() ( std::vector<Triangle>& triangles, std::vector<Node>& nodes ) = 0;
};


struct BVH
{
    std::vector<Triangle> triangles;
    std::vector<Node> nodes;
    int root;
    
    BVH( ) : triangles(), nodes(), root(-1) {}
    BVH( const Mesh& mesh, Builder& builder ) : triangles(), nodes(), root(-1) { build(mesh, builder); }
    
    void build( const Mesh& mesh, Builder& builder )
    {
        triangles.clear();
        nodes.clear();
        
        auto cpu_start= std::chrono::high_resolution_clock::now();
        
        // recupere les triangles
        for(int i= 0; i < mesh.triangle_count(); i++)
        {
            TriangleData triangle= mesh.triangle(i);
            triangles.push_back( Triangle(triangle.a, triangle.b, triangle.c, i) );
        }
        
        // construit l'arbre
        root= builder(triangles, nodes);
        
        auto cpu_stop= std::chrono::high_resolution_clock::now();
        int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
        printf("bvh %d nodes, %d triangles, root %d\n", int(nodes.size()), int(triangles.size()), root);
        printf("  cpu  %ds %03dms\n", int(cpu_time / 1000), int(cpu_time % 1000));
    }
    
    //! renvoie vrai si une intersection valide existe. la position de l'intersection *la plus proche* de l'origine du rayon est renvoyee dans hit.
    Hit intersect( const Ray& ray ) const
    {
        assert(root != -1);
        
        Vector invd= Vector(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        Hit hit;
        hit.t= ray.tmax;
        intersect(root, ray, invd, hit);
        return hit;
    }
    
    //! renvoie vrai si une intersection valide existe. plus rapide que de calculer la plus proche de l'origine du rayon.
    bool visible( const Ray& ray ) const
    {
        assert(root != -1);
        
        Vector invd= Vector(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
        return visible(root, ray, invd);
    }    
    
protected:
    void intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const
    {
        const Node& node= nodes[index];
        if(node.leaf())
        {
            // feuille
            for(int i= node.begin(); i < node.end(); i++)
            {
                // ne renvoie vrai que si l'intersection existe dans l'intervalle [0 tmax]
                if(Hit h= triangles[i].intersect(ray, hit.t))
                    hit= h;
            }
        }
        else 
        {
            // noeud interne
        #if 0
            // parcours simple
            if(node.intersect(ray, invd, hit.t))
            {
                intersect(node.left, ray, invd, hit);
                intersect(node.right, ray, invd, hit);
            }
            
        #else
            // parcours ordonne
            NodeHit left= nodes[node.left].intersect(ray, invd, hit.t);
            NodeHit right= nodes[node.right].intersect(ray, invd, hit.t);
            
            if(left && right)
            {
                // les 2 fils sont touches par le rayon
                if(left.tmin < right.tmin)
                {
                    // le fils gauche est plus pres
                    intersect(node.left, ray, invd, hit);
                    if(hit.t >= right.tmin)
                        // pacourir le fils droit, si necessaire
                        intersect(node.right, ray, invd, hit);
                }
                else
                {
                    intersect(node.right, ray, invd, hit);
                    if(hit.t >= left.tmin)
                        intersect(node.left, ray, invd, hit);
                }
            }
            
            // 1 seul fils est touche par le rayon
            else if(left)
                intersect(node.left, ray, invd, hit);
            else if(right)
                intersect(node.right, ray, invd, hit);
        #endif
        }
    }
    
    bool visible( const int index, const Ray& ray, const Vector& invd ) const
    {
        const Node& node= nodes[index];
        if(node.leaf())
        {
            for(int i= node.begin(); i < node.end(); i++)
                // ne renvoie vrai que si l'intersection existe dans l'intervalle [0 tmax]
                if(triangles[i].intersect(ray, ray.tmax))
                    return false;
        }
        else if(node.intersect(ray, invd, ray.tmax))
        {
            if(visible(node.left, ray, invd) == false)
                return false;
            if(visible(node.right, ray, invd) == false)
                return false;
        }
        
        return true;
    }    
};

#endif
