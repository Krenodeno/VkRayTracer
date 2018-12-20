
#ifndef _PRIMITIVE_BVH_H
#define _PRIMITIVE_BVH_H

#include <cassert>
#include <vector>
#include <algorithm>

#include "vec.h"
#include "ray.h"
#include "bvh.h"


struct Primitive
{
    BBox bounds;
    Point centroid;
    int id;
    
    Primitive( const BBox& _bounds, const int _id ) : bounds(_bounds), centroid(_bounds.centroid()), id(_id) {}
};

static const Point PointMax(FLT_MAX, FLT_MAX, FLT_MAX);
static const Point PointMin(-FLT_MAX, -FLT_MAX, -FLT_MAX);

 
struct PrimitiveBuilder : public Builder
{
    PrimitiveBuilder( ) = default;
    virtual ~PrimitiveBuilder( ) = default;
    
    virtual int operator() ( std::vector<Triangle>& triangles, std::vector<Node>& nodes )
    {
        // construit les primitives a partir des triangles
        primitives.reserve(triangles.size());
        for(int i= 0; i < int(triangles.size()); i++)
            primitives.push_back( Primitive(triangles[i].bounds(), i) );
        
        // trie les primitives et construit les noeuds
        int root= build_node(nodes, 0, int(primitives.size()));
        
        // copie les triangles dans le meme ordre que les primitives
        std::vector<Triangle> remap;
        remap.reserve(primitives.size());
        for(int i= 0; i < int(primitives.size()); i++)
        {
            assert(primitives[i].id < int(triangles.size()));
            remap.push_back( triangles[primitives[i].id] );
        }
        triangles.swap(remap);
        
        // renvoie l'indice de la racine
        return root;
    }
    
protected:
    std::vector<Primitive> primitives;
    
    virtual int build_node( std::vector<Node>& nodes, const int begin, const int end ) = 0;
    
    BBox bounds( const int begin, const int end )
    {
        BBox bounds= primitives[begin].bounds;
        for(int i= begin +1; i < end; i++)
            bounds.insert(primitives[i].bounds);
            
        return bounds;
    }
    
    BBox centroid_bounds( const int begin, const int end )
    {
        BBox bounds= primitives[begin].centroid;
        for(int i= begin+1; i < end; i++)
            bounds.insert(primitives[i].centroid);
            
        return bounds;
    }
};

#endif
