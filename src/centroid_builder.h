
#ifndef _CENTROID_BVH_H
#define _CENTROID_BVH_H

#include <cassert>
#include <vector>
#include <algorithm>

#include "vec.h"
#include "ray.h"
#include "bvh.h"
#include "primitive_builder.h"


struct CentroidBuilder : public PrimitiveBuilder
{
    CentroidBuilder( ) = default;
    virtual ~CentroidBuilder( ) = default;
    
protected:
    int build_node( std::vector<Node>& nodes, const int begin, const int end )
    {
        // construit l'englobant des centres
        BBox cbounds= centroid_bounds(begin, end);
        
        // choisit l'axe le plus etire de l'englobant, par defaut
        int min_axis= cbounds.axis();
        float min_cost= 0;
        
    #if 1
        // compare le cout de la repartition sur chaque axe et garde le meilleur...
        min_cost= FLT_MAX;
        for(int axis= 0; axis < 3; axis++)
        {
            // position de la coupe / centre de l'englobant
            float cut= (cbounds.pmin(axis) + cbounds.pmax(axis)) / 2;
            
            // calcule les englobants des triangles dans la moitiee gauche et dans la moitiee droite
            int left_n= 0, right_n= 0;
            BBox left(PointMax, PointMin), right(PointMax, PointMin);
            for(int i= begin; i < end; i++)
            {
                if(primitives[i].centroid(axis) < cut)
                {
                    left.insert(primitives[i].bounds);
                    left_n++;
                }
                else
                {
                    right.insert(primitives[i].bounds);
                    right_n++;
                }
            }
            assert(left_n + right_n == end - begin);
            
            // evalue le cout de cette repartition
            float area= merge(left, right).area();
            float cost= 1 
                + left.area() / area * left_n 
                + right.area() / area * right_n;
            if(cost < min_cost)
            {
                // garde la meilleure repartition
                min_cost= cost;
                min_axis= axis;
            }
        }
    #endif
    
        // repartir les triangles 
        Primitive *p= std::partition(primitives.data() + begin, primitives.data() + end, centroid_less1(min_axis, (cbounds.pmin(min_axis) + cbounds.pmax(min_axis)) / 2));
        int m= std::distance(primitives.data(), p);
        if(m == begin || m == end 
        || min_cost > end - begin)
        {
            // construit une feuille si la partition est degeneree,
            // ou si repartir les triangles est moins interressant que de construire une feuille de plusieurs triangles...
            // quel est le cout de n triangles ?
            
            int index= int(nodes.size());
            nodes.push_back( make_leaf(bounds(begin, end), begin, end) );
            return index;
        }
        
        // repartir les triangles et construire les fils du noeud
        int left= build_node(nodes, begin, m);
        int right= build_node(nodes, m, end);
        
        // englobant du noeud : union des englobants des fils
        BBox bounds= merge(nodes[left].bounds, nodes[right].bounds);
        
        // construire le noeud
        int index= int(nodes.size());
        nodes.push_back( make_node(bounds, left, right) );
        return index;
    }
    
    // comparaison pour std::partition
    struct centroid_less1
    {
        int axis;
        float cut;
        
        centroid_less1( const int _axis, const float _cut ) : axis(_axis), cut(_cut) {}
        bool operator() ( const Primitive& p ) const { return p.centroid(axis) < cut; }
    };
};

#endif
