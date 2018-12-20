
#include <cfloat>
#include <cassert>

#include <vector>
#include <algorithm>
#include <chrono>

#include "vec.h"
#include "mesh.h"
#include "wavefront.h"


struct Ray
{
    Point o;
    Vector d;
    float tmax;
    
    Ray( const Point& _o, const Point& _e ) : o(_o), d(Vector(_o, _e)), tmax(1) {}
    Ray( const Point& _o, const Vector& _d ) : o(_o), d(_d), tmax(FLT_MAX) {}
};

struct Hit
{
    int triangle_id;
    float t;
    float u, v;
    
    Hit( ) : triangle_id(-1), t(0), u(0), v(0) {}       // pas d'intersection
    Hit( const int _id, const float _t, const float _u, const float _v ) : triangle_id(_id), t(_t), u(_u), v(_v) {}
    
    operator bool( ) const { return (triangle_id != -1); }      // renvoie vrai si l'intersection est initialisee...
};


struct NodeHit
{
    float tmin, tmax;
    
    NodeHit( ) : tmin(0), tmax(-1) {}   // pas d'intersection
    NodeHit( const float _tmin, const float _tmax ) : tmin(_tmin), tmax(_tmax) {}
    
    operator bool( ) const { return (tmin <= tmax); }      // renvoie vrai si l'intersection est initialisee...
};

struct Node
{
    Point pmin;
    Point pmax;
    int left, right;
    
    NodeHit intersect( const Ray& ray, const Vector& invd, const float htmax ) const
    {
        Point rmin= pmin;
        Point rmax= pmax;
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

// utilitaire : renvoie les coordonnees min / max de 2 points / vecteurs
inline Point min( const Point& a, const Point& b ) { return Point( std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) ); }
inline Point max( const Point& a, const Point& b ) { return Point( std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) ); }

struct Triangle
{
    Point p;
    Vector e1, e2;
    int id;
    
    Triangle( const Point& _a, const Point& _b, const Point& _c, const int _id ) : p(_a), e1(Vector(_a, _b)), e2(Vector(_a, _c)), id(_id) {}
    
    //! renvoie l'englobant du triangle.
    void bounds( Point& pmin, Point& pmax ) const
    {
        pmin= p;
        pmin= min(pmin, p+e1);
        pmin= min(pmin, p+e2);
        
        pmax= p;
        pmax= max(pmax, p+e1);
        pmax= max(pmax, p+e2);
    }
    
    //! renvoye le centre de l'englobant du tiangle.
    Point centroid( ) const
    {
        Point pmin, pmax;
        bounds(pmin, pmax);
        return (pmin+pmax)/2;
    }
    
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


struct BVH
{
    std::vector<Triangle> triangles;
    std::vector<Node> nodes;
    int root;
    
    BVH( ) : triangles(), nodes(), root(-1) {}
    
    //! construit l'arbre avec les triangles de mesh.
    void build( const Mesh& mesh )
    {
        auto cpu_start= std::chrono::high_resolution_clock::now();
        
        // recupere les triangles
        for(int i= 0; i < mesh.triangle_count(); i++)
        {
            TriangleData triangle= mesh.triangle(i);
            triangles.push_back( Triangle(triangle.a, triangle.b, triangle.c, i) );
        }
        
        // construit l'arbre
        //~ root= build_node(0, triangles.size());
        root= build_node_centroids(0, triangles.size());
        
        auto cpu_stop= std::chrono::high_resolution_clock::now();
        int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
        printf("cpu  %ds %03dms\n", int(cpu_time / 1000), int(cpu_time % 1000));
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

protected:
    void intersect( const int index, const Ray& ray, const Vector& invd, Hit& hit ) const
    {
        const Node& node= nodes[index];
        if(node.right < 0)
        {
            // feuille
            int begin= -node.left;
            int end= -node.right;
            for(int i= begin; i < end; i++)
            {
                // ne renvoie vrai que si l'intersection existe dans l'intervalle [0 tmax]
                if(Hit h= triangles[i].intersect(ray, hit.t))
                    hit= h;
            }
        }
        else 
            // noeud interne
            if(node.intersect(ray, invd, hit.t))
            {
            #if 0
                // parcours simple
                intersect(node.left, ray, invd, hit);
                intersect(node.right, ray, invd, hit);
            #else
                // parcours ordonne
                NodeHit left= nodes[node.left].intersect(ray, invd, hit.t);
                NodeHit right= nodes[node.right].intersect(ray, invd, hit.t);
                
                if(left && right)
                {
                    // les 2 fils sont touches par le rayon
                    // rappel: visiter les feuilles en s'eloignant de l'origine du rayon.
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
    
    int build_node( const int begin, const int end )
    {
        // solution 2 : trier les triangles sur l'axe le plus etire de l'englobant
            // construire l'englobant des 2 moities des triangles
            // 1. construire l'englobant
            // 2. trouver l'axe le plus etire de l'englobant
            // 3. trier les triangles sur l'axe
            // 4. fils gauche = 1ere moitiee des triangles et fils droit = 2ieme moitiee des triangles
            // 5. construire les 2 fils
            
        // construire l'englobant des triangles
        Point bmin, bmax;
        bounds(begin, end, bmin, bmax);

        // construire une feuille, s'il ne reste plus qu'un triangle
        if(end - begin < 2)
        {
            int index= int(nodes.size());
            nodes.push_back( {bmin, bmax, -begin, -end} );
            
            /* triche un peu avec la representation du noeud pour stocker une feuille
                au lieu de declarer begin et end dans la structure noeud (ou utiliser 2 structures differentes ?), reutilise left et right, 
            
                mais il faut pourvoir faire la difference entre les 2 cas (noeud / feuille) :
                if(node.right < 0)
                    // feuille
                else
                    // noeud
            
                question: if(node.left < 0) ne fonctionne pas pour detecter une feuille, pourquoi ?
             */
            
            return index;
        }
        
        // axe le plus etire de l'englobant
        int axis= bounds_max(bmin, bmax);
        
        // trier les triangles sur l'axe
        std::sort(triangles.data() + begin, triangles.data() + end, triangle_less(axis));
        
        // construire les fils du noeud
        int left= build_node(begin, (begin+end) / 2);
        int right= build_node((begin+end) / 2, end);
        
        // construire le noeud
        int index= nodes.size();
        nodes.push_back( {bmin, bmax, left, right} );
        return index;
    }
    
    int build_node_centroids( const int begin, const int end )
    {
        // solution 1 : couper l'englobant en 2 sur l'axe le plus etire
            // 1. construire l'englobant
            // 2. trouver l'axe le plus etire
            // 3. repartir les triangles dans chaque moitiee de l'englobant
            // 4. fils gauche = triangles dans la 1ere moitiee de l'englobant et fils droit = triangles dans la 2ieme moitiee de l'englobant
            // 5. construire les 2 fils
            
            // comment repartir les triangles a cheval sur les 2 boites ?
            // solution : un seul point des triangles est teste pour la repartition, le centre de l'englobant, ce qui elimine l'ambiguite.
            
            // consequence indirecte : il n'est plus necessaire de trier les triangles / centres pour les repartir entre les 2 moitiees de l'englobant
            // on aurait pu ecrire :
                // trier les centres des triangles (sur un axe)
                // trouver les n premiers triangles dont les centres sont dans la moitiee gauche de l'englobant
            // mais on peut obtenir le meme resultat directement en parcourant 1 fois l'ensemble de triangles, sans faire le tri
            // c'est ce que fait std::partition()
            // et cette solution est plus rapide : partiton en O(n), tri en O(n log n)
        
        // construire l'englobant des triangles, pour initialiser le noeud / la feuille
        Point bmin, bmax;
        bounds(begin, end, bmin, bmax);
        
        // construire l'englobant des centres des triangles, pour repartir les triangles
        Point cmin, cmax;
        centroid_bounds(begin, end, cmin, cmax);
        
        // construire une feuille, s'il ne reste plus qu'un triangle
        if(end - begin < 2)
        {
            int index= int(nodes.size());
            nodes.push_back( {bmin, bmax, -begin, -end} );
            return index;
        }
        
        // axe le plus etire de l'englobant des centres
        int axis= bounds_max(cmin, cmax);
        
        // repartir les triangles par rapport au milieu de l'englobant des centres
        auto *p= std::partition(triangles.data() + begin, triangles.data() + end, centroid_less1(axis, (cmax(axis)+cmin(axis)) / 2));
        int m= std::distance(triangles.data(), p);
        assert(m != begin);
        assert(m != end);
        
        // construire les fils du noeud
        int left= build_node(begin, m);
        int right= build_node(m, end);
        
        // construire le noeud
        int index= int(nodes.size());
        nodes.push_back( {bmin, bmax, left, right} );
        return index;
    }
    
// utilitaires    
    //! renvoie l'englobant des triangles[begin .. end[
    void bounds( const int begin, const int end, Point& bmin, Point& bmax )
    {
        triangles[begin].bounds(bmin, bmax);
        for(int i= begin +1; i < end; i++)
        {
            Point pmin, pmax;
            triangles[i].bounds(pmin, pmax);
            
            bmin= min(bmin, pmin);
            bmax= max(bmax, pmax);
        }
    }
    
    //! renvoie l'englobant des centres des triangles[begin .. end[
    void centroid_bounds( const int begin, const int end, Point& bmin, Point& bmax )
    {
        Point pm= triangles[begin].centroid();
        bmin= pm;
        bmax= pm;
        for(int i= begin +1; i < end; i++)
        {
            Point pm= triangles[i].centroid();
            bmin= min(bmin, pm);
            bmax= max(bmax, pm);
        }
    }
    
    //! renvoie l'axe le plus etire d'un englobant
    int bounds_max( const Point& pmin, const Point& pmax )
    {
        Vector d= pmax - pmin;
        if(d.x > d.y && d.x > d.z) return 0;
        else if(d.y > d.z) return 1;
        else  return 2;
    }
    
    //! comparaison de 2 triangles, cf build_node, std::sort
    struct triangle_less
    {
        int axis;
        triangle_less( const int _axis ) : axis(_axis) {}
        
        bool operator() ( const Triangle& a, const Triangle& b ) const
        {
            Point apmin, apmax;
            a.bounds(apmin, apmax);
            
            Point bpmin, bpmax;
            b.bounds(bpmin, bpmax);
            
            if(apmin(axis) < bpmin(axis)) return true;
            else if(apmin(axis) == bpmin(axis)) return apmax(axis) < bpmax(axis);
            return false;
        }
    };
    
    //! position d'un triangle, cf build_centroids, std::partition
    struct centroid_less1
    {
        int axis;
        float cut;
        
        centroid_less1( const int _axis, const float _cut ) : axis(_axis), cut(_cut) {}
        
        bool operator() ( const Triangle& a ) const
        {
            Point am= a.centroid();
            return am(axis) < cut;
        }
    };
};


    
float node_area( const Point& pmin, const Point& pmax )
{
    Vector d(pmin, pmax);
    return 2 * d.x*d.y + 2 * d.x*d.z + 2 * d.y*d.z;
}

int main( const int argc, const char **argv )
{
    const char *mesh_filename= "cornell.obj";
    if(argc > 1) mesh_filename= argv[1];
    
    Mesh mesh= read_mesh(mesh_filename);
    assert(mesh.triangle_count());
    printf("triangles %d\n", mesh.triangle_count());
    
    BVH bvh;
    bvh.build(mesh);
    
    printf("root %d, nodes %d, triangles %d\n", bvh.root, int(bvh.nodes.size()), int(bvh.triangles.size()));
    
    // evaluer le cout de l'arbre
    const Node& root= bvh.nodes[bvh.root];
    double root_area= node_area(root.pmin, root.pmax);    // aire de l'englobant de la racine
    
    double cost= 0;
    for(const Node& node : bvh.nodes)
    {
        if(node.right < 0)
        {
            // feuille
            int begin= -node.left;
            int end= -node.right;
            int n= end - begin;
            cost+= node_area(node.pmin, node.pmax) / root_area * n;  // n intersections rayon / triangles par visite
        }
        else
        {
            // noeud interne
            cost+= node_area(node.pmin, node.pmax) / root_area * 1;  // 1 intersection rayon / bbox par visite
        }
    }
    printf("SAH cost %lf\n", cost);
    
    return 0;
}

/* que peut-on ameliorer dans ce code de construction ?
    1. les englobants des triangles sont recalcules a la volee, chaque fois que c'est necessaire, 
    2. les centres des englobants aussi
        il est plus efficace de les calculer une seule fois et de les stocker

    3. les englobants des noeuds sont recalcules plusieurs fois,
        cf bmin, bmax= bounds(begin, end)
        en O(n) pour chaque noeud...
        alors que l'on peut le faire en O(1), temps constant, pour chaque noeud

        il suffit de calculer l'englobant des triangles lorsqu'une feuille est cree, en general il n'y a qu'un triangle,
        et l'englobant d'un noeud est l'union de l'englobants de ses 2 fils...

    4. on manipule beaucoup de boites englobantes, il serait plus lisible de les decrire par une classe plutot que par une paire de points.
    5. il serait aussi plus lisible de representer explicitement un noeud et une feuille : cacher les details de representation (comment faire la difference entre noeud et feuille)

    6. et bien sur, trouver une "meilleure" methode de construction... qui construit (rapidement) un arbre plus rapide a parcourir.

et eventuellement, rendre plus facile le choix de la methode de construction, histoire de simplifier les tests ou les comparaisons entre les differentes versions...
 */
