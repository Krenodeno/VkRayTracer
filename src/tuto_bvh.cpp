
#include <cfloat>
#include <vector>

#include "vec.h"
#include "mesh.h"
#include "wavefront.h"

#include "ray.h"
#include "bvh.h"
#include "centroid_builder.h"


int main( const int argc, const char **argv )
{
	const char *mesh_filename = "resources/cornell.obj";
	if(argc > 1) mesh_filename = argv[1];

	Mesh mesh = read_mesh(mesh_filename);
	assert(mesh.triangle_count());
	printf("triangles %d\n", mesh.triangle_count());

	CentroidBuilder builder;

	BVH bvh(mesh, builder);

	printf("root %d, nodes %d, triangles %d\n", bvh.root, int(bvh.nodes.size()), int(bvh.triangles.size()));

	// evaluer le cout de l'arbre
	double root_area = bvh.nodes[bvh.root].bounds.area();	// aire de l'englobant de la racine
	double cost = 0;
	for(const Node& node : bvh.nodes)
	{
		if(node.leaf())
		{
			// feuille
			int begin = node.begin();
			int end = node.end();
			int n = end - begin;
			cost += node.bounds.area() / root_area * n;	// n intersections rayon / triangles par visite
		}
		else
		{
			// noeud interne
			cost += node.bounds.area() / root_area * 1;	// 1 intersection rayon / bbox par visite
		}
	}
	printf("SAH cost %lf\n", cost);

	return 0;
}
