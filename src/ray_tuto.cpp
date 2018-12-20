
#include <cfloat>
#include <random>
#include <chrono>

#include "vec.h"
#include "mesh.h"
#include "wavefront.h"
#include "orbiter.h"

#include "ray.h"
#include "bvh.h"
#include "centroid_builder.h"

#include "image.h"
#include "image_io.h"
#include "image_hdr.h"


Vector normal( const Hit& hit, const TriangleData& triangle )
{
	return normalize((1 - hit.u - hit.v) * Vector(triangle.na) + hit.u * Vector(triangle.nb) + hit.v * Vector(triangle.nc));
}

Point point( const Hit& hit, const TriangleData& triangle )
{
	return (1 - hit.u - hit.v) * Point(triangle.a) + hit.u * Point(triangle.b) + hit.v * Point(triangle.c);
}


Point point( const Hit& hit, const Ray& ray )
{
	return ray.o + hit.t * ray.d;
}


int main( const int argc, const char **argv )
{
	const char *mesh_filename= "resources/cornell.obj";
	const char *orbiter_filename= "resources/cornell_orbiter.txt";

	if(argc > 1) mesh_filename= argv[1];
	if(argc > 2) orbiter_filename= argv[2];

	printf("%s: '%s' '%s'\n", argv[0], mesh_filename, orbiter_filename);

	// creer l'image resultat
	Image image(1024, 640);

	// charger un objet
	Mesh mesh= read_mesh(mesh_filename);
	if(mesh.triangle_count() == 0)
		// erreur de chargement, pas de triangles
		return 1;

	// construire le bvh ou recuperer l'ensemble de triangles du mesh...
	CentroidBuilder builder;
	BVH bvh(mesh, builder);

	if( 1 )
	{
		// evaluer le cout de l'arbre
		double cost= 0;
		double root_area= bvh.nodes[bvh.root].bounds.area();	// aire de l'englobant de la racine
		for(const Node& node : bvh.nodes)
		{
			if(node.leaf())
				// feuille
				cost+= node.bounds.area() / root_area * (node.end() - node.begin());  // n intersections rayon / triangles par visite
			else
				// noeud interne
				cost+= node.bounds.area() / root_area * 1;  // 1 intersection rayon / bbox par visite
		}
		printf("SAH cost %lf\n", cost);
	}

	// charger la camera
	Orbiter camera;
	if(camera.read_orbiter(orbiter_filename))
		// erreur, pas de camera
		return 1;

	// recupere les transformations view, projection et viewport pour generer les rayons
	float fov = 45.f;
	Transform model= Identity();
	Transform view= camera.view();
	Transform projection= camera.projection(image.width(), image.height(), fov);

	auto cpu_start= std::chrono::high_resolution_clock::now();

	// parcourir tous les pixels de l'image
	// en parallele avec openMP, un thread par bloc de 16 lignes
#pragma omp parallel for schedule(dynamic, 16)
	for(int py= 0; py < image.height(); py++)
	{
		for(int px= 0; px < image.width(); px++)
		{
			// generer le rayon pour le pixel (x, y)
			float x= px + .5f;		  // centre du pixel
			float y= py + .5f;

			Point o = camera.position();	// origine

			Point d1;
			Vector dx1, dy1;
			camera.frame(image.width(), image.height(), 1, fov, d1, dx1, dy1);
			Point e = d1 + x*dx1 + y*dy1;	// extremite

			// calculer les intersections
			Ray ray(o, e);
			if(Hit hit= bvh.intersect(ray))
			{
				// recupere les donnees sur l'intersection
				TriangleData triangle= mesh.triangle(hit.triangle_id);
				Point p= point(hit, ray);			// point d'intersection
				Vector pn= normal(hit, triangle);	// normale interpolee du triangle au point d'intersection
				if(dot(pn, ray.d) > 0)				// retourne la normale vers l'origine du rayon
					pn= -pn;

				// couleur du pixel
				Color color= mesh.mesh_material(mesh.materials()[hit.triangle_id]).diffuse;
				image(px, py)= Color(color, 1);
			}
		}
	}

	auto cpu_stop= std::chrono::high_resolution_clock::now();
	int cpu_time= std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
	printf("cpu  %ds %03dms\n", int(cpu_time / 1000), int(cpu_time % 1000));

	// enregistrer l'image resultat
	write_image(image, "render.png");
	write_image_hdr(image, "render.hdr");

	return 0;
}
