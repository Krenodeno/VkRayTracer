
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

#include "ComputeApp.hpp"

struct BBox {
	vec3 min;
	vec3 max;
};

struct BNode {
	vec3 pmin;
	int left;
	vec3 pmax;
	int right;
};

struct Ray {
	vec3 origin;
	vec3 direction;
	float tmax;
};

struct Hit {
	int id;
	float t, u, v;
}


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


// b1, b2, n sont 3 axes orthonormes.
void repere( const Vector &n, Vector &b1, Vector &b2 )
{
	float sign = std::copysign(1.0f, n.z);
	const float a = -1.0f / (sign + n.z);
	const float b = n.x * n.y * a;
	b1 = Vector(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
	b2 = Vector(b, sign + n.y * n.y * a, -n.y);
}


int main( const int argc, const char **argv )
{

	int directionCount = 16;
	const char *mesh_filename = "resources/cornell.obj";
	const char *orbiter_filename = "resources/cornell_orbiter.txt";

	if (argc > 1) mesh_filename = argv[1];
	if (argc > 2) orbiter_filename = argv[2];
	if (argc > 3) directionCount = std::stoi(argv[3]);

	printf("%s: '%s' '%s'\n", argv[0], mesh_filename, orbiter_filename);
	printf("Nombre de directions par rebond : %i\n", directionCount);

	// creer l'image resultat
	Image image(1024, 640);

	// charger un objet
	Mesh mesh = read_mesh(mesh_filename);
	if (mesh.triangle_count() == 0)
		// erreur de chargement, pas de triangles
		return 1;

	// construire le bvh ou recuperer l'ensemble de triangles du mesh...
	CentroidBuilder builder;
	BVH bvh(mesh, builder);

	if ( 1 )
	{
		// evaluer le cout de l'arbre
		double cost = 0;
		double root_area = bvh.nodes[bvh.root].bounds.area();	// aire de l'englobant de la racine
		for (const Node& node : bvh.nodes)
		{
			if (node.leaf())
				// feuille
				cost += node.bounds.area() / root_area * (node.end() - node.begin());  // n intersections rayon / triangles par visite
			else
				// noeud interne
				cost += node.bounds.area() / root_area * 1;  // 1 intersection rayon / bbox par visite
		}
		printf("SAH cost %lf\n", cost);
	}

	// charger la camera
	Orbiter camera;
	if (camera.read_orbiter(orbiter_filename))
		// erreur, pas de camera
		return 1;

	// recupere les transformations view, projection et viewport pour generer les rayons
	float fov = 45.f;
	Transform model = Identity();
	Transform view = camera.view();
	Transform projection = camera.projection(image.width(), image.height(), fov);

	// Créer le pipeline compute
	ComputeApp compute(image.width(), image.height());

	// La taille des buffers est en octets
	// Géométrie : Vertices & Indices
	//compute.addBuffer();
	// BVH
	//compute.addBuffer();
	// Rayons
	//compute.addBuffer();
	// Hits
	//compute.addBuffer();

	auto cpu_start = std::chrono::high_resolution_clock::now();

	// Générer la spirale de Fibonacci
	std::vector<Vector> fiboDistribution;
	// phi = sqrt(5) + 1 / 2
	const float phi = 1.61803398875f;
	for (int i = 0; i < directionCount; i++) {
		float cosThetaI = 1.f - (2.f * i+1) / (float)(2.f * directionCount);
		float phiI = 2.f * M_PI * ((i/phi) - std::floor(i/phi));

		float sinThetaI = std::sqrt(1 - (cosThetaI * cosThetaI));
		float cosPhiI = std::cos(phiI);
		float sinPhiI = std::sin(phiI);
		Vector dir(cosPhiI * sinThetaI, sinPhiI * cosThetaI, cosThetaI);
		fiboDistribution.push_back(dir);
	}

	std::vector<Ray> rays;
	rays.reserve(image.height() * image.width());
	// Générer les rayons
	for (int py = 0; py < image.height(); py++) {
		for (int px = 0; px < image.width(); px++) {
			// generer le rayon pour le pixel (x, y)
			float x = px + .5f;		// centre du pixel
			float y = py + .5f;

			Point o = camera.position();	// origine

			Point d1;
			Vector dx1, dy1;
			camera.frame(image.width(), image.height(), 1, fov, d1, dx1, dy1);
			Point e = d1 + x*dx1 + y*dy1;	// extremite

			rays.push_back(Ray(o, e));
		}
	}

	// Calculer les intersections
	// parcourir tous les pixels de l'image
	// en parallele avec openMP, un thread par bloc de 16 lignes
#pragma omp parallel for schedule(dynamic, 16)
	for (int py = 0; py < image.height(); py++) {
		const int dirs = directionCount;
		for (int px = 0; px < image.width(); px++) {

			int index = px + py * image.width();
			Ray ray = rays[index];
			if (Hit hit = bvh.intersect(ray))
			{
				// recupere les donnees sur l'intersection
				TriangleData triangle = mesh.triangle(hit.triangle_id);
				Point p = point(hit, ray);			// point d'intersection
				Vector pn = normal(hit, triangle);	// normale interpolee du triangle au point d'intersection
				if(dot(pn, ray.d) > 0)				// retourne la normale vers l'origine du rayon
					pn = -pn;

				float sum = 0.f;

				// Relancer des rayon grâce à la spirale de Fibonacci
				for (int i = 0; i < directionCount; i++) {
					Vector b1, b2;
					repere(pn, b1, b2);

					//Vector w(pn * fiboDistribution[i]);
					Vector w = fiboDistribution[i].x * b1 + fiboDistribution[i].y * b2 + fiboDistribution[i].z * pn;

					Ray shadow(p + pn * .001f, p + w * 10.f);

					if (bvh.visible(shadow)) {
						sum += 1.f;
					}
				}

				// couleur du pixel
				float cos_theta = dot(normalize(pn), normalize(-ray.d));
				Color pColor = std::max(0.f, cos_theta) * mesh.mesh_material(mesh.materials()[hit.triangle_id]).diffuse;
				Color color = 1.f/(float)directionCount * sum * pColor;

				image(px, py) = Color(color, 1);
			}
		}
	}

	auto cpu_stop = std::chrono::high_resolution_clock::now();
	int cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(cpu_stop - cpu_start).count();
	printf("cpu  %ds %03dms\n", int(cpu_time / 1000), int(cpu_time % 1000));

	// enregistrer l'image resultat
	write_image(image, "render1.png");
	write_image_hdr(image, "render1.hdr");

	return 0;
}
