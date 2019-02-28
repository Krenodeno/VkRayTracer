
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

#define GPU_COMPUTE 0

struct RayGPU {
	alignas(16) vec3 origin;
	alignas(16) vec3 direction;
	float tmax;

	RayGPU(const Point& o, const Point& e) : origin(o), direction(Vector(o, e)), tmax(1.f) {}
	RayGPU(const Point& o, const Vector& d) : origin(o), direction(d), tmax(std::numeric_limits<float>::max()) {}
};

struct HitGPU {
	int id;
	float t, u, v;

	HitGPU() : id(-1), t(0.f), u(0.f), v(0.f) {}

	HitGPU(int id, float t, float u, float v) : id(id), t(t), u(u), v(v) {}

	operator bool() { return id != -1; }
};

struct TriangleGPU {
	alignas(16) vec3 pos;
	alignas(16) vec3 e1;
	alignas(16) vec3 e2;
	int id;

	TriangleGPU(Triangle tri) : pos(tri.p), e1(tri.e1), e2(tri.e2), id(tri.id) {}

	/* Ray - Triangle Intersection */
	HitGPU intersect( const RayGPU &ray, const float htmax ) const
	{
		Vector pvec = cross(Vector(ray.direction), Vector(e2));
		float det = dot(e1, pvec);

		float inv_det = 1 / det;
		Vector tvec(pos, ray.origin);

		float u = dot(tvec, pvec) * inv_det;
		if (u < 0 || u > 1) return HitGPU();

		Vector qvec = cross(tvec, e1);
		float v = dot(ray.direction, qvec) * inv_det;
		if (v < 0 || u + v > 1) return HitGPU();

		float t = dot(e2, qvec) * inv_det;
		if (t > htmax || t < 0) return HitGPU();

		return HitGPU(id, t, u, v);	// p(u, v)= (1 - u - v) * a + u * b + v * c
	}
};

struct BNodeGPU {
	vec3 pmin;
	int next;
	vec3 pmax;
	int skip;

	bool leaf() const { return next < 0; }

	bool node() const { return !leaf(); }
	/** First triangle of the leaf */
	int begin() const { assert(leaf()); return ( -next >> 4 ); }
	/** Last triangle of the leaf */
	int end() const { assert(leaf()); return ( (-next >> 4) + (-next & 0xF) ); }

	/* Ray - Box Intersection */
	NodeHit intersect( const RayGPU& ray, const Vector& invd, const float htmax ) const
	{
		vec3 rmin = pmin;
		vec3 rmax = pmax;
		if(ray.direction.x < 0) std::swap(rmin.x, rmax.x);
		if(ray.direction.y < 0) std::swap(rmin.y, rmax.y);
		if(ray.direction.z < 0) std::swap(rmin.z, rmax.z);
		Vector dmin = (Point(rmin) - Point(ray.origin)) * invd;
		Vector dmax = (Point(rmax) - Point(ray.origin)) * invd;

		float tmin = std::max(dmin.z, std::max(dmin.y, std::max(dmin.x, 0.f)));
		float tmax = std::min(dmax.z, std::min(dmax.y, std::min(dmax.x, htmax)));
		return NodeHit(tmin, tmax);
	}

	BNodeGPU() : skip(-1) {}

};


std::vector<BNodeGPU> transform(const BVH& bvh) {

	assert(!bvh.nodes[bvh.root].leaf());

	std::vector<BNodeGPU> remap;
	remap.reserve(bvh.nodes.size());
	std::vector<int> stackG;
	int current = 0;	// current node in the remapped bvh

	int id = bvh.root;	// current node in the original bvh
	bool cont = true;

	std::vector<int> stack;
	stack.reserve(bvh.nodes.size());

	// First step : root
	{
		const Node& node = bvh.nodes[bvh.root];
		BNodeGPU root;
		root.pmin = vec3(node.bounds.pmin.x, node.bounds.pmin.y, node.bounds.pmin.z);
		root.pmax = vec3(node.bounds.pmax.x, node.bounds.pmax.y, node.bounds.pmax.z);
		root.skip = -1;
		root.next = current + 1;	// doesn't work is root is a leaf

		id = node.left;
		stackG.push_back(current++);
		stack.push_back(bvh.root);
		remap.push_back(std::move(root));
	}

	while (cont) {
		const Node& node = bvh.nodes[id];
		{
			BNodeGPU res;
			res.pmin = vec3(node.bounds.pmin.x, node.bounds.pmin.y, node.bounds.pmin.z);
			res.pmax = vec3(node.bounds.pmax.x, node.bounds.pmax.y, node.bounds.pmax.z);

			if (node.leaf()) {
				res.next = -( (node.begin() << 4) | node.n() );	// Triangle indices (the offset must be less than 16)
				if (bvh.nodes[stack.back()].left == id)
					res.skip = remap.size() + 1;	// brother of a leaf is next index
				assert(res.leaf());
			}
			else {
				res.next = current + 1;
				assert(res.node());
			}
			remap.push_back(res);
		}

		// tree walking
		if (node.leaf()) {
			int parent = id;
			// go back in the tree until there is a right child node not visited
			while (bvh.nodes[stack.back()].right == parent) {

				// write skip as it is now knowable
				remap[current].skip = remap.size();

				// end of the loop : everything have been visited
				if (stack.empty())
					break;
				// or not : continue to go up
				parent = stack.back();
				stack.pop_back();
				current = stackG.back();
				stackG.pop_back();
			}
			id = bvh.nodes[stack.back()].right;
			current = remap.size();
		}
		else {
			stack.push_back(id);
			id = node.left;
			stackG.push_back(current++);
		}

		// end of the loop
		if (stack.empty())
			cont = false;
	}

	// little cleanup : set node.skip to zero if it is set to bvh.size()
	for (auto node : remap) {
		if (node.skip == remap.size())
			node.skip = -1;
	}

	return remap;
}


HitGPU intersect(const std::vector<BNodeGPU>& pbvh, const std::vector<TriangleGPU>& tris, RayGPU ray) {
	assert(!pbvh.empty());

	Vector invd = Vector(1.f / ray.direction.x, 1.f / ray.direction.y, 1.f / ray.direction.z);

	int id = 0; // 0 is root, but in any node, it means no next node
	do {
		const BNodeGPU& node = pbvh[id];
		if (node.leaf()) {
			for (int i = node.begin(); i != node.end(); i++) {
				HitGPU hit = tris[i].intersect(ray, ray.tmax);
				if (hit)
					return hit;	// intersection is found
			}
			// No intersection
			id = node.skip;
		}
		else {
			if (node.intersect(ray, invd, ray.tmax))
				id = node.next;
			else
				id = node.skip;
		}
	} while (id > -1);
	// No intersection
	return HitGPU();
}


Vector normal( const HitGPU& hit, const TriangleData& triangle )
{
	return normalize((1 - hit.u - hit.v) * Vector(triangle.na) + hit.u * Vector(triangle.nb) + hit.v * Vector(triangle.nc));
}


Point point( const Hit& hit, const TriangleData& triangle )
{
	return (1 - hit.u - hit.v) * Point(triangle.a) + hit.u * Point(triangle.b) + hit.v * Point(triangle.c);
}


Point point( const HitGPU& hit, const RayGPU& ray )
{
	return Point(ray.origin) + hit.t * Vector(ray.direction);
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


void fillVertexBuffer(ComputeApp& compute, uint32_t bufIndex, const Mesh& mesh) {
	std::vector<vec3> positions;
	for (int i = 0; i < mesh.index_count(); i++) {
		positions.push_back(mesh.positions()[mesh.indices()[i]]);
	}
	compute.fillBuffer(bufIndex, positions.data(), positions.size() * sizeof(positions[0]));
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
	ComputeApp compute();

	// Transformer l'arbre pour le parcourir sur GPU
	// l'id 0 est la racine
	std::vector<BNodeGPU> pbvh = transform(bvh);

	std::vector<TriangleGPU> triangles;
	for (const auto& tri : bvh.triangles) {
		triangles.emplace_back(tri);
	}

	assert(!pbvh.empty());
	assert(!triangles.empty());


#if GPU_COMPUTE
	// La taille des buffers est en octets
	// Copie sur le GPU
	// Géométrie : Vertices
	uint32_t vertexBufferSize = mesh.index_count() * sizeof(mesh.positions()[0]);
	compute.addBuffer(vertexBufferSize);
	fillVertexBuffer(compute, 0);
	// BVH
	uint32_t bvhBufferSize = pbvh.size() * sizeof(pbvh[0]);
	compute.addBuffer(bvhBufferSize);
	compute.fillBuffer(1, pbvh.data(), bvhBufferSize);
	// Rayons
	compute.addBuffer();
	// Hits
	compute.addBuffer();
#endif

	auto cpu_start = std::chrono::high_resolution_clock::now();

	std::cout << "Generating Fibonacci Distribution\n";
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

	std::cout << "Generating firsts Rays directions\n";
	std::vector<RayGPU> rays;
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

			rays.push_back(RayGPU(o, e));
		}
	}

#if GPU_COMPUTE
	fillBuffer();
#endif

	std::cout << "Intersecting...\n";
	// Calculer les intersections
	// parcourir tous les pixels de l'image
	// en parallele avec openMP, un thread par bloc de 16 lignes
#pragma omp parallel for schedule(dynamic, 16)
	for (int py = 0; py < image.height(); py++) {
		const int dirs = directionCount;
		for (int px = 0; px < image.width(); px++) {

			int index = px + py * image.width();
			RayGPU ray = rays[index];
			if (HitGPU hit = intersect(pbvh, triangles, ray))
			{
				// recupere les donnees sur l'intersection
				TriangleData triangle = mesh.triangle(hit.id);
				Point p = point(hit, ray);			// point d'intersection
				Vector pn = normal(hit, triangle);	// normale interpolee du triangle au point d'intersection
				if(dot(pn, ray.direction) > 0)				// retourne la normale vers l'origine du rayon
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
				float cos_theta = dot(normalize(pn), normalize(-ray.direction));
				Color pColor = std::max(0.f, cos_theta) * mesh.mesh_material(mesh.materials()[hit.id]).diffuse;
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
